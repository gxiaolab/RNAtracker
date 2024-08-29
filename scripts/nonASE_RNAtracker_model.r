

# EM algorithm functions
#' Expectation Step of the EM Algorithm
#'
#' Calculate the posterior probabilities (soft labels) that each component
#' has to each data point.

deal_with_zeroes = function(sum.of.comps, comp0.post, comp1.post, comp2.post, comp0_df, comp1_df, comp2_df){
    sum.of.comps.ln <- log(sum.of.comps, base = exp(1))

    # 0 exists in sum.of.comps
    if (sum(sum.of.comps == 0) > 0){
        small_prob_gene_id <- which(sum.of.comps == 0)
        small_post_mat <- matrix(0, nrow=length(small_prob_gene_id), ncol=3)
        for (cur_small_prob_gene_id in small_prob_gene_id){
            loggp_vec_prec <- c(exp(Rmpfr::mpfr(comp0_df$logP[cur_small_prob_gene_id],prec=40)), exp(Rmpfr::mpfr(comp1_df$logP[cur_small_prob_gene_id],prec=40)), exp(Rmpfr::mpfr(comp2_df$logP[cur_small_prob_gene_id],prec=40)))
            loggp_vec_prec <- loggp_vec_prec * pi
            small_post_mat[which(small_prob_gene_id==cur_small_prob_gene_id),] <- as.numeric(loggp_vec_prec / sum(loggp_vec_prec))
        }
    comp0.post[small_prob_gene_id] <- small_post_mat[,1]
    comp1.post[small_prob_gene_id] <- small_post_mat[,2]
    comp2.post[small_prob_gene_id] <- small_post_mat[,3]
    sum.of.comps.ln[small_prob_gene_id] <- 0
    }

    sum.of.comps.ln.sum <- sum(sum.of.comps.ln)
    posterior.df <- data.frame(gene_id=comp0_df[,1], comp0.post, comp1.post, comp2.post)
     list("loglik" = sum.of.comps.ln.sum,
       "posterior.df" = posterior.df)
}



#'
#' @param x a List of vectors containing observed reference allele count across timepoints
#' @param n a List of vectors containing observed total count across timepoints
#' @param gene_ind_list a List containing gene ind (length of c(x)) across timepoints
#' @param pi Vector containing the multicategorical distribution mean parameter of each component (length 7)
#' @param a0 Scalar for balanced Beta distribution hyparameter alpha
#' @param b0 Scalar for balanced Beta distribution hyparameter beta
#' @param a1 Scalar for balanced Beta distribution hyparameter alpha
#' @param b1 Scalar for balanced Beta distribution hyparameter beta
#' @return Named list containing the loglik and posterior.df  
e_step <- function(x, n, gene_ind_list, pi, a0, b0, a1, b1){
    gene_ind_vec <- unlist(gene_ind_list)
    # rho = 1 / (1 + a + b); mu = a / (a + b)
    mu0 <- a0 / (a0 + b0)
    rho0 <- 1 / (1 + a0 + b0)
    mu1 <- a1 / (a1 + b1)
    rho1 <- 1 / (1 + a1 + b1)

    # Separate timepoints
    x_time2 <- x[[1]]
    x_time6 <- x[[2]]
    n_time2 <- n[[1]]
    n_time6 <- n[[2]]

    # component complete prob density value P(X,Y = k|theta)
    comp0.dens <- c(dbetabinom(x_time2, n_time2, prob=mu0, rho=rho0), dbetabinom(x_time6, n_time6, prob=mu0, rho=rho0)) # vector of length x
    comp1.dens <- c(dbetabinom(x_time2, n_time2, prob=mu1, rho=rho1), dbetabinom(x_time6, n_time6, prob=mu1, rho=rho1))
    comp2.dens <- c(dbetabinom(x_time2, n_time2, prob=mu0, rho=rho0), dbetabinom(x_time6, n_time6, prob=mu1, rho=rho1))
  
    comp0_df <- data.frame(dens = comp0.dens, gene_id = gene_ind_vec) %>% group_by(gene_id) %>% summarise(product = prod(dens), logP = (sum(log(dens))))
    comp1_df <- data.frame(dens = comp1.dens, gene_id = gene_ind_vec) %>% group_by(gene_id) %>% summarise(product = prod(dens), logP = (sum(log(dens))))
    comp2_df <- data.frame(dens = comp2.dens, gene_id = gene_ind_vec) %>% group_by(gene_id) %>% summarise(product = prod(dens), logP = (sum(log(dens))))

    # Posterior mean Y_g | X_g, theta
    sum.of.comps <- pi[1] * comp0_df$product + pi[2] * comp1_df$product + pi[3] * comp2_df$product 
    comp0.post <- pi[1] * comp0_df$product / sum.of.comps
    comp1.post <- pi[2] * comp1_df$product / sum.of.comps
    comp2.post <- pi[3] * comp2_df$product / sum.of.comps

 
    out = deal_with_zeroes(sum.of.comps, comp0.post, comp1.post, comp2.post, comp0_df, comp1_df, comp2_df)
    sum.of.comps.ln.sum = out[["loglik"]]
    posterior.df = out[["posterior.df"]]



  list("loglik" = sum.of.comps.ln.sum,
       "posterior.df" = posterior.df)
       }




#' Q function
#' @param x Input data.
#' @param posterior.df Posterior probability data.frame.
# posterior.df <- e.step$posterior.df
q_function <- function(x, n, a0, b0, a1, b1, gene_ind_list, posterior.df){
    # rho = 1 / (1 + a + b); mu = a / (a + b)
    mu0 <- a0 / (a0 + b0)
    rho0 <- 1 / (1 + a0 + b0)
    mu1 <- a1 / (a1 + b1)
    rho1 <- 1 / (1 + a1 + b1)

    # Separate timepoints
    x_time2 <- x[[1]]
    x_time6 <- x[[2]]
    n_time2 <- n[[1]]
    n_time6 <- n[[2]]

    # Expand posterior.df according to gene repeat times
    gene_to_complete_id <- lapply(gene_ind_list, function(x){
        return(match(x, posterior.df$gene_id))
    })
    complete.posterior.df <- posterior.df[unlist(gene_to_complete_id),-1]

    # component prob density value P(X,Y = k|theta)
    comp0.log <- c(dbetabinom(x_time2, n_time2, prob=mu0, rho=rho0, log=TRUE), dbetabinom(x_time6, n_time6, prob=mu0, rho=rho0, log=TRUE))
    comp1.log <- c(dbetabinom(x_time2, n_time2, prob=mu1, rho=rho1, log=TRUE), dbetabinom(x_time6, n_time6, prob=mu1, rho=rho1, log=TRUE))
    comp2.log <- c(dbetabinom(x_time2, n_time2, prob=mu0, rho=rho0, log=TRUE), dbetabinom(x_time6, n_time6, prob=mu1, rho=rho1, log=TRUE))
  
    # Wrong here
    comp0.log.sum <- sum(complete.posterior.df[,1] * comp0.log)
    comp1.log.sum <- sum(complete.posterior.df[,2] * comp1.log)
    comp2.log.sum <- sum(complete.posterior.df[,3] * comp2.log)

    return(Q= comp0.log.sum + comp1.log.sum + comp2.log.sum )

}






#' Maximization Step of the EM Algorithm
#'
#' Update the Component Parameters
#'
#' @param x Input data.
#' @param posterior.df Posterior probability data.frame.
#' @return Named list containing the mean (mu), variance (var), and mixing
#'   weights (alpha) for each component.
m_step <- function(x, n, gene_ind_list, posterior.df, previous_theta) {
    # Expand posterior.df according to gene repeat times
    gene_to_complete_id <- lapply(gene_ind_list, function(x){
        return(match(x, posterior.df$gene_id))
    })
    complete.posterior.df <- posterior.df[unlist(gene_to_complete_id),-1]
    # Analytical solution to Pi
    pi_est <- apply(posterior.df[-1], MARGIN = 2, function(x){
        sum(x) / length(x)
    })

    # Numerical solution to a0,a1,b1
    optim_fun <- function(theta){
        a0 <- theta[1]
        a1 <- theta[2]
        return(-1 * q_function(x, n, a0, a0, a1, a1, gene_ind_list, posterior.df))
    }
    cur_theta_est <- tryCatch({
        # "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent"
        # suppressMessages(coef(mle(mle_opt, start=list(a0=previous_theta[1], a1=previous_theta[2], b1=previous_theta[3])), method=c("L-BFGS-B")))
        suppressMessages(c(optim(par=c(a0=previous_theta[1], a1=previous_theta[2]), fn=optim_fun, method=c("Nelder-Mead"))$par, use.names = FALSE))
    }, error =  function(e) {
        print("M-step optimization failed, NA stored instead.")
        NA
    })
    return(list(pi_est = pi_est, theta_est=cur_theta_est))
}
