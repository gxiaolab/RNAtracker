rm(list = ls())
library(VGAM)
library(dplyr)
library(data.table)
source('preprocess_0h_RNAtracker_model.r')

args = commandArgs(trailingOnly=TRUE)


CELL_LINE = args[1] # e.g. GM12878
INPUT_DIR = args[2] # ../data
RESULTS_DIR = args[3] # ../results


estimates_fn = '../data/preprocess_0h.initial_values.txt'
estimated_initial_values = fread(estimates_fn)


A0_INIT = estimated_initial_values$a0[estimated_initial_values$cell_line == CELL_LINE]
B0_INIT = estimated_initial_values$a0[estimated_initial_values$cell_line == CELL_LINE] # intentionally use the same estimate as a0 for b0
A1_INIT = estimated_initial_values$a1[estimated_initial_values$cell_line == CELL_LINE]
B1_INIT = estimated_initial_values$a1[estimated_initial_values$cell_line == CELL_LINE]  # intentionally use the same estimate as a1 for b1


A0_INIT = as.numeric(A0_INIT)
B0_INIT = as.numeric(B0_INIT)
A1_INIT = as.numeric(A1_INIT)
B1_INIT = as.numeric(B1_INIT)

pi_init = c(1/2, 1/2)

#dir.create(RESULTS_DIR, recursive = TRUE)

input_fn = paste0(INPUT_DIR, '/', CELL_LINE, '.snv_info.txt.gz')

all_data = fread(input_fn)

# only testable genes have BEAPR p-values 
pass_genes = all_data %>% filter(!is.na(p_adj_0h)) %>% select(id, gene_name) %>% unique() %>% count(gene_name) 
all_data = filter(all_data, gene_name %in% pass_genes$gene_name)


number_snvs = length(unique(all_data$id))
print(paste("we are testing", length(unique(all_data$id)), 'snvs'))

# Define the Log-likelihood function
time0_ref_allele_count <- all_data$MODEL_COUNT_0h
time0_total_count <- all_data$total_counts_0h



gene_ind_list <- list(all_data$gene_name, all_data$gene_name)
n_gene <- length(unique(unlist(gene_ind_list)))


# Observed data
x <- list(time0_ref_allele_count)
n <- list(time0_total_count)

# EM algorithm
tic <- proc.time()
e.step <- e_step(x, n, gene_ind_list, pi_init, A0_INIT, B0_INIT, A1_INIT, B1_INIT)
m.step <- m_step(x, n, gene_ind_list, e.step[["posterior.df"]], c(A0_INIT, A1_INIT)) 

cur.loglik <- e.step[["loglik"]]
loglik.vector <- e.step[["loglik"]]
for (i in 1:1000) {
    # Repeat E and M steps till convergence
    previous_theta <- m.step$theta_est
    e.step <- e_step(x, n, gene_ind_list, m.step$pi_est, previous_theta[1], previous_theta[1], previous_theta[2], previous_theta[2]) 
    m.step <- m_step(x, n, gene_ind_list, e.step[["posterior.df"]], previous_theta)
    loglik.vector <- c(loglik.vector, e.step[["loglik"]])
    loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
    if (length(m.step$theta_est)==1){
        print(sprintf("Iteration %s M-step failed.", i))
        break
    }
    if(loglik.diff < 1e-6) {
        print(paste('MODEL CONVERGED: finished classifying', n_gene, 'genes'))
        break
    } else {
    cur.loglik <- e.step[["loglik"]]
    }
    lapsed_time <- proc.time() - tic
    print(sprintf("[Iteration %s]: %s seconds used. Log-likelihood: %s.", i, lapsed_time[1], e.step[["loglik"]]))
    print(m.step$pi_est)
}

# Posterior prob for gene states
final_step_post_raw <- e.step$posterior.df
final_step_param <- m.step


summarize_results = function(final_step_post){
    predicted_state = apply(final_step_post[,-1], MARGIN=1, which.max)
    max_probability = apply(final_step_post[,-1], MARGIN=1, max)
    # subtract 1 to account for indexing
    predicted_state = predicted_state - 1
    final_step_post$predicted_state = predicted_state
    final_step_post$max_probability = max_probability
    total_genes = nrow(final_step_post)
    summary_results = final_step_post %>% group_by(predicted_state) %>% count()
    summary_results = summary_results %>% rename(number_genes = n) %>% mutate(prop_genes = number_genes/total_genes)
    return(list(final_step_post, summary_results))
}

final_step_post = summarize_results(final_step_post_raw)[[1]]
summary_results = summarize_results(final_step_post_raw)[[2]]


print('CLASSIFICATION RESULTS:')
summary_results %>% print(n=100)

print('FINAL PARAMETERS')
final_step_param


initial_parameters = c(A0_INIT, A1_INIT, B1_INIT, pi_init)
save(final_step_post, final_step_param, summary_results, initial_parameters, e.step, file= sprintf("%s/%s.preprocess_0h.Rdata", RESULTS_DIR, CELL_LINE))


