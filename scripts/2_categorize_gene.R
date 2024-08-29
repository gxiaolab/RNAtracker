rm(list = ls())
library(VGAM)
library(dplyr)
library(data.table)


args = commandArgs(trailingOnly=TRUE)
CELL_LINE = args[1] # e.g. GM12878
INPUT_DIR = args[2] # ../data
RESULTS_DIR = args[3] # ../results
MAX_PROB_CUTOFF = 0.95
STATE_QUANT_CUTOFF = 0.333333333333333


# read in initial values

estimates_fn = '../data/gene_categorization_initial_values.txt'
estimated_initial_values = fread(estimates_fn)


A0_INIT = estimated_initial_values$a0[estimated_initial_values$cell_line == CELL_LINE]
B0_INIT = estimated_initial_values$a0[estimated_initial_values$cell_line == CELL_LINE] # intentionally use the same estimate as a0 for b0
A1_INIT = estimated_initial_values$a1[estimated_initial_values$cell_line == CELL_LINE]
B1_INIT = estimated_initial_values$a1[estimated_initial_values$cell_line == CELL_LINE]  # intentionally use the same estimate as a1 for b1


A0_INIT = as.numeric(A0_INIT)
B0_INIT = as.numeric(B0_INIT)
A1_INIT = as.numeric(A1_INIT)
B1_INIT = as.numeric(B1_INIT)


### load preprocessed data 

in_fn = paste0(RESULTS_DIR, '/', CELL_LINE, '.preprocess_0h.Rdata')
load(in_fn)
all_df = final_step_post

cutoff_df = all_df %>% group_by(predicted_state) %>% summarise(quantile = scales::percent(c(STATE_QUANT_CUTOFF )),quant_cutoff = quantile(max_probability, c(STATE_QUANT_CUTOFF )))
all_df = left_join(all_df, cutoff_df, by = 'predicted_state')

zero_hr_assignments = all_df %>% filter(max_probability >= MAX_PROB_CUTOFF | max_probability >= quant_cutoff)


print(paste0('before filtering with MAX_PROB_CUTOFF=', MAX_PROB_CUTOFF,  ', STATE_QUANT_CUTOFF=', STATE_QUANT_CUTOFF, ': ', nrow(all_df), ' genes'))
print(paste0('after filtering: ', nrow(zero_hr_assignments), ' genes (', round(nrow(zero_hr_assignments)/nrow(all_df),3), ')'))

filter_label = paste0('MAX_PROB_CUTOFF_', MAX_PROB_CUTOFF, '_STATE_QUANT_CUTOFF_', STATE_QUANT_CUTOFF)
print(table(zero_hr_assignments$predicted_state))
zero_hr_assignments$predicted_state = paste0('state', zero_hr_assignments$predicted_state)




make_predictions = function(WHICH_MODEL){
    if( WHICH_MODEL == 'nonASE_0h'){
        source('nonASE_RNAtracker_model.r')
        pi_init = c(1/3, 1/3, 1/3)
        GENES_TO_TEST = zero_hr_assignments %>% filter(predicted_state == 'state0')
        GENES_TO_TEST = GENES_TO_TEST$gene_id
        #print(GENES_TO_TEST)
    } else if (WHICH_MODEL == 'ASE_0h'){
        source('ASE_RNAtracker_model.r')
        pi_init = c(1/4, 1/4, 1/4, 1/4)
        GENES_TO_TEST = zero_hr_assignments %>% filter(predicted_state == 'state1')
        GENES_TO_TEST = GENES_TO_TEST$gene_id
        #print(GENES_TO_TEST)
    }

    output_fn = paste0(CELL_LINE, '_', WHICH_MODEL, '_model')


    input_fn = paste0(INPUT_DIR, '/', CELL_LINE, '.snv_info.txt.gz')


    all_data = fread(input_fn)


    all_data = all_data %>% filter(gene_name %in% GENES_TO_TEST)

    # in addition to being testable at 0h, gene must also be testable at 2h or 6h
    two_hr_testable = all_data %>% filter(!is.na(p_adj_2h)) # only testable genes have BEAPR p-values
    six_hr_testable = all_data %>% filter(!is.na(p_adj_6h)) # only testable genes have BEAPR p-values
    keep_snvs = unique(c(two_hr_testable$id, six_hr_testable$id))
    all_data = all_data  %>% filter(id %in% keep_snvs)


    print(paste("after filtering for snvs that are in 0h classified genes and testable at 2h or 6h, we have ", length(unique(all_data$id)), 'snvs'))

    if(length(unique(all_data$id))== 0){
        final_step_post = data.frame(cell_line = CELL_LINE, predicted_state = NA)
        save(final_step_post, file= sprintf("%s/%s.est.Rdata", RESULTS_DIR, output_fn))
    } else {

    #print(str(all_data))
    # Define the Log-likelihood function
    time2_ref_allele_count <- all_data$MODEL_COUNT_2h
    time2_total_count <- all_data$total_counts_2h
    time6_ref_allele_count <- all_data$MODEL_COUNT_6h
    time6_total_count <- all_data$total_counts_6h


    gene_ind_list <- list(all_data$gene_name, all_data$gene_name)
    n_gene <- length(unique(unlist(gene_ind_list)))


    # Observed data
    x <- list(time2_ref_allele_count, time6_ref_allele_count)
    n <- list(time2_total_count, time6_total_count)

    # EM algorithm
    tic <- proc.time()
    e.step <- e_step(x, n, gene_ind_list, pi_init, A0_INIT, B0_INIT, A1_INIT, B1_INIT)
    #print(e.step)
    m.step <- m_step(x, n, gene_ind_list, e.step[["posterior.df"]], c(A0_INIT, A1_INIT)) # this step produces a bunch of warnings 

    #1: In lbeta(shape1[okk] + x[okk], shape2[okk] + size[okk] -  ... : NaNs produced
    #2: In lbeta(shape1[okk], shape2[okk]) : NaNs produced

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
    #raw_product_df = e.step[["raw_product.df"]]


    summarize_results = function(final_step_post){
        predicted_state = apply(final_step_post[,-1], MARGIN=1, which.max)
        max_probability = apply(final_step_post[,-1], MARGIN=1, max)
        # subtract 1 to account for indexing
        predicted_state = predicted_state - 1
        final_step_post$predicted_state = predicted_state
        final_step_post$max_probability = max_probability
        #summary_results = combine_dat %>% group_by(state) %>% count(predicted_state)
        #total_genes = summary_results %>% group_by(state) %>% summarize(total_genes = sum(n))
        #summary_results = inner_join(summary_results , total_genes, by = 'state') %>% mutate(prop = n/total_genes) %>% mutate(predicted_state = paste0('state', predicted_state))
        if(WHICH_MODEL == 'ASE_0h'){
            final_step_post$predicted_state = final_step_post$predicted_state + 3
        }
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
    final_step_param$filter_type = paste0('MAX_PROB_CUTOFF_', MAX_PROB_CUTOFF, '| STATE_QUANT_CUTOFF_', STATE_QUANT_CUTOFF)


    initial_parameters = c(A0_INIT, A1_INIT, B1_INIT, pi_init)
    final_step_post$filter_type = paste0('MAX_PROB_CUTOFF_', MAX_PROB_CUTOFF, '| STATE_QUANT_CUTOFF_', STATE_QUANT_CUTOFF)
    final_step_post$cell_line = CELL_LINE
    save(final_step_post, final_step_param, summary_results, initial_parameters, e.step, file= sprintf("%s/%s.est.Rdata", RESULTS_DIR, output_fn))
    }
    return(final_step_post)

}

nonASE_res = make_predictions('nonASE_0h')
ASE_res = make_predictions('ASE_0h')
