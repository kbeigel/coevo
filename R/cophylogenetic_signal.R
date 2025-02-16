# Functions for running tests of cophylogenetic signal (ParaFit, PACo).

#' A single run of the parafit function
#' Based on https://github.com/alixmatthews/cophylogenetic/blob/master/parafit_script.R
#' @export
single_parafit_run = function(i, host_dist, symbiont_dist, assoc_mat, correction_method) {
    
    res = ape::parafit(
        host.D = host_dist, 
        para.D = symbiont_dist,
        HP = assoc_mat,
        correction = correction_method,
        nperm = 999,
        test.links = TRUE
    )

    parafit_global_fit = res$ParaFitGlobal
    parafit_global_p_val = res$p.global
    link_table = res$link.table

    single_run_res = list(
        parafit_global_fit = parafit_global_fit,
        parafit_global_p_val = parafit_global_p_val,
        link_table = link_table   
    )

    return(single_run_res)
}

#' Multiple parafit runs.
#' Use lapply to run the single_parafit_run() multiple times (n_runs)
#' Based on https://github.com/alixmatthews/cophylogenetic/blob/master/parafit_script.R
#' @param host_dist host distance matrix
#' @param symbiont_dist symbiont distance matrix
#' @param assoc_mat association matrix, with hosts as rownames and symbionts as colnames
#' @return a list of results from individual parafit runs (`single_parafit_run()``)
#' @export
run_parafit_loop = function(host_dist, symbiont_dist, assoc_mat, n_runs, correction_method) {

    seed_list = sample(1:10000, n_runs)

    message(paste0("Running parafit ", n_runs, " times."))
    res = parallel::mclapply(
                X = 1:n_runs,
                FUN = single_parafit_run,
                host_dist = host_dist,
                symbiont_dist = symbiont_dist,
                assoc_mat = assoc_mat,
                correction = correction_method,
                mc.preschedule = TRUE,
                mc.cores = 10,
                mc.set.seed = TRUE
            )

    print(res)
    return(res)
}

#' Generate summary stats (pval means, etc.) for the results of run_parafit_loop()
#' Based on https://github.com/alixmatthews/cophylogenetic/blob/master/parafit_script.R
#' @param res results from run_parafit_loop
#' @param assoc_mat association matrix, with hosts as rownames and symbionts as colnames
#' @export
gather_parafit_stats = function(res, assoc_mat) {

    # Global fit metric
    parafit_global_fit_vals = sapply(res, function(x) x$parafit_global_fit)
    parafit_global_mean = mean(parafit_global_fit_vals)

    # Mean global p-value
    parafit_global_pvals = sapply(res, function(x) x$parafit_global_p_val)
    parafit_global_pval_mean = mean(parafit_global_pvals)

    link_order = cbind(
        rownames(assoc_mat)[res[[1]]$link_table[, 1]],
        colnames(assoc_mat)[res[[1]]$link_table[, 2]]
    )

    p_F1_matrix = do.call(rbind, lapply(res, function(x) x$link_table[, 4]))
    p_F2_matrix = do.call(rbind, lapply(res, function(x) x$link_table[, 6]))

    link_adjpvals_F1 = apply(p_F1_matrix, 2, p.adjust, method = 'BH')
    link_adjpvals_F1_means = colMeans(link_adjpvals_F1)

    link_adjpvals_F2 = apply(p_F2_matrix, 2, p.adjust, method = 'BH')
    link_adjpvals_F2_means = colMeans(link_adjpvals_F2)

    mean_link_adjpvals = cbind(link_order, link_adjpvals_F1_means, link_adjpvals_F2_means)
    colnames(mean_link_adjpvals) = c('host', 'symbiont', 'mean_p_val_adj_F1', 'mean_p_val_adj_F2')

    all_res = list(
        parafit_runs = res,
        parafit_global_mean = parafit_global_mean,
        parafit_global_pval_mean = parafit_global_pval_mean,
        link_adjpvals_F1 = link_adjpvals_F1,
        link_adjpvals_F2 = link_adjpvals_F2,
        mean_link_adjpvals = as.data.frame(mean_link_adjpvals)
    )

    return(all_res)
}


#' @export
parafit_analysis = function(host_tree, symbiont_tree, assoc_mat, n_runs, correction_method) {

    message('Taxa names in host tree:')
    message(paste0(host_tree$tip.label, collapse = ', '))

    message('Taxa names for host in association matrix (row names):')
    message(paste0(assoc_mat$host_taxa, collapse = ', '))

    message('Taxa names in symbiont tree:')
    message(paste0(symbiont_tree$tip.label, collapse = ', '))

    message('Taxa names for symbiont in association matrix (row names):')
    message(paste0(assoc_mat$symbiont_taxa, collapse = ', '))
    
    message('Getting cophenetic distances for host tree.')
    host_dist = get_cophenetic(host_tree, assoc_mat$host_taxa)

    message('Getting cophenetic distances for symbiont tree.')
    symbiont_dist = get_cophenetic(symbiont_tree, assoc_mat$symbiont_taxa)

    # Run ParaFit
    res = run_parafit_loop(host_dist, symbiont_dist, assoc_mat$assoc_mat, n_runs, correction_method)
    stats = gather_parafit_stats(res, assoc_mat$assoc_mat)

    return(stats)

}


#' PACO functions
#' Transforms the host and parasite distance matrices into the respective matrices of Principal Coordinates (ape::pcoa)
#' and duplicates taxa (if necessary) to accommodate multiple host-parasite associations
#' @export
get_pcoas = function(host_dist, symbiont_dist, assoc_mat) {
    
    hs_bin = which(assoc_mat > 0, arr.in = TRUE)

    # Host distance Principal COordinates Analysis
    h_pcoa = ape::pcoa(host_dist, correction = "cailliez")$vectors
    h_pcoa = h_pcoa[hs_bin[, 1], ] # adjust host PCoA vectors

    # Host distance Principal COordinates Analysis
    s_pcoa <- ape::pcoa(symbiont_dist, correction = "cailliez")$vectors 
    s_pcoa <- s_pcoa[hs_bin[, 2], ]  ##adjust Parasite PCo vectors

    pcoa_res = list(host_pcoa = h_pcoa, symbiont_pcoa = s_pcoa)

    return(pcoa_res)
}

#' 
#' @export
run_paco = function(host_dist, symbiont_dist, assoc_mat, n_perm) {

    pcoa_list = get_pcoas(host_dist, symbiont_dist, assoc_mat$assoc_mat)

    # Run procrustes once to get the observed m2
    procru = vegan::procrustes(pcoa_list$host_pcoa, pcoa_list$symbiont_pcoa)
    n_links = sum(assoc_mat$assoc_mat)

    # Get sum of squares from procrustes() run
    m2_obs = procru$ss

    # set initial p_value
    p_value = 0
    set.seed(45789)

    m2_perm_list = list()

    hs_mat = assoc_mat$assoc_mat

    for (n in c(1:n_perm)) {
        # if the number of links is less than or equal to the number of hosts or the number of symbionts
        # (avoid all symbionts associating with a single host in the permutation)
        if (n_links <= nrow(hs_mat) | n_links <= ncol(hs_mat)) {
            flag = TRUE
            while (flag == TRUE) { 
                hs_perm = t(apply(hs_mat, 1, sample))
                if (any(colSums(hs_perm) == n_links)) {
                    flag = TRUE 
                } else {
                    flag = FALSE
                }
            }
        } else { 
            hs_perm = t(apply(hs_mat, 1, sample)) # permutes each HP row independently
        }
        
        # Running PACo
        paco_perm = get_pcoas(host_dist, symbiont_dist, hs_perm)
        m2_perm = vegan::procrustes(paco_perm$host_pcoa, paco_perm$symbiont_pcoa)$ss # permuted sum of squares

        m2_perm_list[[n]] = m2_perm

        if (m2_perm <= m2_obs) {
            p_value = p_value + 1
        }
    }

    p_value = p_value / n_perm

    paco_res = list(
        m2_obs = m2_obs,
        m2_perm_list = m2_perm_list,
        p_value = p_value,
        n_perm = n_perm
    )

    # Get jackknife residuals
    jk_res = jackknife_calc(host_dist, symbiont_dist, assoc_mat, procru)

    paco_jk_res = list(paco_res = paco_res, jackknife_res = jk_res)

    return(paco_jk_res)

}

jackknife_calc = function(host_dist, symbiont_dist, assoc_mat, procru) {

    hs_ones = which(assoc_mat$assoc_mat > 0, arr.in = TRUE)
    n_links = sum(assoc_mat$assoc_mat)

    # empty matrix of jackknifed squared residuals
    sqres_jackkn = matrix(rep(NA, n_links**2), n_links)
    colnames(sqres_jackkn) = paste(rownames(procru$X), rownames(procru$Yrot), sep = "-")

    t_crit = qt(0.975, n_links - 1)

    for (i in c(1:n_links)) {
        hs_ind = assoc_mat$assoc_mat
        hs_ind[hs_ones[i, 1], hs_ones[i, 2]] = 0
        paco_ind = get_pcoas(host_dist, symbiont_dist, hs_ind)
        proc_ind = vegan::procrustes(paco_ind$host_pcoa, paco_ind$symbiont_pcoa)
        res_proc_ind = c(residuals(proc_ind))
        res_proc_ind = append(res_proc_ind, NA, after = i - 1)
        sqres_jackkn[i, ] = res_proc_ind
    }

    # Jackknifed residuals are squared
    sqres_jackkn <- sqres_jackkn**2

    # Vector of original square residuals
    sqres <- (residuals(procru))**2

    #jackknife calculations:
    sqres_jackkn = sqres_jackkn * (-(n_links - 1))
    sqres = sqres * n_links
    sqres_jackkn = t(apply(sqres_jackkn, 1, "+", sqres)) # apply jackknife function to matrix
    phi_mean = apply(sqres_jackkn, 2, mean, na.rm = TRUE) # mean jackknife estimate per link
    st_dev = apply(sqres_jackkn, 2, sd, na.rm = TRUE) # standard deviation of estimates
    phi_uci = phi_mean + t_crit * st_dev / sqrt(n_links) # upper 95% confidence interval

    res = list(phi_mean = phi_mean, phi_uci = phi_uci, sqres_jackkn = sqres_jackkn)

    return(res)
}

#' 
#' @export
paco_analysis =  function(host_tree, symbiont_tree, assoc_mat, n_perm) {

    message('Taxa names in host tree:')
    message(paste0(host_tree$tip.label, collapse = ', '))

    message('Taxa names for host in association matrix (row names):')
    message(paste0(assoc_mat$host_taxa, collapse = ', '))

    message('Taxa names in symbiont tree:')
    message(paste0(host_tree$tip.label, collapse = ', '))
    
    message('Getting cophenetic distances for host tree.')
    host_dist = get_cophenetic(host_tree, assoc_mat$host_taxa)

    message('Getting cophenetic distances for symbiont tree.')
    symbiont_dist = get_cophenetic(symbiont_tree, assoc_mat$symbiont_taxa)

    message('Running PACo.')
    paco_res = run_paco(host_dist, symbiont_dist, assoc_mat, n_perm)

    return(paco_res)

}

#' 
#' @export
determine_pval_comp = function(padj_mat, p_alpha) {

    padj_comp = vector()

    # Loop over each column to determine if pvals were all sig, not sig, or a mix
    for (i in 1:ncol(padj_mat)) {

        padj_vals = padj_mat[, i]

        if (all(padj_vals < p_alpha)) {
            padj_comp[[i]] = 'all_significant'
        } else if (all(padj_vals > p_alpha)) {
            padj_comp[[i]] = 'none_significant'
        } else {
            padj_comp[[i]] = 'mixed_significant'
        }
    }

    return(padj_comp)
}

#' 
#' @export
compile_global_results = function(parafit_results, paco_results) {

    res_df = data.frame(
        Analysis = character(),
        Metric = character(),
        Value = numeric(),
        P_Value_Metric = character(),
        P_Value = numeric(),
        N_permutations = numeric(),
        N_runs = numeric()
    )
    
    pf_res = c(
        'ParaFit', 'ParaFit_Global_Fit', parafit_results$parafit_global_mean,
        'ParaFit_Global_Mean_Pvalue', parafit_results$parafit_global_pval_mean,
        999, length(parafit_results$parafit_runs))

    res_df = structure(rbind(res_df, pf_res), .Names = names(res_df))

    paco_res = c(
        'PACo', 'PACo_Observed_m2', paco_results$paco_res$m2_obs,
        'PACo_Pvalue', paco_results$paco_res$p_value,
        paco_results$paco_res$n_perm, 1)

    res_df = structure(rbind(res_df, paco_res), .Names = names(res_df))
    
    return(res_df)

}

#' 
#' @export
compile_link_results = function(parafit_results, paco_results, p_alpha) {

    padj_F1 = parafit_results$link_adjpvals_F1
    padj_F2 = parafit_results$link_adjpvals_F2

    ParaFit_F1_padj_comp = determine_pval_comp(padj_F1, p_alpha)
    ParaFit_F2_padj_comp = determine_pval_comp(padj_F2, p_alpha)

    res_df = data.frame(
        Host = parafit_results$mean_link_adjpvals$host,
        Symbiont = parafit_results$mean_link_adjpvals$symbiont,
        ParaFit_F1_mean_adj_pval = parafit_results$mean_link_adjpvals$mean_p_val_adj_F1,
        ParaFit_F2_mean_adj_pval = parafit_results$mean_link_adjpvals$mean_p_val_adj_F2,
        ParaFit_F1_padj_comp = ParaFit_F1_padj_comp,
        ParaFit_F2_padj_comp = ParaFit_F2_padj_comp,
        PACo_Jackknife_Mean = paco_results$jackknife_res$phi_mean,
        PACo_Jackknife_UCI = paco_results$jackknife_res$phi_uci
    )

    rownames(res_df) = NULL

    res_df = dplyr::mutate(res_df, 
        ParaFit_F1_sig = dplyr::case_when(
            ParaFit_F1_mean_adj_pval < p_alpha ~ TRUE,
            .default = FALSE
        ),
        ParaFit_F2_sig = dplyr::case_when(
            ParaFit_F2_mean_adj_pval < p_alpha ~ TRUE,
            .default = FALSE
        ),
        PACo_Jackknife_sig = dplyr::case_when(
            PACo_Jackknife_UCI < median(paco_results$jackknife_res$phi_mean) ~ TRUE,
            .default = FALSE
        )
    )

    return(res_df)

}