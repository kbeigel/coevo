# Functions for running tests of cophylogenetic signal (ParaFit).

# Calculate cophenetic distances for a phylogenetic tree
get_cophenetic = function(phylo, taxa_sort) {

    cophenetic_dist = as.matrix(cophenetic(phylo))[taxa_sort, taxa_sort]

    return(cophenetic_dist)
}

# A single run of the parafit function
# Based on https://github.com/alixmatthews/cophylogenetic/blob/master/parafit_script.R
single_parafit_run = function(i, host_dist, symbiont_dist, assoc_mat) {
    
    tmp = ape::parafit(
        host.D = host_dist, 
        para.D = symbiont_dist,
        HP = assoc_mat,
        nperm = 999,
        test.links = TRUE,
        correction = 'lingoes'
        )

    return (tmp)
}

# Based on https://github.com/alixmatthews/cophylogenetic/blob/master/parafit_script.R
# Use lapply to run the single_parafit_run() multiple times (n_runs)
run_parafit_loop = function(host_dist, symbiont_dist, assoc_mat, n_runs) {

    message(paste0("Running parafit ", n_runs, " times."))
    res = lapply(1:n_runs, function(i) {
        simplify2array(
            parallel::mclapply(
                X = 1,
                FUN = single_parafit_run,
                host_dist = host_dist,
                symbiont_dist = symbiont_dist,
                assoc_mat = assoc_mat,
                mc.preschedule = TRUE,
                mc.cores = 1,
                mc.set.seed = FALSE
            )
        )
    })

    return(res)
}

# Generate summary stats (pval means, etc.) for the results of run_parafit_loop()
# Based on https://github.com/alixmatthews/cophylogenetic/blob/master/parafit_script.R
gather_stats = function(res, assoc_mat) {
    # Mean global p-value
    mean_global_pval = mean(sapply(X = as.matrix(res), FUN = '[[', n = 2))

    # Get all links
    res_links = c()
    for (run in 1:length(res)){ 
        res_links = cbind(res_links, res[[run]])
    }

    link_tbl = res_links[, 1]$link.table
    link_order = cbind(
        rownames(assoc_mat)[link_tbl[,1]],
        colnames(assoc_mat)[link_tbl[,2]]
    )

    link_pvals = NULL
    for (run in 1:ncol(res_links)) {
        # link.table column 4 = p.F1, link.table column 6 = p.F2
        p_F1 = res_links[, run]$link.table[,4]
        link_pvals = rbind(link_pvals, p_F1)
    }

    link_adjpvals = apply(link_pvals, 2, p.adjust, method = 'BH')
    link_adjpval_means = colMeans(link_adjpvals)
    mean_link_adjpvals = cbind(link_order, link_adjpval_means)

    all_res = list(
        parafit_runs = res,
        mean_global_pval = mean_global_pval,
        mean_link_adjpvals = mean_link_adjpvals
        )

    return(all_res)
}
