# Functions for running tests of cophylogenetic signal (ParaFit, PACo).

#' A single run of the parafit function
#' Based on https://github.com/alixmatthews/cophylogenetic/blob/master/parafit_script.R
#' @export
single_parafit_run = function(i, host_dist, symbiont_dist, assoc_mat) {
    
    tmp = ape::parafit(
        host.D = host_dist, 
        para.D = symbiont_dist,
        HP = assoc_mat,
        nperm = 999,
        test.links = TRUE,
        correction = 'lingoes'
        )

    return(tmp)
}

#' Multiple parafit runs.
#' Use lapply to run the single_parafit_run() multiple times (n_runs)
#' Based on https://github.com/alixmatthews/cophylogenetic/blob/master/parafit_script.R
#' @param host_dist host distance matrix
#' @param symbiont_dist symbiont distance matrix
#' @param assoc_mat association matrix, with hosts as rownames and symbionts as colnames
#' @return a list of results from individual parafit runs (`single_parafit_run()``)
#' @export
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

#' Generate summary stats (pval means, etc.) for the results of run_parafit_loop()
#' Based on https://github.com/alixmatthews/cophylogenetic/blob/master/parafit_script.R
#' @param res results from run_parafit_loop
#' @param assoc_mat association matrix, with hosts as rownames and symbionts as colnames
#' @export
gather_parafit_stats = function(res, assoc_mat) {
    # Mean global p-value
    mean_global_pval = mean(sapply(X = as.matrix(res), FUN = '[[', n = 2))

    # Get all links
    res_links = c()
    for (run in 1:length(res)) { 
        res_links = cbind(res_links, res[[run]])
    }

    link_tbl = res_links[, 1]$link.table
    link_order = cbind(
        rownames(assoc_mat)[link_tbl[, 1]],
        colnames(assoc_mat)[link_tbl[, 2]]
    )

    link_pvals = NULL
    for (run in 1:ncol(res_links)) {
        # link.table column 4 = p.F1, link.table column 6 = p.F2
        p_F1 = res_links[, run]$link.table[, 4]
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


#' @export
parafit_analysis = function(host_tree, symbiont_tree, assoc_mat, n_runs) {

    message('Taxa names in host tree:')
    message(paste0(host_tree$tip.label, collapse = ', '))

    message('Taxa names for host in association matrix (row names):')
    message(paste0(assoc_mat$host_taxa, collapse = ', '))

    message('Taxa names in symbiont tree:')
    message(paste0(host_tree$tip.label, collapse = ', '))
    
    
    message('Getting cophenetic distances for host tree.')
    host_dist = get_cophenetic(host_tree, assoc_mat$host_taxa)

    message('Getting cophenetic distances for host tree.')
    symbiont_dist = get_cophenetic(symbiont_tree, assoc_mat$symbiont_taxa)

    # Run ParaFit
    res = run_parafit_loop(host_dist, symbiont_dist, assoc_mat$assoc_mat, n_runs)
    stats = gather_parafit_stats(res, assoc_mat)

}


# PACO functions
# Transforms the host and parasite distance matrices into the respective matrices of Principal Coordinates (ape::pcoa)
# and duplicates taxa (if necessary) to accommodate multiple host-parasite associations
#' @export
get_pcoas = function(host_dist, symbiont_dist, assoc_mat) {
    
    hs_bin = which(assoc_mat$assoc_mat > 0, arr.in = TRUE)

    # Host distance Principal COordinates Analysis
    h_pcoa = ape::pcoa(host_dist, correction = "cailliez")$vectors
    h_pcoa = h_pcoa[hs_bin[, 1], ] # adjust host PCoA vectors

    # Host distance Principal COordinates Analysis
    s_pcoa <- ape::pcoa(symbiont_dist, correction = "cailliez")$vectors 
    s_pcoa <- s_pcoa[HP.bin[, 2], ]  ##adjust Parasite PCo vectors

    return(list(host_pcoa = h_pcoa, symbiont_pcoa = s_pcoa))
}

#' @export
run_paco = function(pcoa_list, n_perm) {
    # Run procrustes once to get the observed m2
    procru = procrustespcoa_list
}