#' Calculate cophenetic distances for a phylogenetic tree
#' 
#' This function gets the cophenetic distances of a phylo class object.
#' 
#' @param phylo object of class phylo
#' @param taxa_sort character vector of the taxa, in order, for which to get distances.
#' @export
get_cophenetic = function(phylo, taxa_sort) {

    cophenetic_dist = as.matrix(cophenetic(phylo))[taxa_sort, taxa_sort]

    return(cophenetic_dist)
}
