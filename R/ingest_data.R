# Functions for ingesting phylogenetic trees and cleaning up names

#' Read in an association matrix
#' @param filepath path to association matrix (CSV format) with hosts as rownames and symbionts as colnames.
#' @export
read_associations = function(filepath) {

    # Read in the host-symbiont association matrix)
    assoc_mat = as.matrix(
        read.table(
            file = filepath,
            header = TRUE,
            check.names = FALSE
        )
    )

    assoc_info = list(
        assoc_matrix = assoc_mat,
        host_taxa = rownames(assoc_mat),
        symbiont_taxa = colnames(assoc_mat)
    )

    return(assoc_info)
}

#' Read in phylogenetic tree from a .nexus file
#' @param filepath path to the nexus file
#' @param nexus_format exus format to be pass to `phytools::readNexus()`, either 'standard' or 'raxml'
#' @param taxa_to_keep (optional) a list (vector) of taxa in the tree to keep, if left blank, all taxa will be kept
#' @param label_fix (optional) a character vector (length = 2) of strings that can be used to replace string prefixes/suffixes in the data; e.g., c(' ', '_') to replace spaces with underscores
#' @export
read_phylo = function(filepath, nexus_format, taxa_to_keep = NULL, label_fix = NULL) {

    # read in the nexus tree
    message("Reading nexus file.")
    phylo = phytools::readNexus(
        filepath,
        nexus_format
    )

    if (!is.null(label_fix)) {
        message(paste0("Replacing '", label_fix[1], "' in taxa names with '", label_fix[2], "'"))
        phylo[['tip.label']] = stringr::str_replace(
            phylo[['tip.label']],
            label_fix[1],
            label_fix[2]
        )
    }

    if (!is.null(taxa_to_keep)) {
        message(paste0("Keeping only taxa in provided list: ", paste0(taxa_to_keep, collapse = ', ')))
        phylo = ape::keep.tip(
            phylo,
            taxa_to_keep
        )
    }

    return(phylo)
}
