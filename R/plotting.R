# Functions for plotting

#' @export
plot_cophylo = function(h_phy, s_phy, results, filename, plot_h, plot_w, palette) {

    palette = c(parafit = palette[1], paco = palette[2], parafit_and_paco = palette[3])

    line_colors = dplyr::mutate(
        results, 
        color = dplyr::case_when(
            ParaFit_F1_sig == TRUE & PACo_Jackknife_sig == FALSE ~ palette[['parafit']],
            ParaFit_F1_sig == FALSE & PACo_Jackknife_sig == TRUE ~ palette[['paco']],
            ParaFit_F1_sig == TRUE & PACo_Jackknife_sig == TRUE ~ palette[['parafit_and_paco']],
            .default = 'gray30'
            )
        ) |>
        dplyr::pull(color)

    line_style = dplyr::mutate(
        results,
        style = dplyr::case_when(
            ParaFit_F1_padj_comp == 'all_significant'| PACo_Jackknife_sig == TRUE ~ 'solid',
            ParaFit_F1_padj_comp == 'mixed_significant' ~ 'dotdash',
            ParaFit_F1_padj_comp == 'none_significant' | PACo_Jackknife_sig == FALSE ~ 'longdash',
            .default = 'gray20'
            )
        ) |>
        dplyr::pull(style)

    hs_links = data.frame(parafit_results$mean_link_adjpvals$symbiont, parafit_results$mean_link_adjpvals$host)
    cophy = phytools::cophylo(s_phy, h_phy, assoc = hs_links, rotate = TRUE)

    pdf(file = paste0(filename, '.pdf'), height = plot_h, width = plot_w)

    plot(
        cophy, 
        link.col = line_colors, link.lty = line_style, 
        link.lwd = 2.5, fsize = 0.8, lwd = 2, pts = TRUE
    )

    # ADD NODE LABELS (BOOTSTRAP VALUES)
    # METHOD ADAPTED FROM:
    # http://blog.phytools.org/2018/06/preserving-node-edge-labels-after.html
    phytools::nodelabels.cophylo(cophy$trees[[1]]$node.label, 1:cophy$trees[[1]]$Nnode +
                    ape::Ntip(cophy$trees[[1]]),
                    which = "left", frame = "none",
                    cex = 1, adj = c(1.2, 1.2))

    phytools::nodelabels.cophylo(cophy$trees[[2]]$node.label, 1:cophy$trees[[2]]$Nnode +
                    ape::Ntip(cophy$trees[[2]]),
                    which = "right", frame = "none",
                    cex = 1, adj = c(-0.1, 1.3))

    dev.off()

}
