#' Get best methods
#'
#' Provides example of the `best_methods` vector
#'
#' @param type Whether to return keys for RNA or ATAC
#'
#' @return Vector of method names
#'
#' @author Luke Zappia
getBestMethods <- function(type = c("RNA", "ATAC")) {

    type <- match.arg(type)

    best_methods <- list(
        RNA = c(
            "Scanorama_embed_HVG_scaled",
            "scANVI*_embed_HVG_unscaled",
            "fastMNN_embed_HVG_unscaled",
            "scGen*_gene_HVG_unscaled",
            "BBKNN_graph_HVG_unscaled",
            "Scanorama_gene_HVG_scaled",
            "scVI_embed_HVG_unscaled",
            "Seurat v3 RPCA_gene_HVG_unscaled",
            "Harmony_embed_HVG_unscaled",
            "Seurat v3 CCA_gene_HVG_unscaled",
            "fastMNN_gene_HVG_unscaled",
            "Conos_graph_FULL_unscaled",
            "ComBat_gene_HVG_unscaled",
            "MNN_gene_HVG_scaled",
            "trVAE_embed_HVG_unscaled",
            "DESC_embed_FULL_unscaled",
            "LIGER_embed_HVG_unscaled",
            "SAUCIE_embed_HVG_scaled",
            "SAUCIE_gene_HVG_scaled"
        ),
        ATAC = c(
            "Harmony_embed",
            "scANVI*_embed",
            "LIGER_embed",
            "scVI_embed",
            "ComBat_gene",
            "scGen*_gene",
            "Seurat v3 RPCA_gene",
            "Seurat v3 CCA_gene",
            "Scanorama_embed",
            "MNN_gene",
            "trVAE_embed",
            "fastMNN_embed",
            "BBKNN_graph",
            "Scanorama_gene",
            "fastMNN_gene",
            "DESC_embed",
            "SAUCIE_embed",
            "SAUCIE_gene",
            "Conos_graph"
        )
    )

    switch (type,
        RNA = best_methods$RNA,
        ATAC = best_methods$ATAC
    )
}

#' Get dataset key
#'
#' Provides example fo the `dataset_key` vector
#'
#' @param type Whether to return keys for RNA or ATAC
#'
#' @return Named dataset key vector
#'
#' @author Luke Zappia
getDatasetKey <- function(type = c("RNA", "ATAC")) {

    type <- match.arg(type)

    dataset_key <- list(
        RNA  = c(
            pancreas                       = "Pancreas",
            lung_atlas                     = "Lung",
            immune_cell_hum                = "Immune (human)",
            immune_cell_hum_mou            = "Immune (human/mouse)",
            mouse_brain                    = "Mouse brain",
            simulations_1_1                = "Sim 1",
            simulations_2                  = "Sim 2"
        ),
        ATAC = c(
            mouse_brain_atac_genes_large   = "ATAC large (genes)",
            mouse_brain_atac_peaks_large   = "ATAC large (peaks)",
            mouse_brain_atac_windows_large = "ATAC large (windows)",
            mouse_brain_atac_genes_small   = "ATAC small (genes)",
            mouse_brain_atac_peaks_small   = "ATAC small (peaks)",
            mouse_brain_atac_windows_small = "ATAC small (windows)"
        )
    )


    switch (type,
        RNA = dataset_key$RNA,
        ATAC = dataset_key$ATAC
    )
}

#' Get methods palette
#'
#' Example of a vector specifying colour palette for methods used by
#' `plotSummaryScatter()`
#'
#' @return Named vector where names are methods and values are colours
#'
#' @author Luke Zappia
getMethodsPal <- function() {
    c(
        "BBKNN"          = "#5A5156",
        "Conos"          = "#E4E1E3",
        "trVAE"          = "#F6222E",
        "scVI"           = "#FE00FA",
        "ComBat"         = "#16FF32",
        "Harmony"        = "#3283FE",
        "LIGER"          = "#FEAF16",
        "Scanorama"      = "#B00068",
        "Seurat v3 CCA"  = "#1CFFCE",
        "Seurat v3 RPCA" = "#90AD1C",
        "MNN"            = "#2ED9FF",
        "fastMNN"        = "#DEA0FD",
        "scGen*"         = "#AA0DFE",
        "scANVI*"        = "#F8A19F",
        "DESC"           = "#325A9B",
        "SAUCIE"         = "#C4451C",
        "Unintegrated"   = "#66B0FF"
    )
}
