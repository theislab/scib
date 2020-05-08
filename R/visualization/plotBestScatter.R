library(ggplot2)

#' Make best summary scatter
#'
#' This function loads the performance scores for the top performing methods
#' (`loadBestsScores()`), makes the summary scatter plot (`plotBestScatter()`)
#' and saves the results in various formats. See individual function docs for
#'  details.
#'
#' @param score_files Passed to `loadBestScores()`
#'
#' @author Luke Zappia
makeBestScatter <- function(scores_files) {
    scores <- loadBestScores(scores_files)
    summary_plot <- plotBestScatter(scores)

    now_str <- format(lubridate::now(), "%Y%m%d_%H%M%S")

    ggsave(
        paste0(now_str, "_best_scatter.pdf"),
        summary_plot,
        # device = cairo_pdf,
        width = 210, height = 297, units = "mm"
    )

    ggsave(
        paste0(now_str, "_best_scatter.tiff"),
        summary_plot,
        device = "tiff",
        dpi = "retina",
        width = 210, height = 297, units = "mm"
    )

    ggsave(
        paste0(now_str, "_best_scatter.png"),
        summary_plot,
        # device = "png",
        # type = "cairo-png",
        dpi = "retina",
        width = 210, height = 297, units = "mm"
    )
}

#' Plot best summary scatter
#'
#' Produces a plot showing a summary of the performace of top ranked integration
#' methods summarized over a set of test datasets. The overall batch correction
#' score is shown on the x-axis and the overall bio conservation score on the
#' y-axis. Points show individual (versions of) methods.
#'
#' @param scores A tibble containing the scores for various methods on different
#' datasets. Created by `loadBestScores()`.
#'
#' @return A ggplot2 object
#'
#' @author Luke Zappia
plotBestScatter <- function(scores) {

    `%>%` <- magrittr::`%>%`

    scores_summ <- scores %>%
        dplyr::group_by(MethodVersion, Method, OutputFeatures, Output, Features,
                        Scaling) %>%
        dplyr::summarise(
            BatchMean = mean(`Batch Correction`),
            BatchSD   = sd(`Batch Correction`),
            BioMean   = mean(`Bio conservation`),
            BioSD     = sd(`Bio conservation`)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(MethodVersion)

    ggplot(scores_summ) +
        aes(
            x      = BatchMean,
            y      = BioMean,
            colour = `MethodVersion`
        ) +
        geom_errorbarh(
            aes(xmin = BatchMean - BatchSD, xmax = BatchMean + BatchSD)
        ) +
        geom_errorbar(
            aes(ymin = BioMean   - BioSD, ymax = BioMean   + BioSD)
        ) +
        geom_point(size = 3, stroke = 1, fill = "white") +
        scale_color_paletteer_d(
            palette = "ggsci::category20_d3",
            labels = glue::glue(
                "{scores_summ$Method} ({scores_summ$Output}, ",
                "{scores_summ$Features}, {scores_summ$Scaling})"
            )
        ) +
        coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
        labs(
            x = "Batch correction",
            y = "Bio conservation"
        ) +
        guides(
            colour = guide_legend(
                title          = "Method",
                title.position = "top",
                ncol           = 1,
                order          = 10
            )
        ) +
        theme_minimal() +
        theme(
            legend.position  = "right",
            axis.text        = element_text(size = 7),
            panel.border     = element_rect(fill = NA),
            strip.background = element_rect(fill = "black"),
            strip.text       = element_text(size = 10, colour = "white")
        )
}

#' Load best scores
#'
#' Read a set of scores files and return a single tibble with scores for the
#' selected best methods.
#'
#' @param scores_files Vector of paths to CSV files containing scores for each
#' dataset. These files are produced by the `plotSingleAtlas()` function. If
#' they are store in a directory named `data/` this vector can be created using
#' `fs::dir_ls("data")`.
#'
#' @return A tibble
#'
#' @author Luke Zappia
loadBestScores <- function(scores_files) {

    `%>%` <- magrittr::`%>%`

    dataset_key <- c(
        pancreas_jointnorm               = "Pancreas",
        lung_atlas                       = "Lung",
        immune_cell_hum                  = "Immune (human)",
        immune_cell_hum_mou              = "Immune (human/mouse)",
        mouse_brain                      = "Mouse brain",
        simulations_1_1                  = "Sim 1",
        simulations_2                    = "Sim 2",
        mouse_brain_atac_small_3datasets = "ATAC small",
        mouse_brain_atac_large_3datasets = "ATAC large"
    )

    best_methods <- c(
        "BBKNN_graph_HVG_unscaled",
        "Scanorama_embed_HVG_scaled",
        "Conos_graph_HVG_unscaled",
        "Scanorama_gene_HVG_scaled",
        "ComBat_gene_HVG_unscaled",
        "MNN_gene_HVG_scaled",
        "Harmony_embed_HVG_unscaled",
        "Seurat v3_gene_HVG_unscaled",
        "TrVAE_embed_HVG_unscaled",
        "LIGER_embed_HVG_unscaled",
        "scVI_embed_HVG_unscaled"
        # "Unintegrated_gene_FULL_unscaled"
    )

    scores_list <- purrr::map(scores_files, function(.file) {

        dataset <- .file %>%
            fs::path_file() %>%
            fs::path_ext_remove() %>%
            stringr::str_remove("_summary_scores")

        suppressWarnings(
            readr::read_csv(
                .file,
                col_types = readr::cols(
                    .default = readr::col_double(),
                    Method   = readr::col_character(),
                    Output   = readr::col_character(),
                    Features = readr::col_character(),
                    Scaling  = readr::col_character()
                )
            )) %>%
            dplyr::mutate(
                Dataset = dataset,
                Datatype = dplyr::if_else(
                    stringr::str_detect(.file, "atac"),
                    "ATAC",
                    "RNA"
                )
            )
    })

    best_scores <- scores_list %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(
            Dataset = factor(
                Dataset,
                levels = names(dataset_key),
                labels = dataset_key
            ),
            MethodVersion = paste(Method, Output, Features, Scaling, sep = "_"),
            OutputFeatures = paste(Output, Features, sep = "_")
        ) %>%
        dplyr::select(Dataset, Datatype, MethodVersion, Method, OutputFeatures,
                      Output, Features, Scaling, `Overall Score`,
                      `Batch Correction`, `Bio conservation`, -X1,
                      dplyr::everything()) %>%
        dplyr::filter(
            MethodVersion %in% best_methods,
            stringr::str_detect(Dataset, "ATAC", negate = TRUE),
            stringr::str_detect(Dataset, "Sim", negate = TRUE),
            !is.na(`Overall Score`)
        ) %>%
        dplyr::mutate(
            MethodVersion = factor(MethodVersion, levels = best_methods)
        )

    return(best_scores)
}
