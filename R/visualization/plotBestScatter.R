library(ggplot2)

#' Make best summary scatter
#'
#' This function loads the performance scores for the top performing methods
#' (`loadBestScores()`), makes the summary scatter plot (`plotBestScatter()`)
#' and saves the results in various formats. See individual function docs for
#' details.
#'
#' @param score_files Passed to `loadBestScores()`
#' @param best_methods Passed to `loadBestScores()`
#' @param type Whether to plot RNA or ATAC results
#' @param out_dir Path to output directory
#'
#' @author Luke Zappia
makeBestScatter <- function(scores_files, best_methods, type = c("RNA", "ATAC"),
                            out_dir = ".") {

    type <- match.arg(type)

    scores <- loadBestScores(scores_files, best_methods, type)
    summary_plot <- plotBestScatter(scores)

    now_str <- format(lubridate::now(), "%Y%m%d_%H%M%S")

    for (extension in c("pdf", "tiff", "png")) {
        ggsave(
            as.character(glue::glue(
                "{out_dir}/{now_str}_best_scatter_{type}.{extension}"
            )),
            summary_plot,
            width = 210, height = 297, units = "mm"
        )
    }
}

#' Plot best summary scatter
#'
#' Produces a plot showing a summary of the performance of top ranked
#' integration methods summarized over a set of test datasets. The overall batch
#' correction score is shown on the x-axis and the overall bio conservation
#' score on the y-axis. Points show individual (versions of) methods.
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
        dplyr::group_by(MethodVersion, MethodLabel) %>%
        dplyr::summarise(
            BatchMean = mean(`Batch Correction`),
            BatchSD   = sd(`Batch Correction`),
            BioMean   = mean(`Bio conservation`),
            BioSD     = sd(`Bio conservation`),
            .groups   = "drop"
        ) %>%
        dplyr::arrange(MethodVersion)

    ggplot(scores_summ) +
        aes(
            x      = BatchMean,
            y      = BioMean,
            colour = MethodVersion
        ) +
        geom_errorbarh(
            aes(xmin = BatchMean - BatchSD, xmax = BatchMean + BatchSD)
        ) +
        geom_errorbar(
            aes(ymin = BioMean   - BioSD, ymax = BioMean   + BioSD)
        ) +
        geom_point(size = 3, stroke = 1, fill = "white") +
        paletteer::scale_color_paletteer_d(
            palette = "ggsci::category20_d3",
            labels = scores_summ$MethodLabel
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
                ncol           = 2,
                order          = 10
            )
        ) +
        theme_minimal() +
        theme(
            legend.position  = "bottom",
            axis.title       = element_text(size = 20),
            axis.text        = element_text(size = 14),
            panel.border     = element_rect(fill = NA),
            strip.background = element_rect(fill = "black"),
            strip.text       = element_text(size = 10, colour = "white"),
            legend.title     = element_blank(),
            legend.text      = element_text(size = 14)
        )
}

#' Load best scores
#'
#' Read a set of scores files and return a single tibble with scores for the
#' selected best methods.
#'
#' @param scores_files Vector of paths to CSV files containing scores for each
#' dataset. These files are produced by the `plotSingleTaskRNA()` or
#' `plotSingleTaskATAC()` functions. If they are stored in a directory named
#' `data/` this vector can be created using
#' `fs::dir_ls("data", glob = "*_summary_scores.csv")`.
#' @param best_methods Vector of names giving the best methods to plot. See the
#' output of `getBestMethods()` for an example
#' @param type Whether to load RNA or ATAC results
#'
#' @return A tibble
#'
#' @author Luke Zappia
loadBestScores <- function(scores_files, best_methods, type = c("RNA", "ATAC")) {

    type <- match.arg(type)

    `%>%` <- magrittr::`%>%`

    is_rna <- type == "RNA"
    files_keep <- stringr::str_detect(scores_files, "atac", negate = is_rna)
    scores_files <- scores_files[files_keep]

    scores_cols_types <- switch (type,
        RNA = readr::cols(
            .default = readr::col_double(),
            Method   = readr::col_character(),
            Output   = readr::col_character(),
            Features = readr::col_character(),
            Scaling  = readr::col_character()
        ),
        ATAC = readr::cols(
            .default = readr::col_double(),
            Method   = readr::col_character(),
            Output   = readr::col_character()
        )
    )

    scores_list <- purrr::map(scores_files, function(.file) {

        dataset <- .file %>%
            fs::path_file() %>%
            fs::path_ext_remove() %>%
            stringr::str_remove("_summary_scores")

        suppressWarnings(
            readr::read_csv(
                .file,
                col_types = scores_cols_types
        )) %>%
            dplyr::mutate(Dataset = dataset)
    })

    best_scores <- scores_list %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(
            OutputLabels = factor(
                Output,
                levels = c("embed", "gene", "graph"),
                labels = c(
                    "Embedding",
                    ifelse(type == "RNA", "Genes", "Features"),
                    "Graph"
                )
            )
        )

    best_scores <- switch (type,
        RNA = best_scores %>%
            dplyr::mutate(
                FeaturesLabels = factor(
                    Features,
                    levels = c("FULL", "HVG"),
                    labels = c("Full", "HVG")
                ),
                ScalingLabels = factor(
                    Scaling,
                    levels = c("scaled", "unscaled"),
                    labels = c("Scaled", "Unscaled")
                ),
                MethodVersion = paste(
                    Method, Output, Features, Scaling,
                    sep = "_"
                ),
                MethodLabel   = glue::glue(
                    "{Method} ({OutputLabels}, {FeaturesLabels}, {ScalingLabels})"
                )
            ),
        ATAC = best_scores %>%
            dplyr::mutate(
                MethodVersion = paste(Method, Output, sep = "_"),
                MethodLabel   = glue::glue("{Method} ({OutputLabels})")
            )
    )

    best_scores <- best_scores %>%
        dplyr::mutate(
            MethodVersion = factor(MethodVersion, levels = best_methods)
        ) %>%
        dplyr::filter(
            MethodVersion %in% best_methods,
            stringr::str_detect(Dataset, "Sim", negate = TRUE),
            stringr::str_detect(Dataset, "3batches", negate = TRUE),
            !is.na(`Batch Correction`),
            !is.na(`Bio conservation`)
        ) %>%
        dplyr::select(
            Dataset, MethodVersion, MethodLabel, `Batch Correction`,
            `Bio conservation`, Output
        )

    return(best_scores)
}
