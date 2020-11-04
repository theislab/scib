library(ggplot2)

#' Make summary scatter
#'
#' This function loads the performance scores (`loadScores()`), makes the
#' summary scatter plot (`plotSummaryScatter()`) and saves the results in
#' various formats. See individual function docs for details.
#'
#' @param score_files Passed to `loadScores()`
#'
#' @author Luke Zappia
makeSummaryScatter <- function(scores_files) {
    scores <- loadScores(scores_files)
    summary_plot <- plotSummaryScatter(scores)

    now_str <- format(lubridate::now(), "%Y%m%d_%H%M%S")

    ggsave(
        paste0(now_str, "_scatter_summary.pdf"),
        summary_plot,
        # device = cairo_pdf,
        width = 297, height = 210, units = "mm"
    )

    ggsave(
        paste0(now_str, "_scatter_summary.tiff"),
        summary_plot,
        device = "tiff",
        dpi = "retina",
        width = 297, height = 210, units = "mm"
    )

    ggsave(
        paste0(now_str, "_scatter_summary.png"),
        summary_plot,
        # device = "png",
        # type = "cairo-png",
        dpi = "retina",
        width = 297, height = 210, units = "mm"
    )
}

#' Plot summary scatter
#'
#' Produces a plot showing a summary of the performace of integration methods
#' on a set of test datasets. The overall batch correction score is shown on the
#' x-axis and the overall bio conservation score on the y-axis. Points show
#' individual (versions of) methods.
#'
#' @param scores A tibble containing the scores for various methods on different
#' datasets. Created by `loadScores()`.
#'
#' @return A ggplot2 object
#'
#' @author Luke Zappia
plotSummaryScatter <- function(scores) {

    `%>%` <- magrittr::`%>%`

    methods_pal <- c(
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

    scores_int <- scores %>%
        dplyr::filter(Method != "Unintegrated") %>%
        dplyr::mutate(
            OutputFeatures = dplyr::case_when(
                OutputFeatures == "embed_NA" ~ "embed_FULL",
                OutputFeatures == "graph_NA" ~ "graph_FULL",
                OutputFeatures == "gene_NA"  ~ "gene_FULL",
                TRUE                         ~ OutputFeatures
            ),
            Scaling = dplyr::if_else(is.na(Scaling), "unscaled", Scaling)
        )

    medians <- scores_int %>%
        dplyr::group_by(Dataset) %>%
        dplyr::summarise(
            `Batch Correction` = median(`Batch Correction`),
            `Bio conservation` = median(`Bio conservation`)
        ) %>%
        dplyr::mutate(Type = "Median")

    unintegrated <- scores %>%
        dplyr::filter(Method == "Unintegrated") %>%
        dplyr::mutate(Type = "Unintegrated") %>%
        dplyr::select(Dataset, Type, `Batch Correction`, `Bio conservation`)

    ref_lines <- dplyr::bind_rows(unintegrated, medians)

    ggplot(scores_int) +
        aes(
            x      = `Batch Correction`,
            y      = `Bio conservation`,
            colour = `Method`,
            size   = `Overall Score`,
            shape  = `OutputFeatures`
        ) +
        geom_hline(
            data = dplyr::filter(ref_lines, Type == "Unintegrated"),
            aes(yintercept = `Bio conservation`, linetype = Type),
            colour = "red"
        ) +
        geom_vline(
            data = dplyr::filter(ref_lines, Type == "Unintegrated"),
            aes(xintercept = `Batch Correction`, linetype = Type),
            colour = "red"
        ) +
        geom_hline(
            data = dplyr::filter(ref_lines, Type == "Median"),
            aes(yintercept = `Bio conservation`, linetype = Type),
            colour = "blue"
        ) +
        geom_vline(
            data = dplyr::filter(ref_lines, Type == "Median"),
            aes(xintercept = `Batch Correction`, linetype = Type),
            colour = "blue"
        ) +
        geom_point(stroke = 1, fill = "white") +
        geom_point(
            data = dplyr::filter(
                scores_int,
                stringr::str_detect(OutputFeatures, "FULL")
            ),
            aes(alpha = Scaling),
            shape = 4, size = 1.5, colour = "white"
        ) +
        geom_point(
            data = dplyr::filter(
                scores_int,
                stringr::str_detect(OutputFeatures, "HVG")
            ),
            aes(alpha = Scaling),
            shape = 4, size = 1.5
        ) +
        scale_colour_manual(values = methods_pal) +
        scale_size_continuous(range = c(0.5, 3)) +
        scale_shape_manual(
            values = c(16, 21, 15, 22, 17, 24),
            labels = c("Embedding (Full)", "Embedding (HVG)", "Features (Full)",
                       "Features (HVG)", "Graph (Full)", "Graph (HVG)"),
            breaks = c("embed_FULL", "embed_HVG", "gene_FULL",
                       "gene_HVG", "graph_FULL", "graph_HVG"),
            drop = FALSE
        ) +
        scale_alpha_manual(values = c(1, 0)) +
        scale_linetype_manual(values = c(1, 5)) +
        coord_fixed() +
        facet_wrap(~ Dataset) +
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
            ),
            shape = guide_legend(
                title          = "Output (features)",
                title.position = "top",
                ncol           = 2,
                byrow          = TRUE,
                order          = 20
            ),
            alpha = guide_legend(
                title          = "Scaling",
                title.position = "top",
                ncol           = 2,
                order          = 30
            ),
            size = guide_legend(
                title          = "Overall score",
                title.position = "top",
                ncol           = 1,
                order          = 40
            ),
            linetype = guide_legend(
                title          = "",
                title.position = "top",
                ncol           = 1,
                override.aes   = list(colour = c("blue", "red")),
                order          = 90
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

#' Load scores
#'
#' Read a set of scores files and return a single tibble.
#'
#' @param scores_files Vector of paths to CSV files containing scores for each
#' dataset. These files are produced by the `plotSingleAtlas()` function. If
#' they are store in a directory named `data/` this vector can be created using
#' `fs::dir_ls("data")`.
#'
#' @return A tibble
#'
#' @author Luke Zappia
loadScores <- function(scores_files) {

    `%>%` <- magrittr::`%>%`

    dataset_key <- c(
        pancreas                       = "Pancreas",
        lung_atlas                     = "Lung",
        immune_cell_hum                = "Immune (human)",
        immune_cell_hum_mou            = "Immune (human/mouse)",
        mouse_brain                    = "Mouse brain",
        simulations_1_1                = "Sim 1",
        simulations_2                  = "Sim 2",
        mouse_brain_atac_genes_large   = "ATAC large (genes)",
        mouse_brain_atac_peaks_large   = "ATAC large (peaks)",
        mouse_brain_atac_windows_large = "ATAC large (windows)",
        mouse_brain_atac_genes_small   = "ATAC small (genes)",
        mouse_brain_atac_peaks_small   = "ATAC small (peaks)",
        mouse_brain_atac_windows_small = "ATAC small (windows)"
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

    scores <- scores_list %>%
        dplyr::bind_rows() %>%
        dplyr::mutate(
            Dataset = factor(
                Dataset,
                levels = names(dataset_key),
                labels = dataset_key
            ),
            OutputFeatures = paste(Output, Features, sep = "_")
        ) %>%
        dplyr::select(Dataset, Datatype, Method, Output, Features,
                      OutputFeatures, Scaling, `Overall Score`,
                      `Batch Correction`, `Bio conservation`, -X1,
                      dplyr::everything()) %>%
        dplyr::filter(!is.na(`Overall Score`))

    return(scores)
}
