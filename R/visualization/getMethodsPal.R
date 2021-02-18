#' Get methods paletter
#'
#' Example of a vector specifying colour paletter for methods used by
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
