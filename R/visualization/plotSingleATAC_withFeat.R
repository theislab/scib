#' Plot scib single methods - ATAC with features
#' @description Plotting scib-metrics result of the integration for each single ATAC task in an overview table. Multiple feature 
#' spaces are considered simultaneously (genes/windows/peaks).
#' Integration methods are ranked from top to bottom based on an 'Overall Score', calculated as weighted
#' sum of 'Batch correction' (default w:0.4) and 'Bio conservation' (default w:0.6).
#'
#' @return \code{plotSingleATAC_withFeat} saves in `outdir` one overview table for each task in three formats (.pdf/.tiff/.png)
#' and a .csv file for each task containing the ranked summary table scores. 
#'
#' @param csv_metrics_path path to a .csv file output of scib that contains the metrics calculated 
#' across one or multiple ATAC tasks, over multiple feature spaces. 
#' @param outdir output directory where the plots and .csv will be saved.
#' @param weight_batch number in [0,1] to use as weight for the batch correction metrics. Weight for
#' bio conservation is calculated as 1-weight_batch
#' @example
#' plotSingleATAC_withFeat(csv_metrics_path = "./data/ATAC_metrics/metrics_ATAC_large11.csv")
#' 


library(tibble)
library(RColorBrewer)
library(dynutils)
library(stringr)
library(Hmisc)
library(plyr)

source("knit_table.R") # Please put knit_table.R in your working dir


plotSingleATAC_withFeat <- function(csv_metrics_path, outdir = ".", weight_batch = 0.4){
  
  metrics_tab_lab <- read.csv(csv_metrics_path, sep = ",")
  
  
  # get metrics names from columns
  metrics <- colnames(metrics_tab_lab)[-1]
  metrics <- gsub("\\.", "/", metrics)
  metrics <- gsub("_", " ", metrics)
  metrics <- plyr::mapvalues(metrics, from = c("ASW label", "ASW label/batch", "cell cycle conservation", "hvg overlap", "trajectory", "graph conn"), 
                             to = c("Cell type ASW", "Batch ASW", "CC conservation", "HVG conservation", "trajectory conservation", "graph connectivity"))
  
  
  # metrics names as they are supposed to be ordered
  group_batch <- c("PCR batch", "Batch ASW", "iLISI", "graph connectivity", "kBET")
  group_bio <- c("NMI cluster/label", "ARI cluster/label", "Cell type ASW", 
                 "isolated label F1", "isolated label silhouette", "CC conservation", "HVG conservation", "trajectory conservation","cLISI")
  # set original values of number of metrics
  n_metrics_batch_original <- sum(group_batch %in% metrics)
  n_metrics_bio_original <- sum(group_bio %in% metrics)
  
  # order metrics present in the table
  matching.order <- match(c(group_batch, group_bio), metrics)
  metrics.ord <- metrics[matching.order[!is.na(matching.order)]]
  
  # get methods info from rownames
  methods_info_full  <- as.character(metrics_tab_lab[,1])
  
  # in case methods names start with /
  methods_info_full <- sub("^/", "", methods_info_full)
  
  # Remove trvae full
  ind.trvae_full <- grep("trvae_full", methods_info_full)
  if(length(ind.trvae_full) > 0){
    methods_info_full <- methods_info_full[-ind.trvae_full]
    metrics_tab_lab <- metrics_tab_lab[-ind.trvae_full,]
  }
  
  
  # data scenarios to be saved in file name
  data.scenarios <- unique(unlist(sapply(str_split(methods_info_full, "/"), function(x) x[1])))
  
  # get features from data.scenarios
  features <- unlist(sapply(str_split(methods_info_full, "/"), function(x) x[1]))
  features <- unlist(sapply(str_split(features, "_"), function(x) x[4]))
  
  # rename data.scenarios
  data.scenarios <- unique(gsub(pattern = "_genes|_peaks|_windows", "", data.scenarios))
  
  
  
  
  methods_info <- methods_info_full
  metrics_tab_sub <- metrics_tab_lab
    
  methods <- sapply(str_split(methods_info, "/"), function(x) x[5])
  
  methods_name <- sapply(str_split(methods, "_"), function(x) x[1])
  methods_name <- capitalize(methods_name)
  methods_name <- plyr::mapvalues(methods_name, 
                                  from = c("Seurat", "Seuratrpca", "Mnn", "Bbknn", "Trvae", "Scvi", "Liger", "Combat", "Saucie", "Fastmnn", "Desc", "Scanvi", "Scgen"), 
                                  to = c("Seurat v3 CCA", "Seurat v3 RPCA", "MNN", "BBKNN", "trVAE", "scVI", "LIGER", "ComBat", "SAUCIE", "fastMNN", "DESC", "scANVI*", "scGen*"))
  
  
  method_groups <- sapply(str_split(methods, "_"), function(x) x[2])
  method_groups <- plyr::mapvalues(method_groups, 
                                   from = c("knn", "embed", "full"), 
                                   to = c("graph", "embed", "gene"))
  
  features_groups <- features
  
  ##### Create dataframe 
  metrics_tab <- as.data.frame(metrics_tab_sub[, -1])
  metrics_tab[metrics_tab == ""] <- NA
  colnames(metrics_tab) <- metrics
    
  #add Methods column
  metrics_tab <- add_column(metrics_tab, "Method" = methods_name, .before = 1)
  
  # reorder columns by metrics
  col.ordered <- c("Method", metrics.ord)
  metrics_tab <- metrics_tab[, col.ordered]
  
  ## Remove columns that are full NAs
  na.col <- apply(metrics_tab, 2, function(x) sum(is.na(x)) == nrow(metrics_tab))
  # redefine numbers of metrics per group
  if(sum(colnames(metrics_tab)[na.col] %in%  group_batch) > 0){
    n_metrics_batch <- n_metrics_batch_original - sum(colnames(metrics_tab)[na.col] %in%  group_batch)
  } else {
    n_metrics_batch <- n_metrics_batch_original
  }
  
  if(sum(colnames(metrics_tab)[na.col] %in% group_bio) > 0){
    n_metrics_bio <- n_metrics_bio_original - sum(colnames(metrics_tab)[na.col] %in% group_bio)
  } else{
    n_metrics_bio <- n_metrics_bio_original
  }
  
  metrics_tab <- metrics_tab[, !na.col]
    
    
  ## Scores should be already scaled [0,1] - however, we aim to compute the scores based on the min-max scaled metrics
  scaled_metrics_tab <- as.matrix(metrics_tab[, -1])
  scaled_metrics_tab <- apply(scaled_metrics_tab, 2, function(x) scale_minmax(x))
  
  # calculate average score by group and overall
  score_group_batch <- rowMeans(scaled_metrics_tab[, 1:n_metrics_batch], na.rm = T)
  score_group_bio <- rowMeans(scaled_metrics_tab[, (1+n_metrics_batch):ncol(scaled_metrics_tab)], 
                           na.rm = T)
  
  score_all <- (weight_batch*score_group_batch + (1-weight_batch)*score_group_bio)
    
  metrics_tab <- add_column(metrics_tab, "Overall Score" = score_all, .after = "Method")
  metrics_tab <- add_column(metrics_tab, "Batch Correction" = score_group_batch, .after = "Overall Score")
  metrics_tab <- add_column(metrics_tab, "Bio conservation" = score_group_bio, .after = "kBET")
  
  metrics_tab <- add_column(metrics_tab, "Output" = method_groups, .after = "Method")
  metrics_tab <- add_column(metrics_tab, "Feature Space" = features_groups, .after = "Output")
  
  # order methods by the overall score
  metrics_tab <- metrics_tab[order(metrics_tab$`Overall Score`,  decreasing = T), ]
  write.csv(metrics_tab, file = paste0(outdir, "/", data.scenarios, "_summary_scores.csv"), quote = F)
    
    
  # Delete rows that are empty
  rowsNA <- which(is.na(metrics_tab$`Overall Score`))
  if(length(rowsNA) >0){
    metrics_tab <- metrics_tab[-rowsNA, ]
  }
  
  # Defining column_info, row_info and palettes
  row_info <- data.frame(id = metrics_tab$Method)
  
  column_info <- data.frame(id = colnames(metrics_tab),
                            group = c("Text", "Image", "Text", "Score overall", 
                                      rep("Removal of batch effects", (1 + n_metrics_batch)),
                                      rep("Cell type label variance", (1 + n_metrics_bio))), 
                            geom = c("text", "image","text",  "bar", "bar", 
                                     rep("circle", n_metrics_batch), "bar", rep("circle", n_metrics_bio)),
                            width = c(8,2,3,2,2, rep(1,n_metrics_batch), 2, rep(1,n_metrics_bio)),
                            overlay = F)
  
  # defining colors palette
  palettes <- list("Score overall" = "YlGnBu",
                   "Removal of batch effects" = "BuPu",
                   "Cell type label variance" = "RdPu")
  
  
  g <- scIB_knit_table(data = metrics_tab, column_info = column_info, row_info = row_info, palettes = palettes, usability = F, atac = T) 
  now <- Sys.time()
  ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"),  data.scenarios, "_summary_metrics.pdf"), g, device = cairo_pdf, width = 297, height = 420, units = "mm")
  ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"),  data.scenarios, "_summary_metrics.tiff"), g, device = "tiff", dpi = "retina", width = 297, height = 420, units = "mm")
  ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"),  data.scenarios, "_summary_metrics.png"), g, device = "png", dpi = "retina", width = 297, height = 420, units = "mm")
  
  
    
  
  
  
}
