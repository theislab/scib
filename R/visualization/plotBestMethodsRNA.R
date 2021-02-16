#' Plot scib best methods - RNA
#' @description Plotting scib-metrics result of the integration across multiple RNA tasks in an
#'  overview plot that shows the best-performing combination of pre-processing choices for 
#'  each integration method.
#'
#' @return \code{plotBestMethodsRNA} saves in `outdir` the overview table in three formats (.pdf/.tiff/.png)
#'
#' @param csv_metrics_path path to a .csv file output of scib pipeline that contains the metrics calculated 
#' across multiple tasks. All tasks will be considered for ranking best methods, except simulations.
#' @param outdir output directory where the plots will be saved.
#' @param csv_usability_path path to a .csv file containing the results of the usability analysis. 
#' Default to "/data/usability4bestMethods.csv". 
#' These scores will NOT be used for ranking best methods.
#' @param csv_scalability_time_path path to a .csv file containing the results of the scalability analysis, 
#' regarding run time. Default to "/data/scalability_score_time.csv". 
#' These scores will NOT be used for ranking best methods.
#' @param csv_scalability_memory_path path to a .csv file containing the results of the scalability analysis, 
#' regarding memory consumption. Default to "/data/scalability_score_memory.csv". 
#' These scores will NOT be used for ranking best methods.
#' @param ids_RNA character vector of ids for RNA tasks, as they are named in `csv_metrics_path`.
#' @param ids_simulation character vector of ids for simulated tasks, as they are named in `csv_metrics_path`.
#' @param labels_RNA character vector of label names for RNA tasks, to rename ids.These names will be plotted in the summary table.
#' @param labels_simulation character vector of label names for simulated tasks, to rename ids.These names will be plotted in the summary table.
#' @param weight_batch number in [0,1] to use as weight for the batch correction metrics. Weight for
#' bio conservation is calculated as 1-weight_batch
#' 
#' @example 
#' plotBestMethodsRNA(csv_metrics_path = "./data/metrics_RNA_allTasks.csv", 
#' ids_RNA = c("pancreas", "lung_atlas", "immune_cell_hum", "immune_cell_hum_mou", "mouse_brain"),
#' ids_simulation = c("simulations_1_1", "simulations_2"),
#' labels_RNA = c("Pancreas", "Lung", "Immune (hum)", "Immune (hum & mou)", "Brain (mou)"),
#' labels_simulation = c("Sim 1", "Sim 2"))
#' 
#' 
#' 
library(tibble)
library(RColorBrewer)
library(dynutils)
library(stringr)
library(Hmisc)
library(plyr)

source("knit_table.R") # Please put knit_table.R in your working dir

plotBestMethodsRNA <- function(csv_metrics_path,
                               outdir = ".",
                               csv_usability_path = "./data/usability4bestMethods.csv", 
                               csv_scalability_time_path = "./data/scalability_score_time.csv", 
                               csv_scalability_memory_path = "./data/scalability_score_memory.csv", 
                               ids_RNA, ids_simulation,
                               labels_RNA, labels_simulation,
                               weight_batch = 0.4){
  
  metrics_tab_lab <- read.csv(csv_metrics_path, sep = ",")
  n_atlas_RNA <- length(ids_RNA)
  n_simulation <- length(ids_simulation)
  
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
  # order scenarios 
  data.scenarios <- data.scenarios[match(c(ids_RNA, ids_simulation), data.scenarios)]

  methods.table.list <- list()
  
  ###### Get overall score for each data scenario
  for (dt.sc in data.scenarios){
    ind.scen <- grep(paste0(dt.sc, "/"), methods_info_full)
    methods_info <- methods_info_full[ind.scen]
    metrics_tab_sub <- metrics_tab_lab[ind.scen, ]
    
    # info on scaling data
    scaling <- sapply(str_split(methods_info, "/"), function(x) x[3])
    
    
    # info on HVG selection
    hvg <- sapply(str_split(methods_info, "/"), function(x) x[4])
    hvg <- plyr::mapvalues(hvg, from = c("hvg", "full_feature"), to = c("HVG", "FULL"))
    
    methods <- sapply(str_split(methods_info, "/"), function(x) x[5])
    
    methods_name <- sapply(str_split(methods, "_"), function(x) x[1])
    methods_name <- capitalize(methods_name)
    methods_name <- plyr::mapvalues(methods_name, 
                                    from = c("Seurat", "Seuratrpca", "Mnn", "Bbknn", "Trvae", "Scvi", "Liger", "Combat", "Saucie", "Fastmnn", "Desc", "Scanvi", "Scgen"), 
                                    to = c("Seurat v3 CCA", "Seurat v3 RPCA", "MNN", "BBKNN", "trVAE", "scVI", "LIGER", "ComBat", "SAUCIE", "fastMNN", "DESC", "scANVI", "scGen"))
    
    
    method_groups <- sapply(str_split(methods, "_"), function(x) x[2])
    method_groups <- plyr::mapvalues(method_groups, 
                                     from = c("knn", "embed", "full"), 
                                     to = c("graph", "embed", "gene"))
    
    
    
    ##### Create dataframe 
    metrics_tab <- as.data.frame(metrics_tab_sub[, -1])
    metrics_tab[metrics_tab == ""] <- NA
    colnames(metrics_tab) <- metrics
    
    #add Methods column
    metrics_tab <- add_column(metrics_tab, "Method" = methods_name, .before = 1)
    
    # reorder columns by metrics
    col.ordered <- c("Method", metrics.ord)
    metrics_tab <- metrics_tab[, col.ordered]
    
    #####----- Set cLISI scores to NA for BBKNN
    #metrics_tab[metrics_tab$Method == "BBKNN", "cLISI"] <- NA
    
    
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
    metrics_tab <- add_column(metrics_tab, "Bio conservation" = score_group_bio, .after = 3+n_metrics_batch)
    
    metrics_tab <- add_column(metrics_tab, "Output" = method_groups, .after = "Method")
    metrics_tab <- add_column(metrics_tab, "Features" = hvg, .after = "Output")
    metrics_tab <- add_column(metrics_tab, "Scaling" = scaling, .after = "Features")
    
    # order methods by the overall score
    metrics_tab <- metrics_tab[order(metrics_tab$`Overall Score`,  decreasing = T), ]
    
    metrics_tab <- metrics_tab[, c("Method", "Output", "Features", "Scaling", "Overall Score")]
    colnames(metrics_tab)[5] <- dt.sc
    methods.table.list[[dt.sc]] <- metrics_tab
    
  }
  
  methods.table.merged <- merge(methods.table.list[[1]], methods.table.list[[2]], by = c("Method", "Output", "Features", "Scaling"),
                                all = T)
  for(i in 3:length(methods.table.list)){
    methods.table.merged <- merge(methods.table.merged, methods.table.list[[i]], by = c("Method", "Output", "Features", "Scaling"),
                                  all = T)
  }
  
  # Rename columns
  colnames(methods.table.merged) <- plyr::mapvalues(colnames(methods.table.merged), 
                                                    from = c(ids_RNA, ids_simulation),
                                                    to = c(labels_RNA, labels_simulation))
  
  # Assign unintegrated overall score to methods that did not run over one/more atlases
  
  atlas.ranks <- methods.table.merged
  # columns to be ranked on
  rank.cols <- labels_RNA
  
  for(c in rank.cols){
    na.idx <- is.na(atlas.ranks[,c])
    if(sum(na.idx)>0){
      atlas.ranks[na.idx,c] <- atlas.ranks[atlas.ranks$Method == "Unintegrated", c]
    }
  }
 
  
 
  atlas.ranks[, rank.cols] <- apply(atlas.ranks[, rank.cols], 2, function(x) rank(-x, na.last = T, ties.method = "average"))
  avg.ranks <- apply(atlas.ranks[, rank.cols], 1, mean)
  
  
  # order atlas.rank by average rank
  atlas.rank.ord <- atlas.ranks[order(avg.ranks, decreasing = F), ]
  
  # Keep best performing solution for each method
  keep.best <- NULL
  for(met in unique(atlas.ranks$Method)){
    if(met %in% c("Scanorama", "fastMNN", "SAUCIE")){
      keep.best <- c(keep.best, which(atlas.rank.ord$Method == met & atlas.rank.ord$Output == "gene")[1])
      keep.best <- c(keep.best, which(atlas.rank.ord$Method == met & atlas.rank.ord$Output == "embed")[1])
    } else{
      keep.best <- c(keep.best, which(atlas.rank.ord$Method == met)[1])
    }
    
  }
  
  atlas.rank.ord <- atlas.rank.ord[sort(keep.best),]
  atlas.rank.ord[, rank.cols] <- apply(atlas.rank.ord[, rank.cols], 2, function(x) rank(x, na.last = T, ties.method = "average"))
  avg.ranks2 <- apply(atlas.rank.ord[, rank.cols], 1, mean)
  
  
  # order atlas.rank by average rank
  best_methods_tab <- atlas.rank.ord[order(avg.ranks2, decreasing = F), ]
  
  
  
  
  best_methods_tab <- merge(best_methods_tab[, 1:4], methods.table.merged, by = c("Method", "Output", "Features", "Scaling"),
                            all = F, sort = F)
  
  
  
  # Delete rows that are empty
  rowsNA <- which(rowSums(is.na(best_methods_tab[, rank.cols])) == (n_atlas_RNA+n_simulation))
  if(length(rowsNA) >0){
    best_methods_tab <- best_methods_tab[-rowsNA, ]
    }
  
  
  
  
  ########## ADD USABILITY TABLE
  usability_mat <- read.csv(csv_usability_path)
  usability_mat <- usability_mat %>%
    filter(!(methods %in% c("ComBat", "MNN")))
  usability_mat$methods <- mapvalues(usability_mat$methods, from = c("ComBat (Scanpy)", "MNNpy"), 
                                     to = c("ComBat", "MNN"))
  methods_best <- mapvalues(best_methods_tab$Method, from= c("Seurat v3 RPCA", "Seurat v3 CCA"), 
                            to = c("Seurat v3", "Seurat v3"))
  best_methods_tab$`Usability Package` <- usability_mat[match(tolower(methods_best), 
                                                          tolower(usability_mat$methods)), "Package"]
  best_methods_tab$`Usability Paper` <- usability_mat[match(tolower(methods_best), 
                                                          tolower(usability_mat$methods)), "Paper"]
  
  ######### ADD SCALABILITY time and memory
  scalability_time <- read.csv(csv_scalability_time_path)
  scalability_memory <- read.csv(csv_scalability_memory_path)
  
  # create comparable strings out of best_methods
  hvg_full <- mapvalues(best_methods_tab$Features, from = c("HVG", "FULL"), to = c("hvg", "full_feature"))
  methods_string <- tolower(mapvalues(best_methods_tab$Method, from = c("Seurat v3 RPCA", "Seurat v3 CCA"),
                                      to = c("seuratrpca", "seurat")))
  best_method_string <- paste0(best_methods_tab$Scaling, "/", hvg_full, "/", methods_string)
  
  # match best methods with scalability data
  best_methods_tab$`Scalability time` <- scalability_time[match(best_method_string, scalability_time$metrics), "AUC_scaled"]
  best_methods_tab$`Scalability memory` <- scalability_memory[match(best_method_string, scalability_memory$metrics), "AUC_scaled"]
  
  ############## add first column = ranking
  best_methods_tab <- add_column(best_methods_tab, "Ranking" = 1:nrow(best_methods_tab), .before = "Method")
  
  # Adjust scGen and scANVI names
  best_methods_tab$Method <- mapvalues(best_methods_tab$Method, from = c("scGen", "scANVI"), to = c("scGen*", "scANVI*"))
  
  # Defining column_info, row_info and palettes
  row_info <- data.frame(id = best_methods_tab$Method)
  
  column_info <- data.frame(id = colnames(best_methods_tab),
                            group = c("Text", "Text", "Image", "Text", "Text",  
                                      rep("RNA", n_atlas_RNA),
                                      rep("Simulation", n_simulation), 
                                      "Usability","Usability", "Scalability", "Scalability"), 
                            geom = c("text", "text", "image", "text", "text", 
                                     rep("bar", n_atlas_RNA + n_simulation + 4)),
                            width = c(1.5,8,2.5,2,2.5, rep(2,n_atlas_RNA + n_simulation + 4)),
                            overlay = F)
  
  # defining colors palette
  palettes <- list("RNA" = "Blues",
                   "Simulation" = "Greens",
                   "Usability" = "Oranges",
                   "Scalability" = "Greys")
 
  
  g <- scIB_knit_table(data = best_methods_tab, column_info = column_info, row_info = row_info, palettes = palettes, usability = T)
  now <- Sys.time()
  ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), "BestMethods_summary.pdf"), g, device = cairo_pdf, width = 210, height = 297, units = "mm")
  ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), "BestMethods_summary.tiff"), g, device = "tiff", dpi = "retina", width = 210, height = 297, units = "mm")
  ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), "BestMethods_summary.png"), g, device = "png", dpi = "retina", width = 210, height = 297, units = "mm")
  
}
