library(tibble)
library(RColorBrewer)
library(dynutils)
library(stringr)
library(Hmisc)
library(plyr)

source("/home/python_scRNA/Munich/visualization/knit_table.R")# You will need to have in the same folder knit_table.R and this plotSingleAtlas.R


plotBestMethodsATAC<- function(csv_atac_small_path, csv_atac_large_path){
  
  atac_small <- read.csv(csv_atac_small_path, sep = ",")
  atac_large <- read.csv(csv_atac_large_path, sep = ",")
  
  metrics_tab_lab <- rbind(atac_small, atac_large)
    
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
  if(substring(methods_info_full[1], 1, 1) == "/"){
    methods_info_full <- sub("/", "", methods_info_full)
  }
  
  # Remove trvae full
  ind.trvae_full <- grep("trvae_full", methods_info_full)
  if(length(ind.trvae_full) >0){
  methods_info_full <- methods_info_full[-ind.trvae_full]
  metrics_tab_lab <- metrics_tab_lab[-ind.trvae_full,]
  }
  # data scenarios to be saved in file name
  data.scenarios <- unique(unlist(sapply(str_split(methods_info_full, "/"), function(x) x[1])))
  
  
  methods.table.list <- list()
  
  ###### Get overall score for each data scenario
  for (dt.sc in data.scenarios){
    ind.scen <- grep(paste0(dt.sc, "/"), methods_info_full)
    methods_info <- methods_info_full[ind.scen]
    metrics_tab_sub <- metrics_tab_lab[ind.scen, ]
    
    methods <- sapply(str_split(methods_info, "/"), function(x) x[5])
    
    methods_name <- sapply(str_split(methods, "_"), function(x) x[1])
    methods_name <- capitalize(methods_name)
    methods_name <- plyr::mapvalues(methods_name, 
                                    from = c("Seurat", "Mnn", "Bbknn", "Trvae", "Scvi", "Liger", "Combat", "Saucie", "Fastmnn", "Desc", "Scanvi", "Scgen"), 
                                    to = c("Seurat v3", "MNN", "BBKNN", "trVAE", "scVI", "LIGER", "ComBat", "SAUCIE", "fastMNN", "DESC", "scANVI", "scGen"))
    
    
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
    score_group1 <- rowMeans(scaled_metrics_tab[, 1:n_metrics_batch], na.rm = T)
    score_group2 <- rowMeans(scaled_metrics_tab[, (1+n_metrics_batch):ncol(scaled_metrics_tab)], 
                             na.rm = T)
    
    score_all <- (0.4*score_group1 + 0.6*score_group2)
    
    metrics_tab <- add_column(metrics_tab, "Overall Score" = score_all, .after = "Method")
    metrics_tab <- add_column(metrics_tab, "Batch Correction" = score_group1, .after = "Overall Score")
    metrics_tab <- add_column(metrics_tab, "Bio conservation" = score_group2, .after = "kBET")
    
    metrics_tab <- add_column(metrics_tab, "Output" = method_groups, .after = "Method")
    
    
    # order methods by the overall score
    metrics_tab <- metrics_tab[order(metrics_tab$`Overall Score`,  decreasing = T), ]
    
    metrics_tab <- metrics_tab[, c("Method", "Output", "Overall Score")]
    colnames(metrics_tab)[3] <- dt.sc
    methods.table.list[[dt.sc]] <- metrics_tab
    
  }
  
  methods.table.merged <- merge(methods.table.list[[1]], methods.table.list[[2]], by = c("Method", "Output"),
                                all = T)
  
  
  # Rename columns
  colnames(methods.table.merged) <- plyr::mapvalues(colnames(methods.table.merged), 
                                                    from = c("mouse_brain_atac_large_3datasets", "mouse_brain_atac_small_3datasets"),
                                                    to = c("Brain (mou) large", "Brain (mou) small"))
  # columns to be ranked on
  rank.cols <- c("Brain (mou) large", "Brain (mou) small")
  
  # Assign unintegrated overall score to methods that did not run over one/more atlases
  atlas.ranks <- methods.table.merged
  
  for(c in rank.cols){
    na.idx <- is.na(atlas.ranks[,c])
    if(sum(na.idx)>0){
      atlas.ranks[na.idx,c] <- atlas.ranks[atlas.ranks$Method == "Unintegrated", c]
    }
  }
  
  
 
  atlas.ranks[, rank.cols] <- apply(atlas.ranks[, rank.cols], 2, function(x) rank(-x, ties.method = "average"))
  avg.ranks <- apply(atlas.ranks[, rank.cols], 1, mean)
  
  
  # order atlas.rank by average rank
  atlas.rank.ord <- atlas.ranks[order(avg.ranks, decreasing = F), ]
  
  # Keep best performing solution for each method
  keep.best <- NULL
  for(met in unique(atlas.ranks$Method)){
    if(met == "Scanorama"){
      keep.best <- c(keep.best, which(atlas.rank.ord$Method == met & atlas.rank.ord$Output == "gene")[1])
      keep.best <- c(keep.best, which(atlas.rank.ord$Method == met & atlas.rank.ord$Output == "embed")[1])
    } else{
      keep.best <- c(keep.best, which(atlas.rank.ord$Method == met)[1])
    }
    
  }
  
  best_methods_tab <- atlas.rank.ord[sort(keep.best),]
  best_methods_tab <- merge(best_methods_tab[, 1:2], methods.table.merged, by = c("Method", "Output"),
                            all = F, sort = F)
  
  # Delete rows that are empty
  rowsNA <- which(rowSums(is.na(best_methods_tab[, rank.cols])) == 2)
  if(length(rowsNA)>0){
    best_methods_tab <- best_methods_tab[-rowsNA, ]
  }

  
  ############## add first column = ranking
  best_methods_tab <- add_column(best_methods_tab, "Ranking" = 1:nrow(best_methods_tab), .before = "Method")
  
  
  
  # Defining column_info, row_info and palettes
  row_info <- data.frame(id = best_methods_tab$Method)
  
  column_info <- data.frame(id = colnames(best_methods_tab),
                            group = c("Text", "Text", "Image", "ATAC", "ATAC"), 
                            geom = c("text", "text", "image", "bar", "bar"),
                            width = c(1.5,4,2.5,2.5,2.5),
                            overlay = F)
  
  # defining colors palette
  palettes <- list("ATAC" = "Purples")
  
  
  g <- scIB_knit_table(data = best_methods_tab, column_info = column_info, row_info = row_info, palettes = palettes, usability = F, atac_best = T)
  now <- Sys.time()
  ggsave(paste0(format(now, "%Y%m%d_%H%M%S_"), "ATAC_BestMethods_summary.pdf"), g, device = cairo_pdf, width = 100, height = 150, units = "mm")
  ggsave(paste0(format(now, "%Y%m%d_%H%M%S_"), "ATAC_BestMethods_summary.tiff"), g, device = "tiff", dpi = "retina", width = 100, height = 150, units = "mm")
  ggsave(paste0(format(now, "%Y%m%d_%H%M%S_"), "ATAC_BestMethods_summary.jpeg"), g, device = "jpeg", dpi = "retina", width = 100, height = 150, units = "mm")
  
}