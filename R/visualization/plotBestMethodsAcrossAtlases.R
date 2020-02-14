# library(ggplot2)
# library(cowplot)
# library(fmsb)
# 
# library(tidyr)
# library(ggforce)
library(tibble)
library(RColorBrewer)
library(dynutils)
library(stringr)
library(Hmisc)
library(plyr)
# library(R.utils)

source("/home/python_scRNA/Munich/visualization/knit_table.R")# You will need to have in the same folder knit_table.R and this plotSingleAtlas.R

# parameters: 
# - 'csv_atlases_path' would be the full path of the csv file with metrics over simulation + real atlases 
# - 'csv_usability_path' would be the path of the usability sheet (which does not change)
# - 'csv_scalability_path' would be the path of the scalability file

plotBestMethodsAcrossAtlases <- function(csv_atlases_path, 
                                         csv_usability_path = "./scIB Usability  - Sheet4.csv", 
                                         csv_scalability_path = "./scalability_tab.csv", 
                                         n_atlas_RNA = 3, n_simulation = 1){
  
  #csv_atlases_path <- "./test_best.csv"
  
  metrics_tab_lab <- read.csv(csv_atlases_path, sep = ",")
  
  # get metrics names from columns
  metrics <- colnames(metrics_tab_lab)[-1]
  metrics <- gsub("\\.", "/", metrics)
  metrics <- gsub("_", " ", metrics)
  metrics <- plyr::mapvalues(metrics, from = c("ASW label", "ASW label/batch", "cell cycle conservation"), 
                             to = c("Cell type ASW", "Batch ASW", "CC conservation"))
  
  # metrics names as they are supposed to be ordered
  group_batch <- c("PCR batch", "Batch ASW", "iLISI", "kBET")
  group_bio <- c("NMI cluster/label", "ARI cluster/label", "Cell type ASW", 
                 "isolated label F1", "isolated label silhouette", "CC conservation","cLISI")
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
  
  # data scenarios to be saved in file name
  data.scenarios <- unique(unlist(sapply(str_split(methods_info_full, "/"), function(x) x[1])))
  
  
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
                                    from = c("Seurat", "Mnn", "Bbknn", "Trvae", "Scvi", "Liger", "Combat"), 
                                    to = c("Seurat v3", "MNN", "BBKNN", "TrVAE", "scVI", "LIGER", "ComBat"))
    
    
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
    
    score_all <- (score_group1 + score_group2) / 2
    
    metrics_tab <- add_column(metrics_tab, "Overall Score" = score_all, .after = "Method")
    metrics_tab <- add_column(metrics_tab, "Batch Correction" = score_group1, .after = "Overall Score")
    metrics_tab <- add_column(metrics_tab, "Bio conservation" = score_group2, .after = "kBET")
    
    metrics_tab <- add_column(metrics_tab, "Output" = method_groups, .after = "Method")
    metrics_tab <- add_column(metrics_tab, "Features" = hvg, .after = "Output")
    metrics_tab <- add_column(metrics_tab, "Scaling" = scaling, .after = "Features")
    
    # order methods by the overall score
    metrics_tab <- metrics_tab[order(metrics_tab$`Overall Score`,  decreasing = T), ]
    
    metrics_tab <- metrics_tab[, c("Method", "Output", "Features", "Scaling", "Overall Score")]
    colnames(metrics_tab)[5] <- paste("Score", dt.sc)
    methods.table.list[[dt.sc]] <- metrics_tab
    
  }
  
  methods.table.merged <- merge(methods.table.list[[1]], methods.table.list[[2]], by = c("Method", "Output", "Features", "Scaling"),
                                all = T)
  for(i in 3:length(methods.table.list)){
    methods.table.merged <- merge(methods.table.merged, methods.table.list[[i]], by = c("Method", "Output", "Features", "Scaling"),
                                  all = T)
  }
  atlas.ranks <- methods.table.merged
  atlas.ranks[, 5:ncol(atlas.ranks)] <- apply(atlas.ranks[, 5:ncol(atlas.ranks)], 2, function(x) rank(-x, na.last = T, ties.method = "average"))
  avg.ranks <- apply(atlas.ranks[, 5:ncol(atlas.ranks)], 1, mean)
  
  
  # order atlas.rank by average rank
  atlas.rank.ord <- atlas.ranks[order(avg.ranks, decreasing = F), ]
  
  # Keep best performing solution for each method
  keep.best <- NULL
  for(met in unique(atlas.ranks$Method)){
    keep.best <- c(keep.best, which(atlas.rank.ord$Method == met)[1])
  }
  
  best_methods_tab <- atlas.rank.ord[sort(keep.best),]
  best_methods_tab <- merge(best_methods_tab[, 1:4], methods.table.merged, by = c("Method", "Output", "Features", "Scaling"),
                            all = F, sort = F)
  


  ########## ADD USABILITY TABLE
  usability_mat <- read.csv(csv_usability_path)
  avg.usability <- rowMeans(usability_mat[,-1])
  
  best_methods_tab$Usability <- avg.usability[match(best_methods_tab$Method, usability_mat$Method)]
  
  ######### ADD SCALABILITY
  
  scalability_mat <- read.csv(csv_scalability_path)
  best_methods_tab$Scalability <- scalability_mat[match(best_methods_tab$Method, scalability_mat$Method), "Scalability"]
  
  # Defining column_info, row_info and palettes
  row_info <- data.frame(id = best_methods_tab$Method)
  
  column_info <- data.frame(id = colnames(best_methods_tab),
                            group = c("Text", "Text", "Text", "Text",  
                                      rep("RNA", n_atlas_RNA),
                                      rep("Simulation", n_simulation), 
                                      "Usability", "Scalability"), 
                            geom = c("text", "text", "text", "text", 
                                     rep("bar", n_atlas_RNA + n_simulation + 2)),
                            width = c(3.5,2.5,2,2.5, rep(2,n_atlas_RNA + n_simulation + 2)),
                            overlay = F)
  
  # defining colors palette
  palettes <- list("RNA" = colorRampPalette(rev(brewer.pal(9, "Blues")))(nrow(best_methods_tab)),
                   "Simulation" = colorRampPalette(rev(brewer.pal(9, "Greens")))(nrow(best_methods_tab)),
                   "Usability" = colorRampPalette(rev(brewer.pal(9, "Oranges")))(nrow(best_methods_tab)),
                   "Scalability" = colorRampPalette(rev(brewer.pal(9, "Greys")))(nrow(best_methods_tab)))
  
  
  g <- scIB_knit_table(data = best_methods_tab, column_info = column_info, row_info = row_info, palettes = palettes, usability = T)  
  ggsave("BestMethods_summary_metrics.pdf", g, device = cairo_pdf, width = g$width/4, height = g$height/4)

}