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
# - 'csv_file_path' would be the full path of the csv file (or not if you have it in the same folder 
# 'n_metrics_batch' is by default 4 (PCR batch, Batch ASW, iLISI, kBET)
# 'n_metrics_bio' is by default 7 (NMI cluster/label, ARI cluster/label, Cell type ASW, isolated label F1, isolated label sihlouette, CC conservation, cLISI)
# in case the csv file does not contain one or more metrics, just modify the number corresponding to the right group


plotSingleAtlas <- function(csv_file_path){
  
  metrics_tab_lab <- read.csv(csv_file_path, sep = ",")

  # get metrics names from columns
  metrics <- colnames(metrics_tab_lab)[-1]
  metrics <- gsub("\\.", "/", metrics)
  metrics <- gsub("_", " ", metrics)
  metrics <- plyr::mapvalues(metrics, from = c("ASW label", "ASW label/batch", "cell cycle conservation", "hvg overlap", "trajectory"), 
                             to = c("Cell type ASW", "Batch ASW", "CC conservation", "HVG conservation", "trajectory conservation"))
  
  # metrics names as they are supposed to be ordered
  group_batch <- c("PCR batch", "Batch ASW", "iLISI", "kBET")
  group_bio <- c("NMI cluster/label", "ARI cluster/label", "Cell type ASW", 
                 "isolated label F1", "isolated label silhouette", "CC conservation", "HVG conservation", "trajectory conservation", "cLISI")
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
  
  
    
  
  ###### Plot one figure for each data scenario
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
    write.csv(metrics_tab, file = paste0("./", dt.sc, "_summary_scores.csv"), quote = F)
    
    # Defining column_info, row_info and palettes
    row_info <- data.frame(id = metrics_tab$Method)
    
    column_info <- data.frame(id = colnames(metrics_tab),
                              group = c("Text", "Text", "Text", "Text", "Score overall", 
                                        rep("Removal of batch effects", (1 + n_metrics_batch)),
                                        rep("Cell type label variance", (1 + n_metrics_bio))), 
                              geom = c("text", "text", "text", "text", "bar", "bar", 
                                       rep("circle", n_metrics_batch), "bar", rep("circle", n_metrics_bio)),
                              width = c(3.5,2.5,2,2.5,2,2, rep(1,n_metrics_batch), 2, rep(1,n_metrics_bio)),
                              overlay = F)
    
    # defining colors palette
    palette.score.all <- colorRampPalette(rev(brewer.pal(9, "YlGnBu")))(nrow(metrics_tab))
    palette.score.batch <- colorRampPalette(rev(brewer.pal(9, "BuPu")))(nrow(metrics_tab))
    palette.score.celltype <- colorRampPalette(rev(brewer.pal(9, "RdPu")))(nrow(metrics_tab))

    
    palettes <- list("Score overall" = palette.score.all,
                     "Removal of batch effects" = palette.score.batch,
                     "Cell type label variance" = palette.score.celltype)
    
    
    g <- scIB_knit_table(data = metrics_tab, column_info = column_info, row_info = row_info, palettes = palettes, usability = F)  
    ggsave(paste0(dt.sc, "_summary_metrics.pdf"), g, device = cairo_pdf, width = g$width/4, height = g$height/4)
    
    
    
  }
  
  

}