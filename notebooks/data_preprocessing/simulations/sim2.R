library(splatter)

source("sample.R")

base_params <- newSplatParams()

sim2_params <- setParams(
    base_params,
    lib.loc        = 12,
    # Two batches with 1000 cells each. Needs to be more than we want in the
    # final simulation.
    batchCells     = c(1000, 1000, 1000, 1000) * 10,
    batch.facLoc   = c(0.20, 0.25, 0.22, 0.28),
    batch.facScale = c(0.10, 0.08, 0.12, 0.10),
    # Groups with equal probabilities
    group.prob     = rep(1, 4) / 4,
    # Differential expression by group
    de.prob        = c(0.10, 0.12, 0.08, 0.20),
    de.facLoc      = c(0.10, 0.08, 0.12, 0.18),
    de.facScale    = c(0.40, 0.30, 0.45, 0.48),
    # Seed
    seed           = 1
)

# Simulate the full dataset that we will downsample
sim2_full <- splatSimulateGroups(sim2_params)

# Number of cells in each batch in the final simulation
batch_ncells <- c(500, 500, 500, 500) * 10

# Proportions in each group in each batch that we want. Should sum to 1.
batch_props <- list(
    c(0.35, 0.22, 0.25, 0.18),
    c(0.25, 0.40, 0.27, 0.10),
    c(0.30, 0.15, 0.35, 0.20),
    c(0.10, 0.30, 0.20, 0.40)
)

# Downsample cells in our simulation
message("Downsampling cells...")
sim2 <- sample_group_props(sim2_full, batch_ncells, batch_props)

# Set proportion of counts for each batch (relative library size)
batch_count_props <- c(
    Batch1 = 1.00,
    Batch2 = 0.70,
    Batch3 = 0.60,
    Batch4 = 0.30
)
cell_count_props <- batch_count_props[colData(sim2)$Batch]

add_subbatches <- function(sim, subbatch_props) {
    message("Adding subbatches...")

    params <- metadata(sim)$Params
    n_batches <- getParam(params, "nBatches")

    if (!all(sapply(subbatch_props, sum) == 1)) {
        stop("Not all subbatch props sum to 1")
    }

    colData(sim)$SubBatch <- NA

    for (batch_idx in seq_len(n_batches)) {

        message("Processing Batch ", batch_idx, "...")

        batch <- paste0("Batch", batch_idx)
        is_batch <- colData(sim)$Batch == batch
        batch_size <- sum(is_batch)

        subbatches <- BBmisc::chunk(colnames(sim)[is_batch],
                                    props = subbatch_props[[batch_idx]],
                                    shuffle = FALSE)

        for (subbatch_idx in seq_along(subbatches)) {

            message("Adding SubBatch ", subbatch_idx, "...")

            subbatch_cells <- subbatches[[subbatch_idx]]

            noise_sim <- splatSimulate(
                batchCells = length(subbatch_cells),
                lib.loc    = 8.00,
                lib.scale  = 0.50,
                verbose    = FALSE,
                seed       = batch_idx * subbatch_idx
            )

            counts(sim)[, subbatch_cells] <- counts(sim)[, subbatch_cells] +
                counts(noise_sim)
            colData(sim)[subbatch_cells, "Sub"] <- paste0("Sub", subbatch_idx)
            colData(sim)[subbatch_cells, "SubBatch"] <- paste0(batch,
                                                               "Sub", subbatch_idx)
        }
    }

    return(sim)
}

subbatch_props <- list(
    c(0.25, 0.25, 0.25, 0.25),
    c(0.30, 0.30, 0.20, 0.20),
    c(0.40, 0.25, 0.25, 0.10),
    c(0.35, 0.35, 0.15, 0.15)
)

sim2 <- add_subbatches(sim2, subbatch_props)

# Downsample counts
message("Downsampling counts...")
counts(sim2) <- DropletUtils::downsampleMatrix(counts(sim2), cell_count_props,
                                               bycol = TRUE)

message("Calculating QC...")
sim2 <- scater::addPerCellQC(sim2)

# Reorder by SubBatch
subbatch_order <- order(colData(sim2)$SubBatch)
sim2 <- sim2[, subbatch_order]

discard <- lapply(unique(colData(sim2)$SubBatch), function(subbatch) {
    in_batch <- colData(sim2)$SubBatch == subbatch
    scater::quickPerCellQC(colData(sim2)[in_batch, ], nmads = 2)$discard
})
discard <- unlist(discard)
colData(sim2)$Discard <- discard

message("Filtering cells...")
sim2_qc <- sim2[, !discard]
message("Filtering genes...")
sim2_qc <- scater::addPerFeatureQC(sim2_qc)
is_exprs <- rowData(sim2_qc)$detected >= 0.01
sim2_qc <- sim2_qc[is_exprs, ]

message("Normalising...")
sim2 <- scater::logNormCounts(sim2)
sim2_qc <- scater::logNormCounts(sim2_qc)
message("Embedding...")
sim2 <- scater::runTSNE(sim2)
sim2_qc <- scater::runTSNE(sim2_qc)

message("Saving QC plots..")
fs::dir_create("qc_plots/sim2")

scater::plotTSNE(sim2, colour_by = "Batch")
ggplot2::ggsave("qc_plots/sim2/tsne-batch.png")
scater::plotTSNE(sim2, colour_by = "SubBatch")
ggplot2::ggsave("qc_plots/sim2/tsne-subbatch.png")
scater::plotTSNE(sim2, colour_by = "Group")
ggplot2::ggsave("qc_plots/sim2/tsne-group.png")
scater::plotTSNE(sim2, colour_by = "total")
ggplot2::ggsave("qc_plots/sim2/tsne-total.png")
scater::plotTSNE(sim2_qc, colour_by = "Batch")
ggplot2::ggsave("qc_plots/sim2/tsne-batch-qc.png")
scater::plotTSNE(sim2_qc, colour_by = "SubBatch", shape_by = "Batch")
ggplot2::ggsave("qc_plots/sim2/tsne-subbatch-qc.png")
scater::plotTSNE(sim2_qc, colour_by = "Group")
ggplot2::ggsave("qc_plots/sim2/tsne-group-qc.png")
scater::plotTSNE(sim2_qc, colour_by = "total")
ggplot2::ggsave("qc_plots/sim2/tsne-total-qc.png")
scater::plotTSNE(sim2, colour_by = "Discard")
ggplot2::ggsave("qc_plots/sim2/tsne-filt.png")
scater::plotColData(sim2, "sum", "detected", colour_by = "Discard",
                    other_fields = c("Batch", "Sub")) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_grid(Batch ~ Sub, scales = "free")
ggplot2::ggsave("qc_plots/sim2/filtering.png")

message("Converting to AnnData...")
anndata <- reticulate::import("anndata")
sim2_adata <- anndata$AnnData(
    X = t(counts(sim2)),
    obs = data.frame(colData(sim2)),
    var = data.frame(rowData(sim2))
)
sim2_qc_adata <- anndata$AnnData(
    X = t(counts(sim2_qc)),
    obs = data.frame(colData(sim2_qc)),
    var = data.frame(rowData(sim2_qc))
)

message("Saving simulation...")
fs::dir_create("simulations")
# Remove intermediate matrices to reduce file size
assays(sim2) <- assays(sim2)[assayNames(sim2) == "counts"]
assays(sim2_qc) <- assays(sim2_qc)[assayNames(sim2_qc) == "counts"]
saveRDS(sim2, "simulations/sim2.Rds")
saveRDS(sim2_qc, "simulations/sim2_qc.Rds")
sim2_adata$write(filename = "simulations/sim2.h5ad")
sim2_qc_adata$write(filename = "simulations/sim2_qc.h5ad")
