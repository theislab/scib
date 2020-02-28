library(splatter)

source("sample.R")

base_params <- newSplatParams()

sim1_params <- setParams(
    base_params,
    lib.loc        = 12,
    # Six batches with 10000 cells each. Needs to be more than we want in the
    # final simulation.
    batchCells     = c(1000, 1000, 1000, 1000, 1000, 1000) * 10,
    batch.facLoc   = c(0.10, 0.15, 0.12, 0.18, 0.15, 0.25),
    batch.facScale = c(0.10, 0.08, 0.12, 0.10, 0.12, 0.18),
    # Groups with equal probabilities
    group.prob     = rep(1, 7) / 7,
    # Differential expression by group
    de.prob        = c(0.10, 0.12, 0.08, 0.20, 0.12, 0.10, 0.16),
    de.facLoc      = c(0.10, 0.08, 0.12, 0.18, 0.06, 0.20, 0.14),
    de.facScale    = c(0.40, 0.30, 0.45, 0.48, 0.42, 0.38, 0.36),
    # Seed
    seed           = 1
)

# Simulate the full dataset that we will downsample
sim1_full <- splatSimulateGroups(sim1_params)

# Number of cells in each batch in the final simulation
batch_ncells <- c(300, 250, 220, 200, 180, 100) * 10
# Proportions in each group in each batch that we want. Should sum to 1.
batch_props <- list(
    c(0.35, 0.22, 0.15, 0.13, 0.10, 0.05, 0.00),
    c(0.25, 0.33, 0.20, 0.10, 0.00, 0.08, 0.04),
    c(0.30, 0.30, 0.00, 0.20, 0.10, 0.10, 0.00),
    c(0.15, 0.20, 0.35, 0.15, 0.15, 0.00, 0.00),
    c(0.40, 0.05, 0.10, 0.10, 0.17, 0.12, 0.06),
    c(0.45, 0.35, 0.15, 0.05, 0.00, 0.00, 0.00)
)

# Downsample cells in our simulation
message("Downsampling cells...")
sim1 <- sample_group_props(sim1_full, batch_ncells, batch_props)

# Set proportion of counts for each batch (relative library size)
batch_count_props <- c(
    Batch1 = 0.30,
    Batch2 = 0.25,
    Batch3 = 0.32,
    Batch4 = 0.45,
    Batch5 = 0.28,
    Batch6 = 1.00
)
cell_count_props <- batch_count_props[colData(sim1)$Batch]

# Downsample counts
message("Downsampling counts...")
counts(sim1) <- DropletUtils::downsampleMatrix(counts(sim1), cell_count_props,
                                               bycol = TRUE)

message("Calculating QC...")
sim1 <- scater::addPerCellQC(sim1)

discard <- lapply(unique(colData(sim1)$Batch), function(batch) {
    in_batch <- colData(sim1)$Batch == batch
    scater::quickPerCellQC(colData(sim1)[in_batch, ], nmads = 2)$discard
})
discard <- unlist(discard)
colData(sim1)$Discard <- discard

message("Filtering cells...")
sim1_qc <- sim1[, !discard]
message("Filtering genes...")
sim1_qc <- scater::addPerFeatureQC(sim1_qc)
is_exprs <- rowData(sim1_qc)$detected >= 0.01
sim1_qc <- sim1_qc[is_exprs, ]

message("Normalising...")
sim1 <- scater::logNormCounts(sim1)
sim1_qc <- scater::logNormCounts(sim1_qc)
message("Embedding...")
sim1 <- scater::runTSNE(sim1)
sim1_qc <- scater::runTSNE(sim1_qc)

message("Saving QC plots..")
fs::dir_create("qc_plots/sim1")

scater::plotTSNE(sim1, colour_by = "Batch")
ggplot2::ggsave("qc_plots/sim1/tsne-batch.png")
scater::plotTSNE(sim1, colour_by = "Group")
ggplot2::ggsave("qc_plots/sim1/tsne-group.png")
scater::plotTSNE(sim1, colour_by = "total")
ggplot2::ggsave("qc_plots/sim1/tsne-total.png")
scater::plotTSNE(sim1_qc, colour_by = "Batch")
ggplot2::ggsave("qc_plots/sim1/tsne-batch-qc.png")
scater::plotTSNE(sim1_qc, colour_by = "Group")
ggplot2::ggsave("qc_plots/sim1/tsne-group-qc.png")
scater::plotTSNE(sim1_qc, colour_by = "total")
ggplot2::ggsave("qc_plots/sim1/tsne-total-qc.png")
scater::plotTSNE(sim1, colour_by = "Discard")
ggplot2::ggsave("qc_plots/sim1/tsne-filt.png")
scater::plotColData(sim1, "sum", "detected", colour_by = "Discard",
                    other_fields = "Batch") +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~ Batch, scales = "free")
ggplot2::ggsave("qc_plots/sim1/filtering.png")

message("Converting to AnnData...")
anndata <- reticulate::import("anndata")
sim1_adata <- anndata$AnnData(
    X = t(counts(sim1)),
    obs = data.frame(colData(sim1)),
    var = data.frame(rowData(sim1))
)
sim1_qc_adata <- anndata$AnnData(
    X = t(counts(sim1_qc)),
    obs = data.frame(colData(sim1_qc)),
    var = data.frame(rowData(sim1_qc))
)

message("Saving simulation...")
fs::dir_create("simulations")
# Remove intermediate matrices to reduce file size
assays(sim1) <- assays(sim1)[assayNames(sim1) == "counts"]
assays(sim1_qc) <- assays(sim1_qc)[assayNames(sim1_qc) == "counts"]
saveRDS(sim1, "simulations/sim1.Rds")
saveRDS(sim1_qc, "simulations/sim1_qc.Rds")
sim1_adata$write(filename = "simulations/sim1.h5ad")
sim1_qc_adata$write(filename = "simulations/sim1_qc.h5ad")
