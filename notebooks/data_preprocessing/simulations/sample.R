#' Sample group proportions
#'
#' Downsample a Splat simulation to get the desired cell proportions
#'
#' @param sim SingleCellExperiment containing a Splat simulation
#' @param batch_ncells Final number of cells in each batch
#' @param batch_props Final cell proportions in each batch
#'
#' @return Downsampled SingleCellExperiment object
sample_group_props <- function(sim, batch_ncells, batch_props) {

    # Get the cell metadata
    col_data <- as.data.frame(colData(sim))

    # Get the current counts of each group in each batch
    sim_batch_group_count <- table(col_data$Batch, col_data$Group)

    # Get the final number of each group in each batch
    batch_group_ncells <- lapply(seq_along(batch_ncells), function(idx) {
        batch_ncells[idx] * batch_props[[idx]]
    })

    # Check that we aren't trying to get more cells than we have
    lapply(seq_along(batch_group_ncells), function(idx) {
        if (any(batch_group_ncells[[idx]] > sim_batch_group_count[idx, ])) {
            stop("Not enough cells for these proportions in batch ", idx)
        }
    })

    # Downsample cells
    selected <- lapply(seq_along(batch_group_ncells), function(batch) {
        group_cells <- batch_group_ncells[[batch]]
        is_batch <- col_data$Batch == paste0("Batch", batch)
        batch_groups <- col_data$Group[is_batch]
        # Downsample batch
        selected_batch <- sapply(seq_along(group_cells), function(group) {
            is_group <- batch_groups == paste0("Group", group)
            sample(col_data$Cell[is_batch][is_group], group_cells[group])
        })
        unlist(selected_batch)
    })
    selected <- as.character(unlist(selected))

    # Subset SingleCellExperiment
    sim[, selected]
}
