# Plotting Summary Tables - how to

In order to generate different types of summary tables three simple steps are necessary. Please download and place in your working directory:

1. The basic script `scib/R/visualization/knit_table.R`, which will be used by most plotting functions
2. The folder `scib/R/visualization/img`, containing icons used by `knit_table.R`
3. The folder `scib/R/visualization/data`, which contains example metrics files, result files from the paper as well as usability/scalability metrics.

Depending on the data type of interest (RNA/ATAC) different functions will be called to generate the plots.
It would be then necessary to download each specific script and place it in the same working directory.

## RNA

### Single-Task Summary Tables

Here we are interested in plotting one summary table for each RNA task.
This summary table shows all integration methods ranked by **Overall Score**, which is calculated as a weighted sum of **Batch Correction** and **Bio Conservation**.
**Batch Correction** and **Bio Conservation** are average scores calculated over all respective metrics, which are shown in the columns of the summary table.
Please refer to **Supplementary Figure 4** of the paper for an example.

To generate this figure you will need to:

1. Have a `.csv` file computed by scib over one (or multiple) RNA task(s), containing metrics scores. Example CSV files can be found in the `data/` folder:
   * `metrics_RNA_allTasks.csv` contains results shown in the paper, for all five RNA tasks and two simulations
   * `metrics_RNA_immune_human.csv` contains results specific to only one task: "Immune (human)"
2. Download and place the script `plotSingleTaskRNA.R` in your working directory
3. Run `source('plotSingleTaskRNA.R')`
4. Call `plotSingleTaskRNA()` with the following parameters:
   * `csv_metrics_path`: Path to a CSV file output of scib that contains the metrics calculated across one or multiple RNA tasks. 
   * `outdir`: Output directory where the plots and a CSV file for each task containing the ranked summary table scores will be saved.
   * `weight_batch`: Number between 0 and 1 to use as weight for the batch correction metrics.
     Weight for bio conservation is calculated as `1 - weight_batch`.
     Default value used in scib manuscript is 0.4.

**Example function call:**

```r
plotSingleTaskRNA(
    csv_metrics_path = "./data/metrics_RNA_allTasks.csv",
    outdir = ".",
    weight_batch = 0.4
)
```

### Best Methods Summary Table

A second visualization output shows a summary table where only the best-performing combinations of pre-processing choices for each integration method are reported.
The methods are then ranked based on their overall performances across RNA tasks.
Scores related to simulation tasks, usability and scalability are also shown.
Please refer to **Figure 3b** of the paper for an example.

To generate this figure you will need to:

1. Have a CSV file computed by scib over multiple RNA tasks, containing metrics scores.
   An example CSV file can be found in the `data/` folder:
   * `metrics_RNA_allTasks.csv` contains results shown in the paper, for all five RNA tasks and two simulations.
2. Download and place the script `plotBestMethodsRNA.R` in your working directory
3. Run `source('plotBestMethodsRNA.R')`
4. Call `plotBestMethodsRNA()` with the following parameters:
   * `csv_metrics_path`: Path to a CSV file output of scib that contains the metrics calculated across multiple RNA tasks.
   * `outdir`: Output directory where the summary tables (in three formats: PDF/TIFF/PNG) will be saved
   * `csv_usability_path`: Path to a CSV file containing the results of the usability analysis.
     Defaults to `"./data/usability4bestMethods.csv"`.
     These scores will NOT be used for ranking best methods.
   * `csv_scalability_time_path`: Path to a CSV file containing the results of the scalability analysis, regarding run time.
     DefaultS to `"/data/scalability_score_time.csv"`.
     These scores will NOT be used for ranking best methods.
   * `csv_scalability_memory_path`: Path to a CSV file containing the results of the scalability analysis, regarding memory consumption.
     Defaults to `"/data/scalability_score_memory.csv"`.
     These scores will NOT be used for ranking best methods.
   * `ids_RNA`: Character vector of ids for RNA tasks, as they are named in `csv_metrics_path`.
   * `ids_simulation`: Character vector of ids for simulated tasks, as they are named in `csv_metrics_path`.
   * `labels_RNA`: Character vector of label names for RNA tasks, to rename ids.
     These names will be plotted in the summary table.
   * `labels_simulation`: Character vector of label names for simulated tasks, to rename ids. 
     These names will be plotted in the summary table.
   * `weight_batch`: Number between 0 and 1 to use as weight for the batch correction metrics.
     Weight for bio conservation is calculated as `1 - weight_batch`.
     Default value used in scib manuscript is 0.4.

**Example function call:**

```r
plotBestMethodsRNA(
    csv_metrics_path = "./data/metrics_RNA_allTasks.csv",
    outdir = ".",
    csv_usability_path = "./data/usability4bestMethods.csv",
    csv_scalability_time_path = "./data/scalability_score_time.csv", 
    csv_scalability_memory_path = "./data/scalability_score_memory.csv",
    ids_RNA = c("pancreas", "lung_atlas", "immune_cell_hum",
                "immune_cell_hum_mou", "mouse_brain"),
    ids_simulation = c("simulations_1_1", "simulations_2"),
    labels_RNA = c("Pancreas", "Lung", "Immune (hum)", "Immune (hum & mou)",
                   "Brain (mou)"),
    labels_simulation = c("Sim 1", "Sim 2"), weight_batch = 0.4
)
``` 

## ATAC

### Single-Task, single-feature space Summary Tables

Here we are interested in plotting one summary table for each ATAC task, considering different feature spaces separately. This summary table shows all integration methods ranked by **Overall Score**, which is calculated as weighted sum of **Batch Correction** and **Bio Conservation**. **Batch Correction** and **Bio Conservation** are average scores calculated over all respective metrics, which are shown in the columns of the summary table. Please refer to **Supplementary Figure 25** of the paper for an example.

To generate this figure you will need to:

1. Have a CSV file computed by scib over one (or multiple) ATAC task, containing metrics scores.
   Example CSV files can be found in the `data/ATAC_metrics/` folder:
   * `metrics_ATAC_large11.csv` contains results of the Mouse Brain Large scenario (11 batches) over all feature spaces (genes/windows/peaks).
   * `metrics_atac_large_11batches_gene.csv` contains results of the Mouse Brain Large scenario (11 batches), specific to only one feature space (gene). 
2. Download and place the script `plotSingleTaskATAC.R` in your working directory
3. Run `source('plotSingleTaskRNA.R')`
4. Call `plotSingleTaskATAC()`, with the following parameters:
   * `csv_metrics_path`: Path to a CSV file output of scib that contains the metrics calculated across one or multiple ATAC tasks
   * `outdir`: Output directory where the plots and a CSV file for each task containing the ranked summary table scores will be saved
   * `weight_batch`: number between 0 and 1 to use as weight for the batch correction metrics.
   Weight for bio conservation is calculated as `1 - weight_batch`.
   Default value used in the scib manuscript is 0.4.

**Example function call:**

```r
plotSingleTaskATAC(
    csv_metrics_path = "./data/metrics_ATAC_large11.csv",
    outdir = ".",
    weight_batch = 0.4
)
```

### Single-Task, all features Summary Tables

Another way to plot a summary table for ATAC tasks is by considering all feature spaces. 
This summary table shows, as before, all integration methods ranked by **Overall Score**, but adds one column ("Feature Space") for genes/windows/peaks.
Please refer to **Supplementary Figure 23** of the paper for an example.

To generate this figure you will need to:

1. Have a CSV file computed by scib over one ATAC task, containing metrics scores over multiple feature spaces.
   Example csv files can be found in the `data/ATAC_metrics/` folder:
   * `metrics_ATAC_large11.csv` contains results of the Mouse Brain Large scenario (11 batches) over all feature spaces (genes/windows/peaks)
2. Download and place the script `plotSingleATAC_withFeat.R` in your working directory
3. Run `source('plotSingleATAC_withFeat.R')`
4. Call `plotSingleATAC_withFeat()`, with the following parameters:
   * `csv_metrics_path`: Path to a CSV file output of scib that contains the metrics calculated across one ATAC task, over multiple feature spaces
   * `outdir`: Output directory where the plots and a csv file containing the ranked summary table scores will be saved
   * `weight_batch`: Number between 0 and 1 to use as weight for the batch correction metrics.
   Weight for bio conservation is calculated as `1 - weight_batch`.
   Default value used in the scib manuscript is 0.4.

**Example function call:**

```r
plotSingleATAC_withFeat(
    csv_metrics_path = "./data/metrics_ATAC_large11.csv",
    outdir = ".",
    weight_batch = 0.4
)
```

### Best Methods Summary Table

Also for ATAC, it is possible to plot a summary table where only the best-performing combinations of pre-processing choices (here only influenced by Output) for each integration method are reported.
The methods are then ranked based on their overall performances across ATAC tasks.
Please refer to **Figure 4b** of the paper for an example.

To generate this figure you will need to:

1. Have a csv file computed by scib over multiple ATAC tasks, containing metrics scores. 
   Example csv files can be found in the `data/` folder:
   * `metrics_ATAC_small3_large11.csv` contains results shown in the paper, for all six ATAC tasks.
2. Download and place the script `plotBestMethodsATAC.R` in your working directory
3. Run `source('plotBestMethodsATAC.R')`
4. Call `plotBestMethodsATAC()`, with the following parameters:
   * `csv_metrics_path`: Path to a CSV file output of scib that contains the metrics calculated across multiple ATAC tasks. 
   * `outdir`: Output directory where the summary table (in three formats: PDF/TIFF/PNG) will be saved
   * `ids_ATAC`: Character vector of ids for ATAC tasks, as they are named in `csv_metrics_path`.
   * `labels_ATAC`: Character vector of label names for ATAC tasks, to rename ids.
     These names will be plotted in the summary table.
   * `weight_batch`: Number between 0 and 1 to use as weight for the batch correction metrics.
     Weight for bio conservation is calculated as `1 - weight_batch`.
     Default value used in the scib manuscript is 0.4.

**Example function call:**

```r
plotBestMethodsATAC(
    csv_metrics_path = "./data/ATAC_metrics/metrics_ATAC_small3_large11.csv",
    ids_ATAC = c("mouse_brain_atac_windows_small",
                 "mouse_brain_atac_windows_large",
                 "mouse_brain_atac_peaks_small", 
                 "mouse_brain_atac_peaks_large",
                 "mouse_brain_atac_genes_small",
                 "mouse_brain_atac_genes_large"),
    labels_ATAC = c("Brain (mou) Windows small", "Brain (mou) Windows large",
                    "Brain (mou) Peaks small", "Brain (mou) Peaks large",
                    "Brain (mou) Genes small", "Brain (mou) Genes large")
)
```

## Scatter plots

An alternative, higher-level way to summarize integration performance is by plotting the overall batch correction score against the overall bio-conservation score.
These scatter can be made either individually for each dataset or averaged over datasets.

### Summary scatter

The `makeSummaryScatter()` function produces a faceted scatter plot where each facet is a dataset.
It calls the `loadScores()` function to load the metrics data followed by the `plotSummaryScatter()` function to produce the plot.
The figure is then saved to disk in various formats.

**Parameters**

* `scores_files`: A vector of paths to `_summary_scores.csv` files produced by `plotSingleTaskRNA()` or `plotSingleTaskATAC()`
* `dataset_key`: Named character vector giving names for datasets
* `out_dir`: Path to output directory

**Example usage:**

_Assuming `summaries/` is a directory containing `_summary_score.csv` files_

```r
source("plotSummaryScatter.R")
source("exampleKeys.R")

scores_files <- fs::dir_ls("summaries")
dataset_key <- getDatasetKey()
methods_pal <- getMethodsPal()

makeSummaryScatter(
    scores_files = scores_files,
    dataset_key = dataset_key,
    methods_pal = methods_pal,
    out_dir = "."
)
```

### Best methods scatter

The `makeBestScatter()` function produces a plot showing the performance of a selection of methods across several datasets.
It calls the `loadBestScores()` function to load the metrics data followed by the `plotBestScatter()` function to produce the plot.
The figure is then saved to disk in various formats.

**Parameters**

* `scores_files`: A vector of paths to `_summary_scores.csv` files produced by `plotSingleTaskRNA()` or `plotSingleTaskATAC()`
* `best_methods`: Character vector giving the names of methods to plot
* `type`: Either `"RNA"` or `"ATAC"`.
  Used to select datasets when reading `scores_files` if both are present.
* `out_dir`: Path to output directory

**Example usage:**

_Assuming `summaries/` is a directory containing `_summary_score.csv` files_

```r
source("plotBestScatter.R")
source("exampleKeys.R")

scores_files <- fs::dir_ls("summaries")
best_methods <- getBestMethods()

makeBestScatter(
    scores_files = scores_files,
    best_methods = best_methods,
    type = "RNA",
    out_dir = "."
)
```
