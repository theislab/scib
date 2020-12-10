# Plotting Summary Tables - how to

In order to generate different types of summary tables three simple steps are necessary. Please download and place in your working directory:
1. the basic script _scib/R/visualization/knit_table.R_, which will be used by all plotting functions;
2. the folder _scib/R/visualization/img_, containing icons used by _knit_table.R_;
3. the folder _scib/R/visualization/data_, which contains examplary metrics files, result files from the paper as well as usability/scalability metrics.

Depending on the data type of interest (RNA/ATAC) different functions will be called to generate the plots. It would be then necessary to download each specific script and place it in the same _working_dir_.

## RNA
### Single-Task Summary Tables
Here we are interested in plotting one summary table for each RNA task. This summary table shows all integration methods ranked by **Overall Score**, which is calculated as weighted sum of **Batch Correction** and **Bio Conservation**. **Batch Correction** and **Bio Conservation** are average scores calculated over all respective metrics, which are shown in the columns of the summary table. Please refer to **Supplementary Figure 4** of the paper for an example.

To generate this figure you will need to:
1. have a .csv file computed by scib over one (or multiple) RNA task, containing metrics scores. Exemplary .csv can be found in _data_ folder:
   * _metrics_RNA_allTasks.csv_ contains results shown in the paper, for all five RNA tasks and two simulations.
   * _metrics_RNA_immune_human.csv_ contains results specific to only one task: Immune (human). 
2. download and place in your _working_dir_ the script _plotSingleTaskRNA.R_;
3. run `source('plotSingleTaskRNA.R')`;
4. call **plotSingleTaskRNA**, with the following parameters:
   * _csv_metrics_path_: path to a .csv file output of scib that contains the metrics calculated across one or multiple RNA tasks. 
   * _outdir_: output directory where the plots and a .csv for each task containing the ranked summary table scores will be saved.
   * _weight_batch_: number in [0,1] to use as weight for the batch correction metrics. Weight for bio conservation is calculated as 1-weight_batch.

Exemplary call to the function: `plotSingleTaskRNA(csv_metrics_path = "./data/metrics_RNA_allTasks.csv", outdir = ".", weight_batch = 0.4)`

### Best Methods Summary Table
A second visualization output shows a summary table where only the best-performing combinations of pre-processing choices for each integration method are kept. The methods are then ranked based on their overall performances across RNA tasks. Scores related to simulation tasks, usability and scalability are also shown. Please refer to **Figure 3b** of the paper for an example.

To generate this figure you will need to:
1. have a .csv file computed by scib over multiple RNA tasks, containing metrics scores. Exemplary .csv can be found in _data_ folder:
   * _metrics_RNA_allTasks.csv_ contains results shown in the paper, for all five RNA tasks and two simulations.
2. download and place in your _working_dir_ the script _plotBestMethodsRNA.R_;
3. run `source('plotBestMethodsRNA.R')`;
4. call **plotBestMethodsRNA**, with the following parameters:
   * _csv_metrics_path_: path to a .csv file output of scib that contains the metrics calculated across multiple RNA tasks. 
   * _outdir_: output directory where the summary table (in three formats: .pdf/.tiff/.png) will be saved.
   * _csv_usability_path_: path to a .csv file containing the results of the usability analysis. Default to "/data/usability4bestMethods.csv". These scores will NOT be used for ranking best methods.
   * _csv_scalability_time_path_: path to a .csv file containing the results of the scalability analysis, regarding run time. Default to "/data/scalability_score_time.csv". These scores will NOT be used for ranking best methods.
   * _csv_scalability_memory_path_: path to a .csv file containing the results of the scalability analysis, regarding memory consumption. Default to "/data/scalability_score_memory.csv". These scores will NOT be used for ranking best methods.
   * _ids_RNA_: character vector of ids for RNA tasks, as they are named in _csv_metrics_path_.
   * _ids_simulation_: character vector of ids for simulated tasks, as they are named in _csv_metrics_path_.
   * _labels_RNA_: character vector of label names for RNA tasks, to rename ids. These names will be plotted in the summary table.
   * _labels_simulation_: character vector of label names for simulated tasks, to rename ids. These names will be plotted in the summary table.
   * _weight_batch_: number in [0,1] to use as weight for the batch correction metrics. Weight for bio conservation is calculated as 1-weight_batch.

Exemplary call to the function: `plotBestMethodsRNA(csv_metrics_path = "./data/metrics_RNA_allTasks.csv", outdir = ".", csv_usability_path = "./data/usability4bestMethods.csv", csv_scalability_time_path = "./data/scalability_score_time.csv", csv_scalability_memory_path = "./data/scalability_score_memory.csv", ids_RNA = c("pancreas", "lung_atlas", "immune_cell_hum", "immune_cell_hum_mou", "mouse_brain"), ids_simulation = c("simulations_1_1", "simulations_2"), labels_RNA = c("Pancreas", "Lung", "Immune (hum)", "Immune (hum & mou)", "Brain (mou)"), labels_simulation = c("Sim 1", "Sim 2"), weight_batch = 0.4)` 
