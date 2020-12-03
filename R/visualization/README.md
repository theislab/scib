# Plotting Summary Tables - how to

In order to generate different types of summary tables three simple steps are necessary. Please download and place in your working directory:
1. the folder _scib/R/visualization/img_, containing icons used by _knit_table.R_;
2. the basic script _scib/R/visualization/knit_table.R_, which will be used by all plotting functions;
3. the folder _scib/R/visualization/data_, which contains examplary metrics files, result files from the paper as well as usability/scalability metrics.

Depending on the data type of interest (RNA/ATAC) different functions will be called to generate the plots. It would be then necessary to download each specific script and place it in the same _working_dir_.

## RNA
### Single-Atlas Summary Tables
Here we are interested in plotting one summary table for each RNA scenario. This summary table shows all integration methods ranked by **Overall Score**, which is calculated as weighted sum of **Batch Correction** and **Bio Conservation**. **Batch Correction** and **Bio Conservation** are average scores calculated over all respective metrics, which are shown in the columns of the summary table. Please refer to **Supplementary Figure 4** of the paper for an example.

To generate this figure you will need to:
1. have a .csv file computed by scib pipeline over one or multiple RNA scenarios, containing metrics scores. Exemplary .csv can be found in _data_ folder:
   * _metrics_RNA_allAtlases.csv_ contains results shown in the paper, for all five RNA scenarios and two simulations.
   * _metrics_RNA_immune_human.csv_ contains results specific to only one scenario, Immune (human). 
2. download and place in your _working_dir_ the script _plotSingleAtlasRNA.R_;
3. run `source('plotSingleAtlasRNA.R')`;
4. call **plotSingleAtlasRNA**, with the following parameters:
   * _csv_metrics_path_: path to a .csv file output of scib pipeline that contains the metrics calculated across one or multiple RNA scenarios. 
   * _outdir_: output directory where the plots and a .csv for each scenario containing the ranked summary table scores will be saved.
   * _weight_batch_: number in [0,1] to use as weight for the batch correction metrics. Weight for bio conservation is calculated as 1-weight_batch.

Exemplary call to the function: `plotSingleAtlasRNA(csv_metrics_path = "./data/metrics_RNA_allAtlases.csv", outdir = ".", weight_batch = 0.4)`

### Best Methods Summary Table
