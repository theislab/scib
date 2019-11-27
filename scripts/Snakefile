import os

configfile: "config.yaml"
# config file must contain the following keys:
#  DATASETS: list of <file prefix>_<output_type> in yaml format
#  UNINTEGRATED: path to unintegrated adata object
#  INTEGRATED_DIR: path to all integrated adata output of format {method}.h5ad
#  METRICS_DIR: output directory of metrics
#  batch_key: name of batch column
#  label_key: name of cell type column
#  organism: name of organism (mouse|human)

METRICS_DIR = config["METRICS_DIR"]

all_metrics = os.path.join(METRICS_DIR, "all.csv")

rule all_metrics:
    input: 
        expand(os.path.join(METRICS_DIR, "{dataset}_metrics.csv"), dataset=config["DATASETS"])
    output:
        all_metrics
    shell:
        "python all_metrics.py -i {METRICS_DIR} -o {all_metrics} "

rule metrics:
    input:
        u = config["UNINTEGRATED"],
        i = os.path.join(config["INTEGRATED_DIR"], "{method}.h5ad")
    output:
        os.path.join(METRICS_DIR, "{method}_{type}_metrics.csv")
    shell:
        "python metrics.py -u {input.u} -i {input.i} -o {METRICS_DIR} "
        "-b {config[batch_key]} -l {config[label_key]} --type {wildcards.type} "
        "--organism {config[organism]} -v "

