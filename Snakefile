from scripts.snakemake_parse import *
import pathlib

#configfile: "config.yaml"
cfg = ParsedConfig(config)

wildcard_constraints:
    hvg = "hvg|full_feature"

rule all:
    input:
        cfg.get_filename_pattern("metrics", "final")

rule integration:
    input:
        cfg.get_all_file_patterns("integration")
    message: "Integration done"

rule integration_prepare:
    input:
        adata  = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="file"),
        script = "scripts/runPP.py"
    output:
        join_path(cfg.get_filename_pattern("prepare", "directory_by_setting"), "adata_pre.{prep}")
    message:
        """
        Preparing adata
        wildcards: {wildcards}
        parameters: {params}
        output: {output}
        """
    params:
        batch_key = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="batch_key"),
        hvgs      = lambda wildcards: cfg.get_feature_selection(wildcards.hvg),
        scale     = lambda wildcards: "-s" if wildcards.scaling == "scaled" else "",
        rout      = lambda wildcards: "-r" if wildcards.prep == "RDS" else "",
        seurat    = lambda wildcards: "-l" if wildcards.prep == "RDS" else "",
        cmd       = f"conda run -n {cfg.py_env} python"
    benchmark:
        join_path(cfg.get_filename_pattern("prepare", "directory_by_setting"),
                  "prep_{prep}.benchmark")
    shell:
        """
        {params.cmd} {input.script} -i {input.adata} -o {output} -b {params.batch_key} \
        --hvgs {params.hvgs} {params.scale} {params.rout} {params.seurat}
        """

def get_prep_adata(wildcards):
    """
    get R or python adata file depending on integration method
    """
    if cfg.get_from_method(wildcards.method, "R"):
        prep = "RDS"
    else:
        prep = "h5ad"
    return expand(rules.integration_prepare.output, **wildcards, prep=prep)

# ------------------------------------------------------------------------------
# Python specific integration rule.
# TODO: decorate with some detailed information
# ------------------------------------------------------------------------------

def get_celltype_option_for_integration(wildcards):
    if cfg.get_from_method(wildcards.method, "use_celltype"):
        label_key = cfg.get_from_scenario(wildcards.scenario, key="label_key")
        return f"-c {label_key}"
    return ""

rule integration_run_python:
    input:
        adata  = get_prep_adata,
        pyscript = "scripts/runIntegration.py"
    output:
        cfg.get_filename_pattern("integration", "single", "h5ad")
    message:
        """
        Run {wildcards.method} on {wildcards.scaling} data
        feature selection: {wildcards.hvg}
        dataset: {wildcards.scenario}
        command: {params.cmd}
        hvgs: {params.hvgs}
        cell type option: {params.cell_type}
        """
    params:
        batch_key = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="batch_key"),
        cell_type = get_celltype_option_for_integration,
        hvgs      = lambda wildcards, input: cfg.get_hvg(wildcards, input.adata[0]),
        timing    = "-t" if cfg.timing else "",
        cmd       = f"conda run -n {cfg.py_env} python"
    benchmark:
        f'{cfg.get_filename_pattern("integration", "single", "h5ad")}.benchmark'
    shell:
        """
        {params.cmd} {input.pyscript} -i {input.adata} -o {output} \
	      -b {params.batch_key} --method {wildcards.method} {params.hvgs} {params.cell_type} \
	      {params.timing}
        """

# ------------------------------------------------------------------------------
# R specific integration rule.
# TODO: decorate with some detailed information
# ------------------------------------------------------------------------------
rule integration_run_r:
    input:
        adata  = get_prep_adata,
        rscript = "scripts/runMethods.R"
    output:
        cfg.get_filename_pattern("integration", "single", "rds")
    message:
        """
        Run {wildcards.method} on {wildcards.scaling} data
        feature selection: {wildcards.hvg}
        dataset: {wildcards.scenario}
        command: {params.cmd}
        hvgs: {params.hvgs}
        """
    params:
        batch_key = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="batch_key"),
        hvgs      = lambda wildcards, input: cfg.get_hvg(wildcards, input.adata[0]),
        cmd       = f"conda run -n {cfg.r_env} Rscript",
        timing    = "-t" if cfg.timing else ""
    benchmark:
        f'{cfg.get_filename_pattern("integration", "single", "rds")}.benchmark'
    shell:
        """
        {params.cmd} {input.rscript} -i {input.adata} -o {output} -b {params.batch_key} \
            --method {wildcards.method} {params.hvgs} {params.timing}
        """

# ------------------------------------------------------------------------------
# Simply converts the RDS files created by the R scripts to h5ad files for
# further processing with the metrics rule
# ------------------------------------------------------------------------------
rule convert_RDS_h5ad:
    input:
        i = cfg.get_filename_pattern("integration", "single", "rds"),
        script = "scripts/runPost.py"
    output:
        cfg.get_filename_pattern("integration", "single", "rds_to_h5ad")
    message:
        """
        Convert integrated data from {wildcards.method} into h5ad
        """
    params:
        cmd = f"conda run -n {cfg.conv_env} python"
    shell:
        """
        if [ {wildcards.method} == "conos" ]
        then
            {params.cmd} {input.script} -i {input.i} -o {output} -c
        else
            {params.cmd} {input.script} -i {input.i} -o {output}
        fi
        """

# ------------------------------------------------------------------------------
# Compute metrics
# ------------------------------------------------------------------------------

rule metrics_unintegrated:
    input:        cfg.get_all_file_patterns("metrics_unintegrated")
    message:        "Collect all unintegrated metrics"

rule metrics_integrated:
    input: cfg.get_all_file_patterns("metrics")
    message: "Collect all integrated metrics"

all_metrics = rules.metrics_integrated.input
if cfg.unintegrated_m:
    all_metrics.extend(rules.metrics_unintegrated.input)

rule metrics:
    input:
        tables = all_metrics,
        script = "scripts/merge_metrics.py"
    output:
        cfg.get_filename_pattern("metrics", "final")
    message: "Merge all metrics"
    params:
        cmd = f"conda run -n {cfg.py_env} python"
    shell: "{params.cmd} {input.script} -i {input.tables} -o {output} --root {cfg.ROOT}"

def get_integrated_for_metrics(wildcards):
    if wildcards.method == "unintegrated":
        pattern = str(rules.integration_prepare.output)
        file = os.path.splitext(pattern)[0]
        return f"{file}.h5ad"
    elif cfg.get_from_method(wildcards.method, "R"):
        return cfg.get_filename_pattern("integration", "single", "rds_to_h5ad")
    else:
        return cfg.get_filename_pattern("integration", "single", "h5ad")

rule metrics_single:
    input:
        u      = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="file"),
        i      = get_integrated_for_metrics,
        script = "scripts/metrics.py"
    output: cfg.get_filename_pattern("metrics", "single")
    message:
        """
        Metrics {wildcards}
        output: {output}
        """
    params:
        batch_key = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="batch_key"),
        label_key = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="label_key"),
        organism  = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="organism"),
        assay     = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="assay"),
        hvgs      = lambda wildcards: cfg.get_feature_selection(wildcards.hvg),
        cmd       = f"conda run -n {cfg.py_env} python"
    shell:
        """
        {params.cmd} {input.script} -u {input.u} -i {input.i} -o {output} -m {wildcards.method} \
        -b {params.batch_key} -l {params.label_key} --type {wildcards.o_type} \
        --hvgs {params.hvgs} --organism {params.organism} --assay {params.assay} -v
        """

# ------------------------------------------------------------------------------
# Cell cycle score sanity check
# ------------------------------------------------------------------------------

rule cc_variation:
    input:
        tables = cfg.get_all_file_patterns("cc_variance", output_types=["full", "embed"]),
        script = "scripts/merge_cc_variance.py"
    output: cfg.get_filename_pattern("cc_variance", "final")
    params:
        cmd = f"conda run -n {cfg.py_env} python"
    shell: "{params.cmd} {input.script} -i {input.tables} -o {output} --root {cfg.ROOT}"

rule cc_single:
    input:
        u      = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="file"),
        i      = cfg.get_filename_pattern("integration", "single"),
        script = "scripts/cell_cycle_variance.py"
    output: cfg.get_filename_pattern("cc_variance", "single")
    params:
        batch_key = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="batch_key"),
        organism  = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="organism"),
        assay     = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="assay"),
        hvgs      = lambda wildcards: cfg.get_feature_selection(wildcards.hvg),
        cmd       = f"conda run -n {cfg.py_env} python"
    shell:
        """
        {params.cmd} {input.script} -u {input.u} -i {input.i} -o {output} \
        -b {params.batch_key} --assay {params.assay} --type {wildcards.o_type} \
        --hvgs {params.hvgs} --organism {params.organism}
        """

# ------------------------------------------------------------------------------
# Merge benchmark files
#
# Run this after the main pipeline using:
# snakemake --configfile config.yaml --cores 1 benchmarks
# ------------------------------------------------------------------------------

rule benchmarks:
    input:
        script = "scripts/merge_benchmarks.py"
    output:
        cfg.get_filename_pattern("benchmarks", "final")
    message: "Merge all benchmarks"
    params:
        cmd = f"conda run -n {cfg.py_env} python"
    shell: "{params.cmd} {input.script} -o {output} --root {cfg.ROOT}"
