from scripts.snakemake_parse import *
import pathlib

configfile: "config.yaml"
cfg = ParsedConfig(config)


rule all:
    input:
        cfg.get_filename_pattern("metrics", "final")

rule integration:
    input: 
        expand(cfg.get_filename_pattern("integration", "single"),
                scenario = cfg.get_all_scenarios(),
                scaling  = cfg.get_all_scalings(),
                hvg      = cfg.get_all_feature_selections(),
                method   = cfg.get_all_methods())
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
        seurat    = lambda wildcards: "-l" if wildcards.prep == "RDS" else ""
    shell:
        """
        python {input.script} -i {input.adata} -o {output} -b {params.batch_key} \
            --hvgs {params.hvgs} {params.scale} {params.rout} {params.seurat}
        """

def get_prep_adata(wildcards):
    """
    get R or python adata file depending on integration method
    """
    if cfg.get_from_method(wildcards.method, "R") == True:
        prep = "RDS"
    else:
        prep = "h5ad"
    return expand(rules.integration_prepare.output, **wildcards, prep=prep)

rule integration_run:
    input:
        adata  = get_prep_adata,
        pyscript = "scripts/runIntegration.py",
        rscript = "R/runMethods.R"
    output:
        cfg.get_filename_pattern("integration", "single")
    params:
        batch_key = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="batch_key"),
        hvgs      = lambda wildcards: cfg.get_feature_selection(wildcards.hvg),
        cmd       = lambda wildcards: "Rscript" if cfg.get_from_method(wildcards.method, "R")
                                       else "python",
        timing    = "-t" if cfg.timing else ""
    shell:
        """
        if [ {params.cmd} = "Rscript" ]
        then
            SCRIPT={input.rscript}
        else
            SCRIPT={input.pyscript}
        fi
        
        {params.cmd} $SCRIPT -i {input.adata} -o {output} -b {params.batch_key} \
            --method {wildcards.method} --hvgs {params.hvgs} {params.timing}
        """

rule metrics:
    input:
        tables = cfg.get_all_file_patterns("metrics"),
        script = "scripts/merge_metrics.py"
    output: cfg.get_filename_pattern("metrics", "final")
    message: "Merge all metrics"
    shell: "python {input.script} -i {input.tables} -o {output} --root {cfg.ROOT}"

rule metrics_single:
    input:
        u      = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="file"),
        i      = cfg.get_filename_pattern("integration", "single"),
        script = "scripts/metrics.py"
    output: cfg.get_filename_pattern("metrics", "single")
    message: "Metrics {wildcards}"
    params:
        batch_key = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="batch_key"),
        label_key = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="label_key"),
        organism  = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="organism"),
        assay     = lambda wildcards: cfg.get_from_scenario(wildcards.scenario, key="assay"),
        hvgs      = lambda wildcards: cfg.get_feature_selection(wildcards.hvg)
    shell:
        """
        OUT_DIR=$(dirname {output})
        python {input.script} -u {input.u} -i {input.i} -o $OUT_DIR \
        -b {params.batch_key} -l {params.label_key} --type {wildcards.o_type} \
        --hvgs {params.hvgs} --organism {params.organism} --assay {params.assay} -v
        """


rule cc_variation:
    input:
        tables = cfg.get_all_file_patterns("cc_variance", output_types=["full", "embed"]),
        script = "scripts/merge_cc_variance.py"
    output: cfg.get_filename_pattern("cc_variance", "final")
    shell: "python {input.script} -i {input.tables} -o {output} --root {cfg.ROOT}"

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
        hvgs      = lambda wildcards: cfg.get_feature_selection(wildcards.hvg)
    shell:
        """
        python {input.script} -u {input.u} -i {input.i} -o {output} \
        -b {params.batch_key} --assay {params.assay} --type {wildcards.o_type} \
        --hvgs {params.hvgs} --organism {params.organism}
        """

