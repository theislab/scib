import pathlib

#configfile: "config.yaml"

""" Specifications
+ General idea:
    + multiple data scenarios that need to go through the same process
        1. data scenario
        + method to output format mapping
        + scaled/unscaled
        + HVG/full feature
        + output type
    + metrics on all integrated outputs and settings
    + cell cycle variance and other analysis on integrated output
+ input files need to follow a specific structure
    + all files should be under a ROOT directory following the above folder structure
    + use symlinks to avoid copying the data, but if data is rerun, please save according to given folder + + structure: 
    + 1 folder per data scenario, must contain:
        + folder for integrated output
        + folder for metrics output
        + folder for cell cycle variance and other analysis output
"""

## variables and function definition for input files
ROOT = config["ROOT"]
DATA_SCENARIOS = config["DATA_SCENARIOS"]
SCALING = ['scaled'] #, 'unscaled']
FEATURE_SELECTION = {'hvg': 2000, 'full_feature': 0}

def get_from_scenario(scenario, key="file"):
    return DATA_SCENARIOS[scenario][key]

def as_list(x):
    return x if isinstance(x, list) else [x]

def join_path(*args):
    path = pathlib.Path(args[0])
    for d in args[1:]:
        path = path / d
    return str(path)

def get_filename_pattern(file_type, level):
    """
    file_type: one of ['integration', 'metrics', 'cc_variance']
    level: one of ['single', 'final', 'by_method']
    """
    
    if file_type == 'integration':
        suffix = "{method}.h5ad"
    elif file_type in ['metrics', 'cc_variance']:
        suffix = "{method}_{o_type}.csv"
    else:
        raise ValueError(f"{file_type} is not a valid file type")
    
    if level == 'single':
        return join_path(ROOT,  "{scenario}", file_type, "{scaling}", "{hvg}", suffix)
    elif level == 'by_method':
        return join_path(ROOT, "{{scenario}}", file_type, "{{scaling}}", "{{hvg}}", 
                         "{method}_{o_type}.csv")
    elif level == 'final':
        return join_path(ROOT, f"{file_type}.csv")
    else:
        raise ValueError(f"{level} is not a valid level")

def get_all_file_patterns(file_type, exclude_type=''):
    """
    expand file patterns
    file_type: one of ['integrated', 'metrics', 'cc_variance']
    exclude_type: output type to be excluded
    """
    all_files = []
    for m, ot in config["METHODS"].items():
    
        # remove types to be excluded
        ot = [x for x in as_list(ot) if x not in as_list(exclude_type)]
        if not ot:
            continue
        
        file_pattern = get_filename_pattern(file_type, "by_method")
        if file_type == 'integration':
            expanded = expand(file_pattern, method=m)
        elif file_type in ['metrics', 'cc_variance']:
            expanded = expand(file_pattern, method=m, o_type=ot)
        else:
            raise ValueError(f"{file_type} is not a valid file type")
        
        for p in expanded:
            f = expand(p, scenario=DATA_SCENARIOS, scaling=SCALING, hvg=FEATURE_SELECTION.keys())
            all_files.extend(f)
    
    return all_files

## RULES
rule all:
 input:
     get_filename_pattern("metrics", "final"),
     get_filename_pattern("cc_variance", "final")

## INTEGRATION
# TODO: run integration using integration script
# rule integration:
#     output: touch(join_path(ROOT, "integration.done"))

## METRICS
rule metrics:
    input:
        tables = get_all_file_patterns("metrics"),
        script = "scripts/merge_metrics.py"
    output: get_filename_pattern("metrics", "final")
    shell: "python {input.script} -i {input.tables} -o {output} --root {ROOT}"

rule metrics_single:
    input:
        u      = lambda wildcards: get_from_scenario(wildcards.scenario, key="file"),
        i      = get_filename_pattern("integration", "single"),
        script = "scripts/metrics.py"
    output: get_filename_pattern("metrics", "single")
    params:
        batch_key = lambda wildcards: get_from_scenario(wildcards.scenario, key="batch_key"),
        label_key = lambda wildcards: get_from_scenario(wildcards.scenario, key="label_key"),
        organism  = lambda wildcards: get_from_scenario(wildcards.scenario, key="organism"),
        assay     = lambda wildcards: get_from_scenario(wildcards.scenario, key="assay"),
        hvgs      = lambda wildcards: FEATURE_SELECTION[wildcards.hvg]
    shell:
        """
        OUT_DIR=$(dirname {output})
        python {input.script} -u {input.u} -i {input.i} -o $OUT_DIR \
        -b {params.batch_key} -l {params.label_key} --type {wildcards.o_type} \
        --hvgs {params.hvgs} --organism {params.organism} --assay {params.assay} -v
        """

## Cell Cycle Variation
rule cc_variation:
    input:
        tables = get_all_file_patterns("cc_variance", exclude_type='knn'),
        script = "scripts/merge_cc_variance.py"
    output: get_filename_pattern("cc_variance", "final")
    shell: "python {input.script} -i {input.tables} -o {output} --root {ROOT}"

rule cc_single:
    input:
        u      = lambda wildcards: get_from_scenario(wildcards.scenario, key="file"),
        i      = get_filename_pattern("integration", "single"),
        script = "scripts/cell_cycle_variance.py"
    output: get_filename_pattern("cc_variance", "single")
    params:
        batch_key = lambda wildcards: get_from_scenario(wildcards.scenario, key="batch_key"),
        organism  = lambda wildcards: get_from_scenario(wildcards.scenario, key="organism"),
        assay     = lambda wildcards: get_from_scenario(wildcards.scenario, key="assay"),
        hvgs      = lambda wildcards: FEATURE_SELECTION[wildcards.hvg]
    shell:
        """
        python {input.script} -u {input.u} -i {input.i} -o {output} \
        -b {params.batch_key} --assay {params.assay} --type {wildcards.o_type} \
        --hvgs {params.hvgs} --organism {params.organism}
        """
