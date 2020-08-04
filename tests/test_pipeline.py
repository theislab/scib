import scIB
from scIB.tests import utils
import pathlib
import os
import subprocess
import json
import snakemake
import scanpy as sc

METHODS = {
    "mnn": {"output_type": "full"},
    "scanorama": {"output_type": ["embed", "full"]},
    "seurat": { "output_type": ["embed", "full"], "R": True},
    "trvae": {"output_type": ["embed", "full"], "no_scale": True},
    "scgen": {"output_type": ["embed", "full"], "use_celltype": True}
    # TODO: include rest of methods
}

def create_if_missing(dir):
    if not os.path.isdir(dir):
        os.mkdir(dir)

def setup_test_directory(methods):
    """
    create necessary files for a test directory
    TODO: use fixtures
    TODO: create environments
    :params methods: list of method names to be used
    """
    methods = [methods] if isinstance(methods, str) else methods
    data_dir = os.path.abspath(f"./pipeline-{'_'.join(methods)}")
    create_if_missing(data_dir)
    print(f"created {data_dir}")

    # create input and output directories
    input_dir = os.path.join(data_dir, "input")
    create_if_missing(input_dir)
    output_dir = os.path.join(data_dir, "output")
    create_if_missing(output_dir)

    # write data files
    input_adata_file = os.path.join(input_dir, "adata_raw.h5ad")
    if not os.path.isfile(input_adata_file):
        adata = utils.create_adata_dummy(pca=True, n_top_genes=2000, neighbors=True)
        adata.write(input_adata_file)

    # write config file
    config = {
        "ROOT" : output_dir,
        "r_env": "scIB_env",
        "py_env": "scIB_env",
        "conv_env": "scIB_env",
        "timing": False,
        "FEATURE_SELECTION": {
            "hvg": 2000,
            "full_feature": 0
        },
        "SCALING": ["unscaled", "scaled"],
        "METHODS": {k: METHODS[k] for k in methods},
        "DATA_SCENARIOS": {
            "test_data": {
                "batch_key": "batch",
                "label_key": "celltype",
                "organism": "human",
                "assay": "expression",
                "file": input_adata_file
            }
        }
    }
    config_file = os.path.join(data_dir, "config.json")
    with open(config_file, 'w') as f:
        json.dump(config, f)

    workdir = pathlib.Path(scIB.__file__).parent.parent

    return {
        "workdir": workdir,
        "config": config,
        "configfile": config_file,
        "data_dir": data_dir,
        "input_dir": input_dir,
        "output_dir": output_dir
    }


def run_snakemake(snakeroot, configfile, dryrun=False):
    """
    Calls snakemake pipeline given snakefile and config file
    :param snakeroot: directory containing Snakefile
    :param configfile: path to config file
    """

    cmd = ["snakemake", "--configfile", str(configfile)]
    if dryrun:
        cmd.append("-n")
    print(f"Running {cmd}")
    response = subprocess.run(cmd, cwd=snakeroot, stderr=subprocess.STDOUT)
    response.check_returncode()

    # if config is None:
    #     if configfile is None:
    #         raise ValueError("Either config or configfile must be given!")
    #     with open(str(configfile)) as f:
    #         config = json.load(f)
    # success = snakemake.snakemake(
    #     snakefile=str(snakefile),
    #     cores=1,
    #     config=config,
    #     verbose=True,
    #     **params
    # )
    # assert success

def test_snakemake_dryrun():
    paths = setup_test_directory(["scgen", "scanorama"])
    run_snakemake(
        snakeroot=paths["workdir"],
        configfile=paths["configfile"],
        dryrun=True
    )

def test_snakemake():
    paths = setup_test_directory(["scgen", "scanorama"])
    run_snakemake(
        snakeroot=paths["workdir"],
        configfile=paths["configfile"],
        dryrun=False
    )

if __name__=='__main__':
    test_snakemake_dryrun()
    test_snakemake()