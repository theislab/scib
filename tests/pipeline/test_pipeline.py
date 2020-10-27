from tests.common import *
import pathlib
import json

METHODS = {
    "mnn": {"output_type": "full"},
    "scanorama": {"output_type": ["embed", "full"]},
    "seurat": {"output_type": ["embed", "full"], "R": True},
    "trvae": {"output_type": ["embed", "full"], "no_scale": True},
    "scgen": {"output_type": ["full"], "use_celltype": True}
    # TODO: include rest of methods
}


@pytest.fixture(scope="module", params=[["scgen", "scanorama", "trvae"]])
def run_dir(request, tmpdir_factory, adata_factory):
    """
        create necessary files for a test directory
        :params methods: list of method names to be used
        """
    data_dir = tmpdir_factory.mktemp("pipeline")
    # create input and output directories
    input_dir = data_dir.mkdir("input")
    output_dir = data_dir.mkdir("output")

    # write data files
    input_adata_file = data_dir.join("adata_raw.h5ad")
    adata = adata_factory(pca=True, n_top_genes=2000, neighbors=True)
    adata.write(input_adata_file)

    config = {
        "ROOT": str(output_dir),
        "r_env": "scIB-R-integration",
        "py_env": "scIB-python",
        "conv_env": "scIB-python",
        "timing": False,
        "FEATURE_SELECTION": {
            "hvg": 2000,
            "full_feature": 0
        },
        "SCALING": ["unscaled", "scaled"],
        "METHODS": {k: METHODS[k] for k in request.param},
        "DATA_SCENARIOS": {
            "test_data": {
                "batch_key": "batch",
                "label_key": "celltype",
                "organism": "mouse",
                "assay": "expression",
                "file": str(input_adata_file)
            }
        }
    }

    # write config file
    config_file = data_dir.join("config.json")
    with config_file.open('w') as f:
        f.write(json.dumps(config, indent=2))

    yield {
        "workdir": pathlib.Path(scIB.__file__).parent.parent,
        "config": config,
        "configfile": config_file,
        "data_dir": data_dir,
        "input_dir": input_dir,
        "adata_raw": input_adata_file,
        "output_dir": output_dir
    }

    # cleanup
    data_dir.remove()


def run_snakemake(run_dir, args=[]):
    snakeroot = run_dir["workdir"]
    configfile = run_dir["configfile"]
    if isinstance(args, str):
        args = [args]
    cmd = ["snakemake", "--configfile", str(configfile), "-j", CORES]
    cmd.extend(args)
    return run(cmd, snakeroot)


def test_dryrun(run_dir):
    run_snakemake(run_dir, args=["-n"])


@pytest.fixture(scope="package")
def test_integration(run_dir):
    run_snakemake(run_dir, args=["integration"])
    return run_dir


@pytest.fixture(scope="package")
def test_metrics(run_dir, test_integration):
    run_snakemake(run_dir, args=["metrics"])
    return run_dir
