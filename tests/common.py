import logging
import os
import subprocess
import warnings

import scanpy as sc

warnings.filterwarnings("ignore")

CORES = "1"
LOGGER = logging.getLogger(__name__)


def assert_near_exact(x, y, diff=1e-5):
    assert abs(x - y) <= diff, f"{x} != {y} with error margin {diff}"


def create_if_missing(dir):
    if not os.path.isdir(dir):
        os.mkdir(dir)


def run(
    cmd,
    dir_path,
    report_stdout=True,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    **kwargs,
):
    LOGGER.info(cmd)
    shell = isinstance(cmd, str)
    response = subprocess.run(
        cmd,
        cwd=dir_path,
        shell=shell,
        stdout=stdout,
        stderr=stderr,
        universal_newlines=True,
        **kwargs,
    )
    try:
        response.check_returncode()
    except subprocess.CalledProcessError:
        if report_stdout:
            LOGGER.error(response.stdout)
        LOGGER.error(response.stderr)
        raise
    return response


def runR(r_cmd, dir_path, report_stdout=False):
    return run(f"Rscript -e '{r_cmd}'", dir_path, report_stdout=report_stdout)


def add_embed(adata, type_):
    if type_ == "pca":
        if "X_pca" in adata.obsm:
            mtx = adata.obsm["X_pca"]
        else:
            mtx = sc.tl.pca(adata, copy=True).obsm["X_pca"]
    elif type_ == "full":
        mtx = adata.X
    else:
        raise ValueError(f"'{type_}' not a valid embedding type")
    adata.obsm["X_emb"] = mtx
    return adata
