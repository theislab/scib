import numpy as np
import scanpy as sc

import scib
from tests.common import LOGGER


def test_scale():
    adata = sc.datasets.blobs()
    adata.obs["blobs"] = adata.obs["blobs"].astype("category")

    scib.pp.scale_batch(adata, "blobs")
    split = scib.utils.split_batches(adata, "blobs")
    for i in split:
        assert np.allclose(i.X.mean(0), np.zeros((0, adata.n_vars)))


def test_merge_adatas(adata_paul15_template):
    adata = adata_paul15_template

    # create ambiguous column names
    adata.var[["1", "2"]] = (0, 2)
    adata.var = adata.var.rename(
        columns={name: "ambig_var" for name in adata.var.columns.values}
    )
    adata.obs = adata.obs.rename(
        columns={name: "ambig_obs" for name in adata.obs.columns.values}
    )

    LOGGER.debug(adata)

    scib.utils.merge_adata(adata, adata, index_unique=None)
