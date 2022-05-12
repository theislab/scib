import numpy as np
import scanpy as sc

import scib


def test_scale():
    adata = sc.datasets.blobs()
    adata.obs["blobs"] = adata.obs["blobs"].astype("category")

    scib.pp.scale_batch(adata, "blobs")
    split = scib.utils.split_batches(adata, "blobs")
    for i in split:
        assert np.allclose(i.X.mean(0), np.zeros((0, adata.n_vars)))
