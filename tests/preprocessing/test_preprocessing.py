import scanpy as sc
import numpy as np
import scib


def test_scale():
    adata = sc.datasets.blobs()
    scib.pp.scale_batch(adata, 'blobs')
    split = scib.utils.splitBatches(adata, 'blobs')
    for i in split:
        assert np.allclose(i.X.mean(0), np.zeros((0,adata.n_vars)))
