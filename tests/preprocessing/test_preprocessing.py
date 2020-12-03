import scanpy as sc
import numpy as np
import scIB


def test_scale():
    adata = sc.datasets.blobs()
    scIB.pp.scale_batch(adata, 'blobs')
    split = scIB.utils.splitBatches(adata, 'blobs')
    for i in split:
        assert np.allclose(i.X.mean(0), np.zeros((0,adata.n_vars)))
