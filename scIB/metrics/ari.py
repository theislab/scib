import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score

from scIB.utils import checkAdata, checkBatch


def ari(adata, group1, group2):
    """
    params:
        adata: anndata object
        group1: ground-truth cluster assignments (e.g. cell type labels)
        group2: "predicted" cluster assignments
    The function is symmetric, so group1 and group2 can be switched
    """

    checkAdata(adata)

    if isinstance(group1, str):
        checkBatch(group1, adata.obs)
        group1 = adata.obs[group1].tolist()
    elif isinstance(group1, pd.Series):
        group1 = group1.tolist()

    if isinstance(group2, str):
        checkBatch(group2, adata.obs)
        group2 = adata.obs[group2].tolist()
    elif isinstance(group2, pd.Series):
        group2 = group2.tolist()

    if len(group1) != len(group2):
        raise ValueError(f'different lengths in group1 ({len(group1)}) and group2 ({len(group2)})')

    return adjusted_rand_score(group1, group2)
