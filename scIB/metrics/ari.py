import numpy as np
import pandas as pd
import scipy.special
from sklearn.metrics.cluster import adjusted_rand_score

from scIB.utils import checkAdata, checkBatch


def ari(adata, group1, group2, implementation='sklearn'):
    """
    params:
        adata: anndata object
        group1: ground-truth cluster assignments (e.g. cell type labels)
        group2: "predicted" cluster assignments
    The function is symmetric, so group1 and group2 can be switched
    """

    checkAdata(adata)
    checkBatch(group1, adata.obs)
    checkBatch(group2, adata.obs)

    group1 = adata.obs[group1].to_numpy()
    group2 = adata.obs[group2].to_numpy()

    if len(group1) != len(group2):
        raise ValueError(
            f'different lengths in group1 ({len(group1)}) and group2 ({len(group2)})'
        )

    if implementation == 'sklearn':
        return adjusted_rand_score(group1, group2)

    def binom_sum(x, k=2):
        return scipy.special.binom(x, k).sum()

    n = len(group1)
    contingency = pd.crosstab(group1, group2)

    ai_sum = binom_sum(contingency.sum(axis=0))
    bi_sum = binom_sum(contingency.sum(axis=1))

    index = binom_sum(np.ravel(contingency))
    expected_index = ai_sum * bi_sum / binom_sum(n, 2)
    max_index = 0.5 * (ai_sum + bi_sum)

    return (index - expected_index) / (max_index - expected_index)
