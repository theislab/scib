import numpy as np
import scanpy as sc
from tqdm import tqdm

import scib
from tests.common import LOGGER, assert_near_exact


def test_kbet(adata_pca):
    score = scib.me.kBET(
        adata_pca, batch_key="batch", label_key="celltype", embed="X_pca", type_="full"
    )
    LOGGER.info(f"score: {score}")
    assert_near_exact(score, 0.55, diff=1e-01)


def test_kbet_random(adata_pca):
    scores = []
    for _ in tqdm(range(5)):
        adata_pca.obs["batch"] = adata_pca.obs["batch"].sample(frac=1).values
        sc.pp.pca(adata_pca, n_comps=20, use_highly_variable=True)
        score = scib.me.kBET(
            adata_pca,
            batch_key="batch",
            label_key="celltype",
            embed="X_pca",
            type_="full",
        )
        LOGGER.info(f"score: {score}")
        scores.append(score)
    assert_near_exact(np.mean(scores), 0.8, diff=1e-01)
