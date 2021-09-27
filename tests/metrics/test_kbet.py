import numpy as np

from tests.common import *


def test_kbet(adata_pca):
    score = scIB.me.kBET(
        adata_pca,
        batch_key='batch',
        label_key='celltype',
        embed='X_pca'
    )
    LOGGER.info(f"score: {score}")
    assert np.isnan(score)
