import numpy as np

from tests.common import *


def test_kbet(adata_pca):
    score = scIB.me.kBET(
        adata_pca,
        batch_key='batch',
        label_key='celltype',
        embed='X_pca'
    )
    score = 1 - np.nanmean(score['kBET'])
    LOGGER.info(f"score: {score}")
    assert np.isnan(score)
