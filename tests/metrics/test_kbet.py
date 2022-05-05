import numpy as np

import scib
from tests.common import *


def test_kbet(adata_pca):
    score = scib.me.kBET(
        adata_pca, batch_key="batch", label_key="celltype", embed="X_pca"
    )
    LOGGER.info(f"score: {score}")
    assert np.isnan(score)
