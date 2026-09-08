import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import scib


@pytest.fixture
def cycle_data():
    rng = np.random.default_rng(8)
    scores = rng.normal(size=(60, 2))
    return ad.AnnData(
        scores,
        obs=pd.DataFrame(
            {
                "batch": pd.Categorical(np.repeat(["a", "b"], 30)),
                "S_score": scores[:, 0],
                "G2M_score": scores[:, 1],
            },
            index=[f"cell_{i}" for i in range(60)],
        ),
    )


@pytest.mark.parametrize("representation", ["dense", "sparse", "embedding"])
def test_cell_cycle_aligns_observations(cycle_data, representation):
    before = cycle_data
    after = before.copy()
    embed = None
    if representation == "sparse":
        before.X = csr_matrix(before.X)
        after.X = csr_matrix(after.X)
    elif representation == "embedding":
        embed = "X_integrated"
        after.obsm[embed] = after.X.copy()
    after = after[np.random.default_rng(61).permutation(after.n_obs)].copy()
    original_names = after.obs_names.copy()
    score = scib.me.cell_cycle(
        before,
        after,
        batch_key="batch",
        embed=embed,
        recompute_cc=False,
        n_comps=2,
        linreg_method="numpy",
        n_threads=1,
    )
    assert score == pytest.approx(1.0, abs=1e-12)
    assert after.obs_names.equals(original_names)


@pytest.mark.parametrize(
    "invalid",
    ["different_cells", "different_batches", "duplicate_pre", "duplicate_post"],
)
def test_cell_cycle_rejects_ambiguous_pairing(cycle_data, invalid):
    before = cycle_data
    after = before.copy()
    if invalid == "different_cells":
        after.obs_names = [f"other_{i}" for i in range(after.n_obs)]
    elif invalid == "different_batches":
        after.obs["batch"] = after.obs["batch"].iloc[::-1].to_numpy()
    else:
        target = before if invalid == "duplicate_pre" else after
        target.obs_names = ["duplicate"] * target.n_obs
    with pytest.raises(ValueError, match="obs_names|same cells"):
        scib.me.cell_cycle(
            before,
            after,
            batch_key="batch",
            recompute_cc=False,
            n_comps=2,
            linreg_method="numpy",
            n_threads=1,
        )
