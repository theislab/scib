import pytest
import scanpy as sc

import scib


def test_mouse(adata_paul15):

    assert "S_score" not in adata_paul15.obs.columns
    assert "G2M_score" not in adata_paul15.obs.columns
    assert "phase" not in adata_paul15.obs.columns

    scib.pp.score_cell_cycle(
        adata_paul15,
        organism="mouse",
    )
    assert "S_score" in adata_paul15.obs.columns
    assert "G2M_score" in adata_paul15.obs.columns
    assert "phase" in adata_paul15.obs.columns

    scib.pp.score_cell_cycle(
        adata_paul15,
        organism="mus musculus",
    )


def test_human(adata_paul15):
    scib.pp.score_cell_cycle(
        sc.datasets.pbmc68k_reduced(),
        organism="human",
    )
    with pytest.raises(ValueError):
        scib.pp.score_cell_cycle(
            adata_paul15,
            organism="human",
        )


@pytest.mark.parametrize("adata_from_url", ["c_elegans", "zebrafish"], indirect=True)
def test_organism(adata_from_url):
    scib.pp.score_cell_cycle(
        adata_from_url,
        organism=adata_from_url.uns["dataset_name"],
    )
