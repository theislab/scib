import pytest
import scanpy as sc

import scib


def test_mouse(adata_paul15):
    scib.pp.score_cell_cycle(
        adata_paul15,
        organism="mouse",
    )
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
