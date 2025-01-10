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


def test_zebrafish(adata_paul15):
    with pytest.raises(AssertionError):
        scib.pp.score_cell_cycle(
            adata_paul15,
            organism="zebrafish",
        )
    # TODO
