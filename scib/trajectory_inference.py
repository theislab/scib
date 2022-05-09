import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

from . import utils


def paga(adata, groups="louvain"):
    """ """
    utils.check_adata(adata)

    sc.pp.neighbors(adata)
    sc.tl.paga(adata, groups=groups)
    _ = sc.pl.paga_compare(adata, show=False)

    fig1, ax1 = plt.subplots()
    sc.pl.umap(adata, size=40, ax=ax1, show=False)
    sc.pl.paga(
        adata,
        pos=adata.uns["paga"]["pos"],
        show=False,
        node_size_scale=10,
        node_size_power=1,
        ax=ax1,
        text_kwds={"alpha": 0},
    )
    plt.show()


def dpt(adata, group, root, opt="min", comp=0):
    utils.check_adata()

    # TODO compute diffmap before

    # get root
    stem_mask = np.isin(adata.obs[group], root)
    if opt == "min":
        opt_stem_id = np.argmin(adata.obsm["X_diffmap"][stem_mask, comp])
    elif opt == "max":
        opt_stem_id = np.argmax(adata.obsm["X_diffmap"][stem_mask, comp])
    else:
        raise ("invalid optimum", opt)
    root_id = np.arange(len(stem_mask))[stem_mask][opt_stem_id]
    adata.uns["iroot"] = root_id
    # compute pseudotime
    sc.tl.dpt(adata)
