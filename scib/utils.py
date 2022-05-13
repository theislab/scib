import anndata


# checker functions for data sanity
def check_adata(adata):
    if type(adata) is not anndata.AnnData:
        raise TypeError("Input is not a valid AnnData object")


def check_batch(batch, obs, verbose=False):
    if batch not in obs:
        raise ValueError(f"column {batch} is not in obs")
    elif verbose:
        print(f"Object contains {obs[batch].nunique()} batches.")


def check_hvg(hvg, adata_var):
    if type(hvg) is not list:
        raise TypeError("HVG list is not a list")
    else:
        if not all(i in adata_var.index for i in hvg):
            raise ValueError("Not all HVGs are in the adata object")


def check_sanity(adata, batch, hvg):
    check_adata(adata)
    check_batch(batch, adata.obs)
    if hvg is not None:
        check_hvg(hvg, adata.var)


def split_batches(adata, batch, hvg=None, return_categories=False):
    """Split batches and preserve category information

    :param adata:
    :param batch: name of column in ``adata.obs``. The data type of the column must be of ``Category``.
    :param hvg: list of highly variable genes
    :param return_categories: whether to return the categories object of ``batch``
    """
    split = []
    batch_categories = adata.obs[batch].cat.categories
    if hvg is not None:
        adata = adata[:, hvg]
    for i in batch_categories:
        split.append(adata[adata.obs[batch] == i].copy())
    if return_categories:
        return split, batch_categories
    return split


def merge_adata(*adata_list, **kwargs):
    """Merge adatas from list while remove duplicated ``obs`` and ``var`` columns

    :param adata_list: ``anndata`` objects to be concatenated
    :param kwargs: arguments to be passed to ``anndata.AnnData.concatenate``
    """

    if len(adata_list) == 1:
        return adata_list[0]

    # Make sure that adatas do not contain duplicate columns
    for _adata in adata_list:
        for attr in ("obs", "var"):
            df = getattr(_adata, attr)
            dup_mask = df.columns.duplicated()
            if dup_mask.any():
                print(
                    f"Deleting duplicated keys `{list(df.columns[dup_mask].unique())}` from `adata.{attr}`."
                )
                setattr(_adata, attr, df.loc[:, ~dup_mask])

    return anndata.AnnData.concatenate(*adata_list, **kwargs)


def todense(adata):
    import scipy

    if isinstance(adata.X, scipy.sparse.csr_matrix):
        adata.X = adata.X.todense()
