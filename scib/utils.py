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


def merge_adata(adata_list, sep="-"):
    """
    merge adatas from list and remove duplicated obs and var columns
    """

    if len(adata_list) == 1:
        return adata_list[0]

    adata = adata_list[0].concatenate(
        *adata_list[1:], index_unique=None, batch_key="tmp"
    )
    del adata.obs["tmp"]

    if len(adata.obs.columns) > 0:
        # if there is a column with separator
        if sum(adata.obs.columns.str.contains(sep)) > 0:
            columns_to_keep = [
                name.split(sep)[1] == "0" for name in adata.var.columns.values
            ]
            clean_var = adata.var.loc[:, columns_to_keep]
        else:
            clean_var = adata.var

    if len(adata.var.columns) > 0:
        if sum(adata.var.columns.str.contains(sep)) > 0:
            adata.var = clean_var.rename(
                columns={name: name.split("-")[0] for name in clean_var.columns.values}
            )

    return adata


def todense(adata):
    import scipy

    if isinstance(adata.X, scipy.sparse.csr_matrix):
        adata.X = adata.X.todense()
