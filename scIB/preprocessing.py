import scanpy as sc
from scIB.utils import *


def normalize(adata, color_col='batch', min_mean = 0.1):
    import_rpy2()
    ro.r('library("scran")')
    
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15, svd_solver='arpack')
    sc.pp.neighbors(adata_pp)
    sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5)
    
    ro.globalenv['data_mat'] = adata.X.T
    ro.globalenv['input_groups'] = adata_pp.obs['groups']
    size_factors = ro.r(f'computeSumFactors(data_mat, clusters = input_groups, min.mean = {min_mean})')
    
    # modify adata
    adata.obs['size_factors'] = size_factors
    adata.layers["counts"] = adata.X.copy()

    adata.X /= adata.obs['size_factors'].values[:,None]
    sc.pp.log1p(adata)
    adata.raw = adata # Store the full data set in 'raw' as log-normalised data for statistical testing    
    utils.summarize_counts(adata, mito=False) # average counts for normalised counts
    
def reduce_data(adata, subset=False,
                hvg=True, flavor='cell_ranger', n_top_genes=4000, bins=20,
                pca=True,
                neighbors=True, 
                paga=False, paga_groups='batch', 
                umap=True,
                tsne=False,
                diffmap=False,
                draw_graph=False):
    
    
    
    if hvg:
        sc.pp.highly_variable_genes(adata, flavor=flavor, n_top_genes=n_top_genes, n_bins=bins, subset=subset)
        n_hvg = np.sum(adata.var["highly_variable"])
        print(f'\nNumber of highly variable genes: {n_hvg}')
    if pca:
        sc.pp.pca(adata, n_comps=50, use_highly_variable=hvg, svd_solver='arpack')
    sc.pp.neighbors(adata)
    if tsne:
        sc.tl.tsne(adata, n_jobs=12) # n_jobs works for MulticoreTSNE, but not regular implementation
    if umap:
        sc.tl.umap(adata)
    if paga:
        print(f'Compute PAGA by group "{paga_groups}"')
        sc.tl.paga(adata, groups=paga_groups)
    if diffmap:
        sc.tl.diffmap(adata)
    if draw_graph:
        sc.tl.draw_graph(adata)

