import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy import sparse

# rpy2 for running R code
import rpy2.rinterface_lib.callbacks
import logging
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR) # Ignore R warning messages
import rpy2.robjects as ro
import anndata2ri

# access to other methods of this module
from scIB.utils import *


def summarize_counts(adata, count_matrix=None, mt_gene_regex='^MT-'):
    
    checkAdata(adata)
    
    if count_matrix is None:
        count_matrix = adata.X
    adata.obs['n_counts'] = count_matrix.sum(1)
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    adata.obs['n_genes'] = (count_matrix > 0).sum(1)

    if mt_gene_regex != None:
        # for each cell compute fraction of counts in mito genes vs. all genes
        mito_genes = adata.var_names.str.match(mt_gene_regex)
        mt_sum = np.sum(adata[:, mito_genes].X, axis=1)
        total_sum = np.sum(adata.X, axis=1)
        # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
        if sparse.issparse(adata.X):
            mt_sum = mt_sum.A1
            total_sum = total_sum.A1
        adata.obs['percent_mito'] =  mt_sum / total_sum 
        
        #mt_gene_mask = [gene.startswith('mt-') for gene in adata.var_names]
        #mt_count = count_matrix[:, mt_gene_mask].sum(1)
        #if mt_count.ndim > 1:
        #    mt_count = np.squeeze(np.asarray(mt_count))
        #adata.obs['mt_frac'] = mt_count/adata.obs['n_counts']

def plot_scatter(adata, count_threshold=0, gene_threshold=0,
                 color=None, title='', lab_size=15, tick_size=11, legend_loc='right margin',
                palette=None):
    
    checkAdata(adata)
    if color:
        checkBatch(color, adata.obs)
    
    ax = sc.pl.scatter(adata, 'n_counts', 'n_genes', color=color, show=False,
                       legend_fontweight=50, legend_loc=legend_loc, palette=palette)
    ax.set_title(title, fontsize=lab_size)
    ax.set_xlabel("Count depth",fontsize=lab_size)
    ax.set_ylabel("Number of genes",fontsize=lab_size)
    ax.tick_params(labelsize=tick_size)
    
    if gene_threshold > 0:
        ax.axhline(gene_threshold, 0,1, color='red')
    if count_threshold > 0:
        ax.axvline(count_threshold, 0,1, color='red')
    
    fig = plt.gcf()
    cbar_ax = fig.axes[-1]
    cbar_ax.tick_params(labelsize=tick_size)
    f1 = ax.get_figure()
    plt.show()
    
#Data quality summary plots        
def plot_count_filter(adata, obs_col='n_counts', bins=60, lower=0, upper=np.inf, filter_lower=0, filter_upper=np.inf):
    
    plot_data = adata.obs[obs_col]
    
    sns.distplot(plot_data, kde=False, bins=bins)
    
    if lower > 0:
        plt.axvline(lower, linestyle = '--', color = 'g')
    if filter_lower > 0:
        plt.axvline(filter_lower, linestyle = '-', color = 'r')
    if not np.isinf(upper):
        plt.axvline(upper, linestyle = '--', color = 'g')
    if not np.isinf(upper):
        plt.axvline(filter_upper, linestyle = '-', color = 'r')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,2))
    plt.show()
    
    # determine lower bound of total, look at points below lower bound
    if filter_lower > 0:
        print(f"lower threshold: {filter_lower}")
        sns.distplot(plot_data[plot_data < lower], kde=False, bins=bins)
        plt.axvline(filter_lower, linestyle = '-', color = 'r')
        plt.axvline(lower, linestyle = '--', color = 'g')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,2))
        plt.show()

    # determine upper bound of total
    if not np.isinf(filter_upper) and not np.isinf(upper):
        print(f"upper threshold: {filter_upper}")
        sns.distplot(plot_data[plot_data > upper], kde=False, bins=bins)
        plt.axvline(filter_upper, linestyle = '-', color = 'r')
        plt.axvline(upper, linestyle = '--', color = 'g')
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,2))
        plt.show()
        
def plot_QC(adata, color=None, bins=60, legend_loc='right margin', histogram=True,
            gene_threshold=(0,np.inf), 
            gene_filter_threshold=(0,np.inf),
            count_threshold=(0,np.inf), 
            count_filter_threshold=(0,np.inf), 
            palette=sc.pl.palettes.godsnot_64):
    
    if count_filter_threshold == (0, np.inf):
        count_filter_threshold = count_threshold
    if gene_filter_threshold == (0, np.inf):
        gene_filter_threshold = gene_threshold
    
    # 2D scatter plot
    plot_scatter(adata, color=color, title=color,
                 gene_threshold=gene_filter_threshold[0], 
                 count_threshold=count_filter_threshold[0],
                 legend_loc=legend_loc, palette=palette)
    
    if not histogram:
        return

    if count_filter_threshold != (0, np.inf):
        print(f"Counts Threshold: {count_filter_threshold}")
        # count filtering
        plot_count_filter(adata, obs_col='n_counts', bins=bins,
                          lower = count_threshold[0],
                          filter_lower = count_filter_threshold[0],
                          upper = count_threshold[1],
                          filter_upper = count_filter_threshold[1])
    
    if gene_filter_threshold != (0, np.inf):
        print(f"Gene Threshold: {gene_filter_threshold}")
        # gene filtering
        plot_count_filter(adata, obs_col='n_genes', bins=bins,
                          lower = gene_threshold[0],
                          filter_lower = gene_filter_threshold[0],
                          upper = gene_threshold[1],
                          filter_upper = gene_filter_threshold[1])
    

def normalize(adata, min_mean = 0.1):
    
    checkAdata(adata)
    
    # massive speedup when working with sparse matrix
    if not sparse.issparse(adata.X): # quick fix: HVG doesn't work on dense matrix
        adata.X = sparse.csr_matrix(adata.X)
    
    anndata2ri.activate()
    ro.r('library("scran")')
    
    # keep raw counts
    adata.layers["counts"] = adata.X.copy()
    
    # Preliminary clustering for differentiated normalisation
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15, svd_solver='arpack')
    sc.pp.neighbors(adata_pp)
    sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5)
    
    ro.globalenv['data_mat'] = adata.X.T
    ro.globalenv['input_groups'] = adata_pp.obs['groups']
    size_factors = ro.r(f'computeSumFactors(data_mat, clusters = input_groups, min.mean = {min_mean})')
    del adata_pp
    
    # modify adata
    adata.obs['size_factors'] = size_factors
    adata.X /= adata.obs['size_factors'].values[:,None]
    sc.pp.log1p(adata)
    # convert to sparse, bc operation always converts to dense
    adata.X = sparse.csr_matrix(adata.X)
    adata.raw = adata # Store the full data set in 'raw' as log-normalised data for statistical testing

def subsetHVG(adata, batch, number):
    ## does not work yet, use hvg_intersect
    import scanpy as sc
    hvg = sc.pp.highly_variable_genes(adata, n_top_genes=number, batch_key=batch, flavor='cell_ranger', inplace=False)
    return hvg

def hvg_intersect(adata, batch, num=4000):
    split = splitBatches(adata, batch)
    hvg = []
    for i in split:
        tmp = sc.pp.highly_variable_genes(i, flavor='cell_ranger', n_top_genes=num, inplace=False)
        hvg.append(set(i.var[[j[0] for j in tmp]].index))
    return list(hvg[0].intersection(*hvg[1:]))

def reduce_data(adata, subset=False,
                hvg=True, flavor='cell_ranger', n_top_genes=4000, bins=20,
                pca=True,
                neighbors=True, 
                paga=False, paga_groups='batch', 
                umap=True,
                tsne=False,
                diffmap=False,
                draw_graph=False):
    
    checkAdata(adata)
    print(f"paga_groups: {paga_groups}")
    checkBatch(paga_groups, adata.obs)
    
    if hvg:
        if not sparse.issparse(adata.X): # quick fix: HVG doesn't work on dense matrix
            adata.X = sparse.csr_matrix(adata.X)
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

def cc_tirosh(marker_gene_file):
    """
    Tirosh et al. cell cycle marker genes downloaded from
    https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt
    return: (s_genes, g2m_genes)
        s_genes: S-phase genes
        g2m_genes: G2- and M-phase genes
    """
    cell_cycle_genes = [x.strip().lower().capitalize() for x in open(marker_gene_file)]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    # split list into S-phase and G2/M-phase genes
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    
    return s_genes, g2m_genes
