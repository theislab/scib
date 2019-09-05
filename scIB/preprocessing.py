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

### Quality Control        
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

### Normalisation
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


def hvg_intersect(adata, batch, target_genes=2000, flavor='cell_ranger', n_bins=20, adataOut=False, n_stop=8000, min_genes=500, step_size=1000):
### Feature Selection
    """
    params:
        adata:
        batch: adata.obs column
        target_genes: maximum number of genes (intersection reduces the number of genes)
        min_genes: minimum number of intersection HVGs targeted
        step_size: step size to increase HVG selection per dataset
    return:
        list of highly variable genes less or equal to `target_genes`
    """
    
    checkAdata(adata)
    checkBatch(batch, adata.obs)
    
    intersect = None
    enough = False
    n_hvg = target_genes
    
    split = splitBatches(adata, batch) 
    hvg_res = []

    for i in split:
        sc.pp.filter_genes(i, min_cells=1) # remove genes unexpressed (otherwise hvg might break)
        hvg_res.append(sc.pp.highly_variable_genes(i, flavor='cell_ranger', n_top_genes=n_hvg, inplace=False))

    while not enough:
        genes = []

        for i in range(len(split)):

            dispersion_norm = hvg_res[i]['dispersions_norm']
            dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
            dispersion_norm[::-1].sort()
            disp_cut_off = dispersion_norm[n_hvg-1]
            gene_subset = np.nan_to_num(hvg_res[i]['dispersions_norm']) >= disp_cut_off

            genes.append(set(split[i].var[gene_subset].index))

        intersect = genes[0].intersection(*genes[1:])
        if len(intersect)>=target_genes:
            enough=True
        else:
            if n_hvg>n_stop:
                if len(intersect) < min_genes:
                    raise Exception(f'Only {len(intersect)} HVGs were found in the intersection.\n'
                                    f'This is fewer than {min_genes} HVGs set as the minimum.\n'
                                    'Consider raising `n_stop` or reducing `min_genes`.')
                break
            n_hvg=int(n_hvg+step_size)

    if adataOut:
        return adata[:,list(intersect)].copy()

    return list(intersect)


def hvg_batch(adata, batch_key=None, target_genes=2000, flavor='cell_ranger', n_bins=20, adataOut=False):
    """

    Method to select HVGs based on mean dispersions of genes that are highly 
    variablein at least all but one batch. Using a the top target_genes per 
    batch. If target genes still hasn't been reached, then HVGs in all but two 
    batches are used to fill up. This is continued until HVGs in a single batch
    are considered.
    """
    
    checkAdata(adata)
    if batch_key:
        checkBatch(batch_key, adata.obs)
    
    adata_hvg = adata if adataOut else adata.copy()

    n_batches = len(adata_hvg.obs[batch_key].cat.categories)

    # Calculate double target genes per dataset
    sc.pp.highly_variable_genes(adata_hvg,
                                flavor=flavor, 
                                n_top_genes=target_genes,
                                n_bins=n_bins, 
                                batch_key=batch_key)

    nbatch1_dispersions = adata_hvg.var['dispersions_norm'][adata_hvg.var.highly_variable_nbatches >
                                                           len(adata_hvg.obs[batch_key].cat.categories)-1]
    
    nbatch1_dispersions.sort_values(ascending=False, inplace=True)

    if len(nbatch1_dispersions) > target_genes:
        hvg = nbatch1_dispersions.index[:target_genes]
    
    else:
        enough = False
        hvg = nbatch1_dispersions.index[:]
        not_n_batches = 1
        
        while not enough:
            print('Got in here')
            target_genes_diff = target_genes - len(hvg)

            tmp_dispersions = adata_hvg.var['dispersions_norm'][adata_hvg.var.highly_variable_nbatches ==
                                                                (len(adata_hvg.obs[batch_key].cat.categories)-not_n_batches)]

            if len(tmp_dispersions) < target_genes_diff:
                print(f'Using {len(tmp_dispersions)} HVGs from n_batch-{not_n_batches} set')
                hvg.append(tmp_dispersions.index)
                not_n_batches += 1

            else:
                print(f'Using {target_genes_diff} HVGs from n_batch-{not_n_batches} set')
                tmp_dispersions.sort_values(ascending=False, inplace=True)
                hvg = hvg.append(tmp_dispersions.index[:target_genes_diff])
                enough=True

    print(f'Using {len(hvg)} HVGs')

    if not adataOut:
        del adata_hvg
        return hvg.tolist()
    else:
        return adata_hvg[:,hvg].copy()


### Feature Reduction
def reduce_data(adata, subset=False,
                hvg=True, filter=True, batch_key=None, flavor='cell_ranger', n_top_genes=2000, n_bins=20,
                pca=True, pca_comps=50,
                neighbors=True, use_rep='X_pca', 
                paga=False, paga_groups='cell_type', 
                umap=True,
                tsne=False,
                diffmap=False,
                draw_graph=False):
    
    checkAdata(adata)
    if batch_key:
        checkBatch(batch_key, adata.obs)
    
    if hvg:
        print("HVG")
        # quick fix: HVG doesn't work on dense matrix
        if not sparse.issparse(adata.X):
            adata.X = sparse.csr_matrix(adata.X)
        if filter:
            sc.pp.filter_genes(adata, min_cells=1)
        sc.pp.highly_variable_genes(adata, flavor=flavor, n_top_genes=n_top_genes, n_bins=n_bins, batch_key=batch_key)
        n_hvg = np.sum(adata.var["highly_variable"])
        print(f'Computed {n_hvg} highly variable genes')
        
    if pca:
        print("PCA")
        sc.tl.pca(adata,
                  n_comps=pca_comps, 
                  use_highly_variable=hvg, 
                  svd_solver='arpack', 
                  return_info=True)
    
    if neighbors:
        print("Nearest Neigbours")
        sc.pp.neighbors(adata, use_rep=use_rep)

    if tsne:
        print("tSNE")
        sc.tl.tsne(adata, n_jobs=12) # n_jobs works for MulticoreTSNE, but not regular implementation
    
    if umap:
        print("UMAP")
        sc.tl.umap(adata)
    
    if paga:
        print(f'PAGA by group "{paga_groups}"')
        checkBatch(paga_groups, adata.obs)
        sc.tl.paga(adata, groups=paga_groups)
    
    if diffmap:
        sc.tl.diffmap(adata)
    
    if draw_graph:
        sc.tl.draw_graph(adata)
        
### Cell Cycle
def cc_tirosh(marker_gene_file, adata=None):
    """
    Tirosh et al. cell cycle marker genes downloaded from
    https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt
    return: (s_genes, g2m_genes)
        s_genes: S-phase genes
        g2m_genes: G2- and M-phase genes
    """
    cell_cycle_genes = [x.strip().lower().capitalize() for x in open(marker_gene_file)]
    if adata:
        cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    # split list into S-phase and G2/M-phase genes
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    
    return s_genes, g2m_genes
