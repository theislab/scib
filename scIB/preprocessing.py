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
            count_filter_threshold=(0,np.inf)):
    
    if count_filter_threshold == (0, np.inf):
        count_filter_threshold = count_threshold
    if gene_filter_threshold == (0, np.inf):
        gene_filter_threshold = gene_threshold
    
    # 2D scatter plot
    plot_scatter(adata, color=color, title=color,
                 gene_threshold=gene_filter_threshold[0], 
                 count_threshold=count_filter_threshold[0],
                 legend_loc=legend_loc)
    
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
def normalize(adata, min_mean = 0.1, log=True, precluster=True, sparsify=True):
    
    checkAdata(adata)

    # Check for 0 count cells
    if np.any(adata.X.sum(axis=1) == 0):
        raise ValueError('found 0 count cells in the AnnData object.'
                         ' Please filter these from your dataset.')

    # Check for 0 count genes
    if np.any(adata.X.sum(axis=0) == 0):
        raise ValueError('found 0 count genes in the AnnData object.'
                         ' Please filter these from your dataset.')

    if sparsify:
        # massive speedup when working with sparse matrix
        if not sparse.issparse(adata.X): # quick fix: HVG doesn't work on dense matrix
            adata.X = sparse.csr_matrix(adata.X)

    anndata2ri.activate()
    ro.r('library("scran")')
    
    # keep raw counts
    adata.layers["counts"] = adata.X.copy()

    is_sparse=False
    X = adata.X.T
    # convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
    if sparse.issparse(X):
        is_sparse = True
        
        if X.nnz > 2**31-1:
            X = X.tocoo()
        else:
            X = X.tocsc()

    ro.globalenv['data_mat'] = X

    if precluster:
        # Preliminary clustering for differentiated normalisation
        adata_pp = adata.copy()
        sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
        sc.pp.log1p(adata_pp)
        sc.pp.pca(adata_pp, n_comps=15, svd_solver='arpack')
        sc.pp.neighbors(adata_pp)
        sc.tl.louvain(adata_pp, key_added='groups', resolution=0.5)

        ro.globalenv['input_groups'] = adata_pp.obs['groups']
        size_factors = ro.r('sizeFactors(computeSumFactors(SingleCellExperiment('
                            'list(counts=data_mat)), clusters = input_groups,'
                            f' min.mean = {min_mean}))')

        del adata_pp

    else:
        size_factors = ro.r('sizeFactors(computeSumFactors(SingleCellExperiment('
                            f'list(counts=data_mat)), min.mean = {min_mean}))')
        
    # modify adata
    adata.obs['size_factors'] = size_factors
    adata.X /= adata.obs['size_factors'].values[:,None]
    if log:
        print("Note! Performing log1p-transformation after normalization.")
        sc.pp.log1p(adata)
    else:
        print("No log-transformation performed after normalization.")

    if is_sparse:
        # convert to sparse, bc operation always converts to dense
        adata.X = sparse.csr_matrix(adata.X)

    adata.raw = adata # Store the full data set in 'raw' as log-normalised data for statistical testing

    # Free memory in R
    ro.r('rm(list=ls())')
    ro.r('lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)')
    ro.r('invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)), '
         'detach, character.only=TRUE, unload=TRUE))')
    ro.r('gc()')

    anndata2ri.deactivate()


def scale_batch(adata, batch):
    """
    Function to scale the gene expression values of each batch separately.
    """

    checkAdata(adata)
    checkBatch(batch, adata.obs)

    # Store layers for after merge (avoids vstack error in merge)
    adata_copy = adata.copy()
    tmp = dict()
    for lay in list(adata_copy.layers):
        tmp[lay] = adata_copy.layers[lay]
        del adata_copy.layers[lay]

    split = splitBatches(adata_copy, batch)

    for i in split:
        sc.pp.scale(i)

    adata_scaled = merge_adata(split)

    # Reorder to original obs_name ordering
    adata_scaled = adata_scaled[adata.obs_names]

    # Add layers again
    for key in tmp:
        adata_scaled.layers[key] = tmp[key]

    del tmp
    del adata_copy
    
    return adata_scaled

    
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
    variable genes in all batches. Using a the top target_genes per batch by
    average normalize dispersion. If target genes still hasn't been reached, 
    then HVGs in all but one batches are used to fill up. This is continued 
    until HVGs in a single batch are considered.
    """
    
    checkAdata(adata)
    if batch_key is not None:
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
        print(f'Using {len(nbatch1_dispersions)} HVGs from full intersect set')
        hvg = nbatch1_dispersions.index[:]
        not_n_batches = 1
        
        while not enough:
            target_genes_diff = target_genes - len(hvg)

            tmp_dispersions = adata_hvg.var['dispersions_norm'][adata_hvg.var.highly_variable_nbatches ==
                                                                (n_batches-not_n_batches)]

            if len(tmp_dispersions) < target_genes_diff:
                print(f'Using {len(tmp_dispersions)} HVGs from n_batch-{not_n_batches} set')
                hvg = hvg.append(tmp_dispersions.index)
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
def reduce_data(adata, batch_key=None, subset=False,
                filter=True, flavor='cell_ranger', n_top_genes=2000, n_bins=20,
                pca=True, pca_comps=50, overwrite_hvg=True,
                neighbors=True, use_rep='X_pca', 
                umap=True):
    """
    overwrite_hvg:
        if True, ignores any pre-existing 'highly_variable' column in adata.var
        and recomputes it if `n_top_genes` is specified else calls PCA on full features.
        if False, skips HVG computation even if `n_top_genes` is specified and uses
        pre-existing HVG column for PCA
    """
    
    checkAdata(adata)
    if batch_key:
        checkBatch(batch_key, adata.obs)
    
    if n_top_genes is not None and overwrite_hvg:
        print("HVG")
        
        overwrite_hvg = False
        
        ## quick fix: HVG doesn't work on dense matrix
        if not sparse.issparse(adata.X):
            adata.X = sparse.csr_matrix(adata.X)
            
        if batch_key is not None:
            hvg_list = hvg_batch(adata, batch_key=batch_key, target_genes=n_top_genes, n_bins=n_bins)
            adata.var['highly_variable'] = np.in1d(adata.var_names, hvg_list)

        else:
            print(f"Calculating {n_top_genes} HVGs for reduce_data.")
            sc.pp.highly_variable_genes(adata,
                                        n_top_genes=n_top_genes,
                                        n_bins=n_bins,
                                        flavor=flavor)

        n_hvg = np.sum(adata.var["highly_variable"])
        print(f'Computed {n_hvg} highly variable genes')
        
    if pca:
        print("PCA")
        use_hvgs = not overwrite_hvg and "highly_variable" in adata.var
        sc.tl.pca(adata,
                  n_comps=pca_comps, 
                  use_highly_variable=use_hvgs, 
                  svd_solver='arpack', 
                  return_info=True)
    
    if neighbors:
        print("Nearest Neigbours")
        sc.pp.neighbors(adata, use_rep=use_rep)
    
    if umap:
        print("UMAP")
        sc.tl.umap(adata)
    
        
### Cell Cycle
def score_cell_cycle(adata, organism='mouse'):
    """
    Tirosh et al. cell cycle marker genes downloaded from
    https://raw.githubusercontent.com/theislab/scanpy_usage/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt
    return: (s_genes, g2m_genes)
        s_genes: S-phase genes
        g2m_genes: G2- and M-phase genes
    """
    import pathlib
    root = pathlib.Path(__file__).parent
    
    cc_files = {'mouse': [root / 'resources/s_genes_tirosh.txt',
                          root / 'resources/g2m_genes_tirosh.txt'],
                'human': [root / 'resources/s_genes_tirosh_hm.txt',
                          root / 'resources/g2m_genes_tirosh_hm.txt']}

    with open(cc_files[organism][0], "r") as f:
        s_genes = [x.strip() for x in f.readlines() if x.strip() in adata.var.index]
    with open(cc_files[organism][1], "r") as f:
        g2m_genes = [x.strip() for x in f.readlines() if x.strip() in adata.var.index]

    if (len(s_genes) == 0) or (len(g2m_genes) == 0):
        rand_choice = np.random.randint(1,adata.n_vars,10)
        rand_genes = adata.var_names[rand_choice].tolist()
        raise ValueError(f"cell cycle genes not in adata\n organism: {organism}\n varnames: {rand_genes}")
    
    sc.tl.score_genes_cell_cycle(adata, s_genes, g2m_genes)

    
def saveSeurat(adata, path, batch, hvgs=None):
    import re
    ro.r('library(Seurat)')
    ro.r('library(scater)')
    anndata2ri.activate()

    if sparse.issparse(adata.X):
        if not adata.X.has_sorted_indices:
            adata.X.sort_indices()

    for key in adata.layers:
        if sparse.issparse(adata.layers[key]):
            if not adata.layers[key].has_sorted_indices:
                adata.layers[key].sort_indices()

    ro.globalenv['adata'] = adata
    
    ro.r('sobj = as.Seurat(adata, counts="counts", data = "X")')

    # Fix error if levels are 0 and 1
    # ro.r(f'sobj$batch <- as.character(sobj${batch})')
    ro.r(f'Idents(sobj) = "{batch}"')
    ro.r(f'saveRDS(sobj, file="{path}")') 
    if hvgs is not None:
        hvg_out = re.sub('\.RDS$', '', path)+'_hvg.RDS'
        #hvg_out = path+'_hvg.rds'
        ro.globalenv['hvgs']=hvgs
        ro.r('unlist(hvgs)')
        ro.r(f'saveRDS(hvgs, file="{hvg_out}")')


    anndata2ri.deactivate()
    
       
def readSeurat(path):
    anndata2ri.activate()
    ro.r('library(Seurat)')
    ro.r('library(scater)')
    ro.r(f'sobj <- readRDS("{path}")')
    adata = ro.r('as.SingleCellExperiment(sobj)')
    anndata2ri.deactivate()

    #Test for 'X_EMB'
    if 'X_EMB' in adata.obsm:
        if 'X_emb' in adata.obsm:
            print('overwriting existing `adata.obsm["X_emb"] in the adata object')
        adata.obsm['X_emb'] = adata.obsm['X_EMB']
        del adata.obsm['X_EMB']
    
    return(adata)
    
def readConos(inPath):
    from time import time
    from shutil import rmtree
    from scipy.io import mmread
    from os import mkdir, path
    import pandas as pd
    
    dir_path = "/localscratch/conos"+str(int(time()))
    while path.isdir(dir_path):
        dir_path += '2'
    dir_path += '/'
    mkdir(dir_path)
    
    ro.r('library(conos)')
    ro.r(f'con <- readRDS("{inPath}")')
    ro.r('meta <- function(sobj) {return(sobj@meta.data)}')
    ro.r('metalist <- lapply(con$samples, meta)')
    ro.r('library(data.table)')
    ro.r('metaM <- do.call(rbind,unname(metalist))')
    ro.r(f'saveConosForScanPy(con, output.path="{dir_path}", pseudo.pca=TRUE, pca=TRUE, metadata.df=metaM)')
    gene_df = pd.read_csv(dir_path + "genes.csv")

    metadata = pd.read_csv(dir_path + "metadata.csv")
    metadata.index = metadata.CellId
    del metadata["CellId"]

    embedding_df = pd.read_csv(dir_path + "embedding.csv")
    # Decide between using PCA or pseudo-PCA
    pseudopca_df = pd.read_csv(dir_path + "pseudopca.csv")
    #pca_df = pd.read_csv(dir_path + "pca.csv")

    graph_conn_mtx = mmread(dir_path + "graph_connectivities.mtx")
    graph_dist_mtx = mmread(dir_path + "graph_distances.mtx")
    
    adata = sc.read_mtx(dir_path+ "raw_count_matrix.mtx")
    
    
    adata.var_names = gene_df["gene"].values
    adata.obs_names = metadata.index.values

    adata.obs = metadata.copy()

    # Depends on which PCA you loaded
    adata.X_pca = pseudopca_df.values
    adata.obsm['X_pca'] = pseudopca_df.values

    # Name according to embedding you saved
    adata.X_umap = embedding_df.values
    adata.obsm['X_umap'] = embedding_df.values

    adata.uns['neighbors'] = dict(connectivities=graph_conn_mtx.tocsr(), distances=graph_dist_mtx.tocsr())

    # Assign raw counts to .raw slot, load in normalised counts
    #adata.raw = adata
    #adata_temp = sc.read_mtx(DATA_PATH + "count_matrix.mtx")
    #adata.X = adata_temp.X

    rmtree(dir_path)
    
    return adata

    
