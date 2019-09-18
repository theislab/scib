import numpy as np
from scipy import sparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata
from scIB.utils import *
from scIB.preprocessing import hvg_intersect, score_cell_cycle
from scIB.clustering import opt_louvain

import rpy2.rinterface_lib.callbacks
import logging
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR) # Ignore R warning messages
import rpy2.robjects as ro
import anndata2ri


### Silhouette score
def silhouette(adata, group_key='cell_type', metric='euclidean', embed='X_pca', scale=True):
    """
    wrapper for sklearn silhouette function values range from [-1, 1] with 1 being an ideal fit, 0 indicating overlapping clusters and -1 indicating misclassified cells
    """
    import sklearn.metrics as scm
    
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')
    asw = scm.silhouette_score(adata.obsm[embed], adata.obs[group_key])
    if scale:
        asw = (asw + 1)/2
    return asw

def silhouette_batch(adata, batch_key, group_key, metric='euclidean', 
                     embed='X_pca', verbose=True, scale=True):
    """
    Silhouette score of batch labels subsetted for each group.
    params:
        batch_key: batches to be compared against
        group_key: group labels to be subsetted by e.g. cell type
        metric: see sklearn silhouette score
        embed: name of column in adata.obsm
    returns:
        all scores: absolute silhouette scores per group label
        group means: if `mean=True`
    """
    import sklearn.metrics as scm

    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')
    
    sil_all = pd.DataFrame(columns=['group', 'silhouette_score'])
    
    for group in adata.obs[group_key].unique():
        adata_group = adata[adata.obs[group_key] == group]
        if adata_group.obs[batch_key].nunique() == 1:
            continue
        sil_per_group = scm.silhouette_samples(adata_group.obsm[embed],
                                               adata_group.obs[batch_key],
                                               metric=metric)
        # take only absolute value
        sil_per_group = [abs(i) for i in sil_per_group]
        if scale:
            # scale s.t. highest number is optimal
            sil_per_group = [1 - i for i in sil_per_group]
        d = pd.DataFrame({'group' : [group]*len(sil_per_group), 'silhouette_score' : sil_per_group})
        sil_all = sil_all.append(d)    
    sil_all = sil_all.reset_index(drop=True)
    sil_means = sil_all.groupby('group').mean()
    
    if verbose:
        print(f'mean silhouette per cell: {sil_means}')
    return sil_all, sil_means

def plot_silhouette_score(adata_dict, batch_key, group_key, metric='euclidean', 
                     embed='X_pca', palette='Dark2', per_group=False, verbose=True):
    """
    params:
        adata_dict: dictionary of adata objects, each labeled by e.g. integration method name
    """
    
    with sns.color_palette(palette):
        for label, adata in adata_dict.items():
            checkAdata(adata)
            sil_scores = silhouette(adata, 
                                          batch_key=batch_key,
                                          group_key=group_key,
                                          metric=metric,
                                          embed=embed,
                                          means=False,
                                          verbose=verbose)
            sns.distplot(sil_scores['silhouette_score'], label=label, hist=False)
        plt.title('Silhouette scores per cell for all groups')
        plt.show()
        
        if per_group:
            for data_set, adata in adata_dict.items():
                sil_scores = silhouette(adata,
                                              batch_key=batch_key,
                                              group_key=group_key,
                                              metric=metric,
                                              embed=embed,
                                              means=False,
                                              verbose=verbose)
                # plot for all groups
                for group in sil_scores['group'].unique():
                    group_scores = sil_scores[sil_scores['group'] == group]
                    sns.distplot(group_scores['silhouette_score'], label=group, hist=False)
                plt.title(f'Silhouette scores per cell for {data_set}')
                plt.show()

### NMI normalised mutual information
def nmi(adata, group1, group2, method="arithmetic", nmi_dir=None):
    """
    Normalized mutual information NMI based on 2 different cluster assignments `group1` and `group2`
    params:
        adata: Anndata object
        group1: column name of `adata.obs` or group assignment
        group2: column name of `adata.obs` or group assignment
        method: NMI implementation
            'max': scikit method with `average_method='max'`
            'min': scikit method with `average_method='min'`
            'geometric': scikit method with `average_method='geometric'`
            'arithmetic': scikit method with `average_method='arithmetic'`
            'Lancichinetti': implementation by A. Lancichinetti 2009 et al.
            'ONMI': implementation by Aaron F. McDaid et al. (https://github.com/aaronmcdaid/Overlapping-NMI) Hurley 2011
        nmi_dir: directory of compiled C code if 'Lancichinetti' or 'ONMI' are specified as `method`. Compilation should be done as specified in the corresponding README.
    return:
        normalized mutual information (NMI)
    """
    
    checkAdata(adata)
    
    if isinstance(group1, str):
        checkBatch(group1, adata.obs)
        group1 = adata.obs[group1].tolist()
    elif isinstance(group1, pd.Series):
        group1 = group1.tolist()
        
    if isinstance(group2, str):
        checkBatch(group2, adata.obs)
        group2 = adata.obs[group2].tolist()
    elif isinstance(group2, pd.Series):
        group2 = group2.tolist()
    
    if len(group1) != len(group2):
        raise ValueError(f'different lengths in group1 ({len(group1)}) and group2 ({len(group2)})')
    
    # choose method
    if method in ['max', 'min', 'geometric', 'arithmetic']:
        from sklearn.metrics import normalized_mutual_info_score
        nmi_value = normalized_mutual_info_score(group1, group2, average_method=method)
    elif method == "Lancichinetti":
        nmi_value = nmi_Lanc(group1, group2, nmi_dir=nmi_dir)
    elif method == "ONMI":
        nmi_value = onmi(group1, group2, nmi_dir=nmi_dir)
    else:
        raise ValueError(f"Method {method} not valid")
    
    return nmi_value

def onmi(group1, group2, nmi_dir=None, verbose=True):
    """
    Based on implementation https://github.com/aaronmcdaid/Overlapping-NMI
    publication: Aaron F. McDaid, Derek Greene, Neil Hurley 2011
    params:
        nmi_dir: directory of compiled C code
    """
    
    if nmi_dir is None:
        raise FileNotFoundError("Please provide the directory of the compiled C code from https://sites.google.com/site/andrealancichinetti/mutual3.tar.gz")
    
    import subprocess
    import os
    
    group1_file = write_tmp_labels(group1, to_int=False)
    group2_file = write_tmp_labels(group2, to_int=False)
    
    nmi_call = subprocess.Popen(
        [nmi_dir+"onmi", group1_file, group2_file], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT)
    
    stdout, stderr = nmi_call.communicate()
    if stderr:
        print(stderr)
    
    nmi_out = stdout.decode()
    if verbose:
        print(nmi_out)
    
    nmi_split = [x.strip().split('\t') for x in nmi_out.split('\n')]
    nmi_max = float(nmi_split[0][1])
    
    # remove temporary files
    os.remove(group1_file)
    os.remove(group2_file)
    
    return nmi_max


def nmi_Lanc(group1, group2, nmi_dir="external/mutual3/", verbose=True):
    """
    paper by A. Lancichinetti 2009
    https://sites.google.com/site/andrealancichinetti/mutual
    recommended by Malte
    """
    
    if nmi_dir is None:
        raise FileNotFoundError("Please provide the directory of the compiled C code from https://sites.google.com/site/andrealancichinetti/mutual3.tar.gz")
    
    import subprocess
    import os
    
    group1_file = write_tmp_labels(group1, to_int=False)
    group2_file = write_tmp_labels(group2, to_int=False)
    
    nmi_call = subprocess.Popen(
        [nmi_dir+"mutual", group1_file, group2_file], 
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT)
    
    stdout, stderr = nmi_call.communicate()
    if stderr:
        print(stderr)
    nmi_out = stdout.decode().strip()
    
    return float(nmi_out.split('\t')[1])

def write_tmp_labels(group_assignments, to_int=False, delim='\n'):
    """
    write the values of a specific obs column into a temporary file in text format
    needed for external C NMI implementations (onmi and nmi_Lanc functions), because they require files as input
    params:
        to_int: rename the unique column entries by integers in range(1,len(group_assignments)+1)
    """
    import tempfile
    
    if to_int:
        label_map = {}
        i = 1
        for label in set(group_assignments):
            label_map[label] = i
            i += 1
        labels = delim.join([str(label_map[name]) for name in group_assignments])
    else:
        labels = delim.join([str(name) for name in group_assignments])
        
    clusters = {label:[] for label in set(group_assignments)}
    for i, label in enumerate(group_assignments):
        clusters[label].append(str(i))
    
    output = '\n'.join([' '.join(c) for c in clusters.values()])
    output = str.encode(output)
    
    # write to file
    with tempfile.NamedTemporaryFile(delete=False) as f:
        f.write(output)
        filename = f.name
    
    return filename

### ARI adjusted rand index
def ari(adata, group1, group2):
    """
    params:
        adata: anndata object
        group1: ground-truth cluster assignments (e.g. cell type labels)
        group2: "predicted" cluster assignments
    The function is symmetric, so group1 and group2 can be switched
    """
    
    checkAdata(adata)
    
    if isinstance(group1, str):
        checkBatch(group1, adata.obs)
        group1 = adata.obs[group1].tolist()
    elif isinstance(group1, pd.Series):
        group1 = group1.tolist()
        
    if isinstance(group2, str):
        checkBatch(group2, adata.obs)
        group2 = adata.obs[group2].tolist()
    elif isinstance(group2, pd.Series):
        group2 = group2.tolist()
    
    if len(group1) != len(group2):
        raise ValueError(f'different lengths in group1 ({len(group1)}) and group2 ({len(group2)})')
    
    from sklearn.metrics.cluster import adjusted_rand_score
    return adjusted_rand_score(group1, group2)
    
### Highly Variable Genes conservation
def hvg_overlap(adata_post, adata_pre, batch, n_hvg=500):
    hvg_post = adata_post.var.index
    
    adata_pre_list = scIB.utils.splitBatches(adata_pre, batch, hvg=hvg_post)
    adata_post_list = scIB.utils.splitBatches(adata_post, batch)
    overlap = []
    
    for i in range(len(adata_pre_list)):#range(len(adata_pre_list)):
        n_hvg_tmp = np.minimum(n_hvg, int(0.5*len(adata_pre_list[i].var)))
        sc.pp.filter_genes(adata_pre_list[i], min_cells=1) # remove genes unexpressed (otherwise hvg might break)
        sc.pp.filter_genes(adata_post_list[i], min_cells=1)
        hvg_pre = sc.pp.highly_variable_genes(adata_pre_list[i], flavor='cell_ranger', n_top_genes=n_hvg_tmp, inplace=False)
        tmp_pre = adata_pre_list[i].var.index[hvg_pre['highly_variable']]
        hvg_post = sc.pp.highly_variable_genes(adata_post_list[i], flavor='cell_ranger', n_top_genes=n_hvg_tmp, inplace=False)
        tmp_post = adata_post_list[i].var.index[hvg_post['highly_variable']]
        #print(len(set(tmp_pre).intersection(set(tmp_post))))
        overlap.append((len(set(tmp_pre).intersection(set(tmp_post))))/n_hvg_tmp)
    return np.mean(overlap)

### Cell cycle effect
def cell_cycle(adata_pre, adata_post, batch_key, hvgs=2000, flavor='cell_ranger', embed=None, agg_func=np.mean, organism='mouse', n_comps=50):
    """
    params:
        adata_pre, adata_post: adatas before and after integration
        organism: 'mouse' or 'human' for choosing cell cycle genes
        agg_func: any function that takes a list of numbers and aggregates them into a single number
    """
    checkAdata(adata_pre)
    checkAdata(adata_post)
    
    if embed == 'X_pca':
        embed = None
    
    score_cell_cycle(adata_pre, organism=organism)

    scores = []
    for batch in adata_pre.obs[batch_key].unique():
        # perform PCR per batch
        raw_sub = adata_pre[adata_pre.obs[batch_key] == batch]
        int_sub = adata_post[adata_post.obs[batch_key] == batch]
        
        # select highly variable genes from non-integrated
        sc.pp.highly_variable_genes(raw_sub, n_top_genes=hvgs, flavor=flavor)
        raw_sub = raw_sub[:, raw_sub.var['highly_variable']]
        # select HVG or take embedding from integrated
        if embed is None:
            if 'highly_variable' in int_sub.var:
                int_sub = int_sub.X[:, int_sub.var['highly_variable']]
            else:
                int_sub = int_sub.X
        else:
            int_sub = int_sub.obsm[embed]
        
        covariate = raw_sub.obs[['S_score', 'G2M_score']]
        before = pc_regression(raw_sub.X, covariate, pca_sd=None, n_comps=n_comps, verbose=False)
        after =  pc_regression(int_sub, covariate, pca_sd=None, n_comps=n_comps, verbose=False)
        
        score = 1 - abs(after - before)/before # scaled result
        scores.append(score)
    
    return agg_func(scores)

### PC Regression
def get_hvg_indices(adata, verbose=True):
    if "highly_variable" not in adata.var.columns:
        if verbose:
            print(f"No highly variable genes computed, continuing with full matrix {adata.shape}")
        return np.array(range(adata.n_vars))
    return np.where((adata.var["highly_variable"] == True))[0]
        
def pcr_comparison(adata_pre, adata_post, covariate, hvgs=2000, flavor='cell_ranger', embed=None, n_comps=None, scale=True, verbose=False):
    """
    Compare the effect before and after integration
    params:
        adata_pre: uncorrected adata
        adata_post: integrated adata
    return:
        difference of R2Var value of PCR
    """
    
    if embed == 'X_pca':
        embed = None
    
    pcr_before = pcr(adata_pre, covariate=covariate, hvgs=hvgs, flavor=flavor, n_comps=n_comps, verbose=verbose)
    if (embed is None) and ('highly_variable' in adata_post.var):
        adata_post = adata_post[:, adata_post.var['highly_variable']]
    pcr_after = pcr(adata_post, covariate=covariate, hvgs=None, embed=embed, n_comps=n_comps, verbose=verbose)

    if scale:
        return abs(pcr_after - pcr_before)/pcr_before
    else:
        return pcr_after - pcr_before

def pcr(adata, covariate, embed=None, n_comps=None, hvgs=2000, flavor='cell_ranger', verbose=False):
    """
    PCR for Adata object
    params:
        adata: Anndata object
        embed: name of embedding in adata.obsm to use. No HVG selection, PCA will be computed on the embedding
        hvgs: specifies number of highly variable genes to compute. If None, no HVGs will be computed, else HVG selection and PCA on these HVGs will be computed
        n_comps: number of PCs if PCA should be computed. None will assume that PCA has been already computed
        covariate: key for adata.obs column to regress against
    return:
        R2Var of PCR
    """
    
    checkAdata(adata)
    checkBatch(covariate, adata.obs)
    
    if (embed is None) and (hvgs is not None):
        sc.pp.highly_variable_genes(adata, n_top_genes=hvgs, flavor=flavor)
        adata = adata[:, adata.var['highly_variable']]
        # note: HVG selection will lead to recomputation of PCA
    
    if verbose:
        print(f"covariate: {covariate}")
    batch = adata.obs[covariate]
    
    if (embed is not None) and (embed in adata.obsm):
        if verbose:
            print("PCR on embedding")
        return pc_regression(adata.obsm[embed], batch, n_comps=n_comps)
    elif ('X_pca' in adata.obsm) and (n_comps is None):
        return pc_regression(adata.obsm['X_pca'], batch, pca_sd=adata.uns['pca']['variance'])
    else:
        return pc_regression(adata.X, batch, n_comps=n_comps)

def pc_regression(data, covariate, pca_sd=None, n_comps=None, svd_solver='arpack', verbose=False):
    """
    params:
        data: expression or PCA matrix. Will be assumed to be PCA values, if pca_sd is given
        covariate: series or list of batch assignemnts
        n_comps: number of PCA components for computing PCA, only when pca_sd is not given. If no pca_sd is given and n_comps=None, comute PCA and don't reduce data
        pca_sd: iterable of variances for `n_comps` components. If `pca_sd` is not `None`, it is assumed that the matrix contains PCA values, else PCA is computed
    PCA is only computed, if variance contribution is not given (pca_sd).
    """
    
    if isinstance(data, (np.ndarray, sparse.csr_matrix)):
        matrix = data
    else:
        raise TypeError(f'invalid type: {data.__class__} is not a numpy array or sparse matrix')
    
    # perform PCA if no variance contributions are given
    if pca_sd is None:
        
        if n_comps is None or n_comps > min(matrix.shape):
            n_comps = min(matrix.shape)

        if n_comps == min(matrix.shape):
            svd_solver = 'full'
    
        if verbose:
            print("PCA")
        pca = sc.tl.pca(matrix, n_comps=n_comps, use_highly_variable=False,
                        return_info=True, svd_solver=svd_solver, copy=True)
        X_pca = pca[0].copy()
        pca_sd = pca[3].copy()
        del pca
    else:
        X_pca = matrix
        n_comps = matrix.shape[1]

    ## PC Regression
    if verbose:
        print("PC regression")
    
    # one-hot encode categorical values
    covariate = pd.get_dummies(covariate)
    
    # fit linear model for n_comps PCs
    from sklearn.linear_model import LinearRegression
    r2 = []
    for i in range(n_comps):
        lm = LinearRegression()
        lm.fit(X_pca[:, [i]], covariate)
        r2.append(lm.score(X_pca[:, [i]], covariate))
    
    Var = pca_sd**2 / sum(pca_sd**2) * 100
    R2Var = sum(r2*Var)/100
    
    return R2Var

### lisi score
def lisi(adata, matrix=None, knn=None, batch_key='sample', label_key='louvain', type_ = None,
         subsample = 0.5, hvg=False, verbose=False):
    """
    Compute lisi score (after integration)
    params:
        matrix: matrix from adata to calculate on
        covariate_key: variable to compute iLISI on
        cluster_key: variable to compute cLISI on
    return:
        pd.DataFrame with median cLISI and median iLISI scores (following the harmony paper)
    """
    
    checkAdata(adata)
    checkBatch(batch_key, adata.obs)
    checkBatch(label_key, adata.obs)
    
    anndata2ri.activate()
    ro.r("library(lisi)")

    if type_ =='knn': #get knn index matrix
        if verbose:
            print("Convert nearest neighbor matrix and distances for LISI.")
        dist_mat = sparse.find(adata.uns['neighbors']['distances'])
        nn_index = np.empty(shape=(adata.uns['neighbors']['distances'].shape[0],
                                   adata.uns['neighbors']['params']['n_neighbors']-1))
        nn_dists = nn_index
        for cell_id in np.arange(np.min(dist_mat[0]), np.max(dist_mat[0])):
            get_idx = dist_mat[0] == cell_id
            nn_index[cell_id,:] = dist_mat[1][get_idx][np.argsort(dist_mat[2][get_idx])]
            nn_dists[cell_id,:] = np.sort(dist_mat[2][get_idx])

        #turn metadata into numeric values (not categorical)
        meta_tmp = adata.obs[[batch_key, label_key]]
        meta_tmp[batch_key] = meta_tmp[batch_key].cat.codes
        meta_tmp[label_key] = meta_tmp[label_key].cat.codes

        if verbose:
            print("importing knn-graph")
        ro.globalenv['nn_indx'] = nn_index.T
        ro.globalenv['nn_dst'] = nn_dists.T
        ro.globalenv['batch'] = adata.obs[batch_key].cat.codes.values
        ro.globalenv['n_batches'] = len(np.unique(adata.obs[batch_key]))
        ro.globalenv['label'] = adata.obs[label_key].cat.codes.values
        ro.globalenv['n_labels'] = len(np.unique(adata.obs[label_key]))
        ro.globalenv['perplexity'] = 30 #LISI default
        
        if verbose:
            print("LISI score estimation")
        simpson_estimate_batch = ro.r(f"simpson.estimate_batch <- compute_simpson_index(nn_indx, nn_dst, batch, n_batches, perplexity)") #batch_label_keys)")
        simpson_estimate_label = ro.r(f"simpson.estimate_label <- compute_simpson_index(nn_indx, nn_dst, label, n_labels, perplexity)") #batch_label_keys)")
        simpson_est_batch = 1/np.squeeze(ro.r("simpson.estimate_batch"))
        simpson_est_label = 1/np.squeeze(ro.r("simpson.estimate_label"))
        d = {batch_key : simpson_est_batch, label_key : simpson_est_label}
        lisi_estimate = pd.DataFrame(data=d, index=np.arange(0,len(simpson_est_label)))    
    
    else:
        if matrix is None:
            matrix = adata.X

        if hvg:
            hvg_idx = get_hvg_indices(adata)
            if verbose:
                print(f"subsetting to {len(hvg_idx)} highly variable genes")
            matrix = matrix[:, hvg_idx]
    
        if sparse.issparse(matrix): #lisi score runs only on dense matrices (knn search)
            matrix = matrix.todense()
        
        if verbose:
            print("importing expression matrix")
        ro.globalenv['data_mtrx'] = matrix

        if verbose:
            print(f"covariates: {batch_key} and {label_key}")
        metadata = adata.obs[[batch_key, label_key]]
        ro.globalenv['metadata'] = metadata
    
        batch_label_keys = ro.StrVector([batch_key, label_key])
        ro.globalenv['batch_label_keys'] = batch_label_keys

        if verbose:
            print("LISI score estimation")
        lisi_estimate = ro.r(f"lisi.estimate <- compute_lisi(data_mtrx, metadata, batch_label_keys)") #batch_label_keys)")
    
    print(lisi_estimate)
    anndata2ri.deactivate()
    return lisi_estimate



### kBET
def kBET_single(matrix, batch, type_ = None, knn=None, subsample=0.5, heuristic=True, verbose=False):
    """
    params:
        matrix: expression matrix (at the moment: a PCA matrix, so do.pca is set to FALSE
        batch: series or list of batch assignemnts
        subsample: fraction to be subsampled. No subsampling if `subsample=None`
    returns:
        kBET p-value
    """
    if type_ != 'knn':
        if isinstance(subsample, float):
            #perform subsampling for large clusters, but not for small clusters (boundary is currently set to 50)
            if np.floor(matrix.shape[0]*subsample) >= 50:
                matrix, indices = sc.pp.subsample(matrix, fraction=subsample, copy=True)
                batch = batch[indices]
    
    anndata2ri.activate()
    ro.r("library(kBET)")
    
    if verbose:
        print("importing expression matrix")
    ro.globalenv['data_mtrx'] = matrix
    ro.globalenv['batch'] = batch
    #print(matrix.shape)
    #print(len(batch))

    if type_ == 'knn':
        ro.globalenv['knn_graph'] = knn
        ro.globalenv['k0'] = np.min([knn.shape[1], matrix.shape[0]])

    if verbose:
        print("kBET estimation")
    #k0 = len(batch) if len(batch) < 50 else 'NULL'
    if type_ == 'knn':
        batch_estimate = ro.r(f"batch.estimate <- kBET(data_mtrx, batch, knn=knn_graph, k0=k0, plot=FALSE, do.pca=FALSE, heuristic=FALSE, adapt=FALSE, verbose={str(verbose).upper()})")
    else:
         batch_estimate = ro.r(f"batch.estimate <- kBET(data_mtrx, batch, plot=FALSE, do.pca=FALSE, heuristic={str(heuristic).upper()}, verbose={str(verbose).upper()})")

    anndata2ri.deactivate()
    try:
        ro.r("batch.estimate$average.pval")[0]
    except rpy2.rinterface_lib.embedded.RRuntimeError:
        return np.nan
    else:
        return ro.r("batch.estimate$average.pval")[0]

def kBET(adata, batch_key, label_key, embed='X_pca', type_ = None,
                    hvg=False, subsample=0.5, heuristic=False, verbose=False):
    """
    Compare the effect before and after integration
    params:
        matrix: matrix from adata to calculate on
    return:
        pd.DataFrame with kBET p-values per cluster for batch
    """
    
    checkAdata(adata)
    checkBatch(batch_key, adata.obs)
    checkBatch(label_key, adata.obs)
    if type_ =='knn':
        if verbose:
            print("Convert nearest neighbor matrix for kBET.")
        dist_mat = sparse.find(adata.uns['neighbors']['distances'])
        nn_index = np.empty(shape=(adata.uns['neighbors']['distances'].shape[0], 
                                   adata.uns['neighbors']['params']['n_neighbors']-1))
        for cell_id in np.arange(np.min(dist_mat[0]), np.max(dist_mat[0])):
            get_idx = dist_mat[0] == cell_id
            nn_index[cell_id,:] = dist_mat[1][get_idx][np.argsort(dist_mat[2][get_idx])]
    
    matrix = adata.obsm[embed]
    
    if verbose:
        print(f"batch: {batch_key}")
    batch = adata.obs[batch_key]
    
    kBET_scores = {'cluster': [], 'kBET': []}
    for clus in adata.obs[label_key].unique():
        idx = np.where((adata.obs[label_key] == clus))[0]
        if type_ == 'knn':
            nn_index_tmp = nn_index[idx,:] #reduce nearest neighbor matrix to the desired indices
            nn_index_tmp[np.invert(np.isin(nn_index_tmp, idx))] = np.nan #set the rest nan
            score = kBET_single(
                matrix[idx, :],
                batch[idx],
                knn = nn_index_tmp,
                subsample=subsample,
                verbose=verbose,
                heuristic=False,
                type_ = type_
                )
        else:
            score = kBET_single(
                matrix[idx, :],
                batch[idx],
                subsample=subsample,
                verbose=verbose,
                heuristic=heuristic
                )
        kBET_scores['cluster'].append(clus)
        kBET_scores['kBET'].append(score)
    
    kBET_scores = pd.DataFrame.from_dict(kBET_scores)
    kBET_scores = kBET_scores.reset_index(drop=True)
    
    return kBET_scores

### Time and Memory
def measureTM(*args, **kwargs):
    """
    params:
        *args: function to be tested for time and memory
        **kwargs: list of function paramters
    returns:
        tuple : (memory (MB), time (s), list of *args function outputs)
    """
    import cProfile
    from pstats import Stats
    import memory_profiler
    
    prof = cProfile.Profile()
    out = memory_profiler.memory_usage((prof.runcall, args, kwargs), retval=True) 
    mem = np.max(out[0])- out[0][0]
    print(f'memory usage:{round(mem,0) } MB')
    print(f'runtime: {round(Stats(prof).total_tt,0)} s')
    return mem, Stats(prof).total_tt, out[1:]


def metrics(adata, adata_int, batch_key, label_key,
            hvgs=True, cluster_nmi=None,
            nmi_=False, ari_=False, nmi_method='arithmetic', nmi_dir=None, 
            silhouette_=False,  embed='X_pca', si_metric='euclidean',
            pcr_=False, cell_cycle_=False, organism='mouse', verbose=False,
            kBET_=False, kBET_sub=0.5, lisi_=False, type_ = None
           ):
    """
    summary of all metrics for one Anndata object
    """
    
    checkAdata(adata)
    checkBatch(batch_key, adata.obs)
    checkBatch(label_key, adata.obs)
    
    checkAdata(adata_int)
    checkBatch(batch_key, adata_int.obs)
    checkBatch(label_key, adata_int.obs)
    
    
    # clustering
    if nmi_ or ari_:
        print('clustering...')
        cluster_key = 'cluster'
        res_max, nmi_max, nmi_all = opt_louvain(adata_int, label_key=label_key, cluster_key=cluster_key,
                plot=False, verbose=verbose, inplace=True, force=True)
        if cluster_nmi is not None:
            nmi_all.to_csv(cluster_nmi, header=False)
            print(f'saved clustering NMI values to {cluster_nmi}')

    results = {}
    
    if nmi_:
        print('NMI...')
        nmi_score = nmi(adata_int, group1=cluster_key, group2=label_key, method=nmi_method, nmi_dir=nmi_dir)
    else:
        nmi_score = np.nan
    results['NMI_cluster/label'] = nmi_score

    if ari_:
        print('ARI...')
        ari_score = ari(adata_int, group1=cluster_key, group2=label_key)
    else:
        ari_score = np.nan
    results['ARI_cluster/label'] = ari_score
    
    if silhouette_:
        print('silhouette score...')
        # global silhouette coefficient
        sil_global = silhouette(adata_int, group_key=label_key, embed=embed)
        # silhouette coefficient per batch
        _, sil_clus = silhouette_batch(adata_int, batch_key=batch_key, group_key=label_key,
                embed=embed, verbose=False)
        sil_clus = sil_clus['silhouette_score'].mean()
    else:
        sil_global = np.nan
        sil_clus =np.nan
    results['ASW_label'] = sil_global
    results['ASW_label/batch'] = sil_clus

    if cell_cycle_:
        print('cell cycle effect...')
        cc_score = cell_cycle(adata, adata_int, batch_key=batch_key, embed=embed, agg_func=np.mean, organism=organism)
    else:
        cc_score = np.nan
    results['cell_cycle_conservation'] = cc_score
    
    if pcr_:
        print('PC regression...')
        pcr_score = pcr_comparison(adata, adata_int, embed=embed, covariate=batch_key, verbose=verbose)
    else:
        pcr_score = np.nan
    results['PCR_batch'] = pcr_score
    
    if kBET_:
        print('kBET...')
        kbet_score = 1-np.nanmean(kBET(adata_int, batch_key=batch_key, label_key=label_key, type_ = type_,
                           subsample=kBET_sub, heuristic=True, verbose=verbose)['kBET'])
    else:
        kbet_score = np.nan
    results['kBET'] = kbet_score

    if lisi_:
        print('LISI score...')
        lisi_score = np.nanmedian(lisi(adata_int, batch_key=batch_key, label_key=label_key, type_ = type_,
                                       verbose=verbose), axis=1)
        ilisi_score = lisi_score[0] - 1 #LISI scores operate on 1 - 2 (for iLISI: 2 good, 1 bad)
        clisi_score = 2 - lisi_score[1] #LISI scores operate on 1 - 2 (for cLISI: 1 good, 2 bad)

    else:
        ilisi_score = np.nan
        clisi_score = np.nan
    results['iLISI'] = ilisi_score
    results['cLISI'] = clisi_score
    
    return pd.DataFrame.from_dict(results, orient='index')
