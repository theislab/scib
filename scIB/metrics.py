import numpy as np
from scipy import sparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata
import networkx as nx
from scIB.utils import *
from scIB.lisi_py.compute_simpson_index_graph_cy import compute_simpson_index_graph_cy #cythonized lisi
from scIB.preprocessing import score_cell_cycle
from scIB.clustering import opt_louvain
from scipy.sparse.csgraph import connected_components

import rpy2.rinterface_lib.callbacks
import logging
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR) # Ignore R warning messages
import rpy2.robjects as ro
import anndata2ri


### Silhouette score
def silhouette(adata, group_key, metric='euclidean', embed='X_pca', scale=True):
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

### Isolated label score
def isolated_labels(adata, label_key, batch_key, cluster_key="iso_cluster", 
                    cluster=True, n=None, all_=False, verbose=True):
    """
    score how well labels of isolated labels are distiguished in the dataset by
        1. clustering-based approach
        2. silhouette score
    params:
        cluster: if True, use clustering approach, otherwise use silhouette score approach
        n: max number of batches per label for label to be considered as isolated.
            if n=1, take labels that are present for a single batch
            if n=None, consider any label that is missing at least 1 batch
        all_: return scores for all isolated labels instead of aggregated mean
    return:
        by default, mean of scores for each isolated label
        retrieve dictionary of scores for each label if `all_` is specified
    """
    
    scores = {}
    isolated_labels = get_isolated_labels(adata, label_key, batch_key, cluster_key,
                                          n=n, verbose=verbose)
    for label in isolated_labels:
        score = score_isolated_label(adata, label_key, batch_key, cluster_key,
                                     label, cluster=cluster, verbose=verbose)
        scores[label] = score
    
    if all_:
        return scores
    return np.mean(list(scores.values()))

def get_isolated_labels(adata, label_key, batch_key, cluster_key, n, verbose):
    """
    get labels that are considered isolated by the number of batches
    """
    
    tmp = adata.obs[[label_key, batch_key]].drop_duplicates()
    batch_per_lab = tmp.groupby(label_key).agg({batch_key: "count"})
    
    # threshold for determining when label is considered isolated
    n_batch = adata.obs[batch_key].nunique()
    if n is None:
        n = batch_per_lab.min().tolist()[0]
    
    if verbose:
        print(f"isolated labels: no more than {n} batches per label")
    
    labels = batch_per_lab[batch_per_lab[batch_key] <= n].index.tolist()
    if len(labels) == 0 and verbose:
        print(f"no isolated labels with less than {n} batches")
    return labels

def score_isolated_label(adata, label_key, batch_key, cluster_key,
                         label, cluster=True, verbose=False, **kwargs):
    """
    compute label score for a single label
    params:
        cluster: if True, use clustering approach, otherwise use silhouette score approach
    """
    
    import sklearn.metrics as scm
    adata_tmp = adata.copy()
    
    # cluster optimizing over cluster with largest number of isolated label per batch
    def max_label_per_batch(adata, label_key, cluster_key, label, argmax=False):
        sub = adata.obs[adata.obs[label_key] == label].copy()
        if argmax:
            return sub[cluster_key].value_counts().argmax()
        return sub[cluster_key].value_counts().max()
    
    if cluster:
        opt_louvain(adata_tmp, label_key, cluster_key, function=max_label_per_batch,
                    label=label, verbose=False)
    
        largest_cluster = max_label_per_batch(adata_tmp, label_key, 
                                              cluster_key, label, argmax=True)
        y_pred = adata_tmp.obs[cluster_key] == largest_cluster
        y_true = adata_tmp.obs[label_key] == label
        score = scm.f1_score(y_pred, y_true)
    else:
        adata_tmp.obs['group'] = adata_tmp.obs[label_key] == label
        score = silhouette(adata_tmp, group_key='group', **kwargs)
    
    del adata_tmp
    
    if verbose:
        print(f"{label}: {score}")
    
    return score
    
    
### Highly Variable Genes conservation
def hvg_overlap(adata_pre, adata_post, batch, n_hvg=500):
    hvg_post = adata_post.var.index
    
    adata_pre_list = splitBatches(adata_pre, batch, hvg=hvg_post)
    adata_post_list = splitBatches(adata_post, batch)
    overlap = []
    
    for i in range(len(adata_pre_list)):#range(len(adata_pre_list)):
        sc.pp.filter_genes(adata_pre_list[i], min_cells=1) # remove genes unexpressed (otherwise hvg might break)
        sc.pp.filter_genes(adata_post_list[i], min_cells=1)
        
        ov = list(set(adata_pre_list[i].var_names).intersection(adata_post_list[i].var_names))
        adata_pre_list[i] = adata_pre_list[i][:,ov]
        adata_post_list[i] = adata_post_list[i][:,ov]
        
        n_hvg_tmp = np.minimum(n_hvg, int(0.5*adata_pre_list[i].n_vars))
        if n_hvg_tmp<n_hvg:
            print(adata_pre_list[i].obs[batch][0]+' has less than the specified number of genes')
            print('Number of genes: '+str(adata_pre_list[i].n_vars))
        hvg_pre = sc.pp.highly_variable_genes(adata_pre_list[i], flavor='cell_ranger', n_top_genes=n_hvg_tmp, inplace=False)
        tmp_pre = adata_pre_list[i].var.index[hvg_pre['highly_variable']]
        hvg_post = sc.pp.highly_variable_genes(adata_post_list[i], flavor='cell_ranger', n_top_genes=n_hvg_tmp, inplace=False)
        tmp_post = adata_post_list[i].var.index[hvg_post['highly_variable']]
        n_hvg_real = np.minimum(len(tmp_pre),len(tmp_post))
        overlap.append((len(set(tmp_pre).intersection(set(tmp_post))))/n_hvg_real)
    return np.mean(overlap)

### Cell cycle effect
def cell_cycle(adata_pre, adata_post, batch_key, embed=None, agg_func=np.mean,
               organism='mouse', n_comps=50, verbose=False):
    """
    Compare the variance contribution of S-phase and G2/M-phase cell cycle scores before and
    after integration. Cell cycle scores are computed per batch on the unintegrated data set,
    eliminatimg the batch effect confounded by the `batch_key` variable. This function
    returns a score between 1 and 0. The larger the score, the stronger the cell cycle
    variance is conserved.
    This score can be calculated on full corrected feature spaces and latent embeddings as 
    variance contributions of a fixed score can be obtained via PC regression here.
    params:
        adata_pre, adata_post: adatas before and after integration
        embed   : if `embed=None`, use the full expression matrix (`adata.X`), otherwise
                  use the embedding provided in `adata_post.obsm[embed]`
        agg_func: any function that takes a list of numbers and aggregates them into a single number. 
                  If `agg_func=None`, all results will be returned
        organism: 'mouse' or 'human' for choosing cell cycle genes
    """
    checkAdata(adata_pre)
    checkAdata(adata_post)
    
    if embed == 'X_pca':
        embed = None
    
    batches = adata_pre.obs[batch_key].unique()
    scores_final = []
    scores_before = []
    scores_after = []
    for batch in batches:
        raw_sub = adata_pre[adata_pre.obs[batch_key] == batch]
        int_sub = adata_post[adata_post.obs[batch_key] == batch]
        int_sub = int_sub.obsm[embed] if embed is not None else int_sub.X
        
        if raw_sub.shape[0] != int_sub.shape[0]:
            message = f'batch "{batch}" of batch_key "{batch_key}" '
            message += 'has unequal number of entries before and after integration.'
            message += f'before: {raw_sub.shape[0]} after: {int_sub.shape[0]}'
            raise ValueError(message)
        
        if verbose:
            print("score cell cycle")
        score_cell_cycle(raw_sub, organism=organism)
        covariate = raw_sub.obs[['S_score', 'G2M_score']]
        
        before = pc_regression(raw_sub.X, covariate, pca_sd=None, n_comps=n_comps, verbose=verbose)
        scores_before.append(before)
        
        after =  pc_regression(int_sub, covariate, pca_sd=None, n_comps=n_comps, verbose=verbose)
        scores_after.append(after)
        
        score = 1 - abs(after - before)/before # scaled result
        if score < 0:
            # Here variance contribution becomes more than twice as large as before
            if verbose:
                print("Variance contrib more than twice as large after integration.")
                print("Setting score to 0.")
            score = 0
        
        scores_final.append(score)
        
        if verbose:
            print(f"batch: {batch}\t before: {before}\t after: {after}\t score: {score}")
        
    if agg_func is None:
        return pd.DataFrame([batches, scores_before, scores_after, scores_final],
                            columns=['batch', 'before', 'after', 'score'])
    else:
        return agg_func(scores_final)

### PC Regression        
def pcr_comparison(adata_pre, adata_post, covariate, embed=None, n_comps=50, scale=True, verbose=False):
    """
    Compare the effect before and after integration
    Return either the difference of variance contribution before and after integration
    or a score between 0 and 1 (`scaled=True`) with 0 if the variance contribution hasn't 
    changed. The larger the score, the more different the variance contributions are before 
    and after integration.
    params:
        adata_pre: uncorrected adata
        adata_post: integrated adata
        embed   : if `embed=None`, use the full expression matrix (`adata.X`), otherwise
                  use the embedding provided in `adata_post.obsm[embed]`
        scale: if True, return scaled score
    return:
        difference of R2Var value of PCR
    """
    
    if embed == 'X_pca':
        embed = None
    
    pcr_before = pcr(adata_pre, covariate=covariate, recompute_pca=True,
                     n_comps=n_comps, verbose=verbose)
    pcr_after = pcr(adata_post, covariate=covariate, embed=embed, recompute_pca=True,
                    n_comps=n_comps, verbose=verbose)

    if scale:
        score = (pcr_before - pcr_after)/pcr_before
        if score < 0:
            print("Variance contribution increased after integration!")
            print("Setting PCR comparison score to 0.")
            score = 0
        return score
    else:
        return pcr_after - pcr_before

def pcr(adata, covariate, embed=None, n_comps=50, recompute_pca=True, verbose=False):
    """
    PCR for Adata object
    Checks whether to
        + compute PCA on embedding or expression data (set `embed` to name of embedding matrix e.g. `embed='X_emb'`)
        + use existing PCA (only if PCA entry exists)
        + recompute PCA on expression matrix (default)
    params:
        adata: Anndata object
        embed   : if `embed=None`, use the full expression matrix (`adata.X`), otherwise
                  use the embedding provided in `adata_post.obsm[embed]`
        n_comps: number of PCs if PCA should be computed
        covariate: key for adata.obs column to regress against
    return:
        R2Var of PCR
    """
    
    checkAdata(adata)
    checkBatch(covariate, adata.obs)
    
    if verbose:
        print(f"covariate: {covariate}")
    batch = adata.obs[covariate]
    
    # use embedding for PCA
    if (embed is not None) and (embed in adata.obsm):
        if verbose:
            print(f"compute PCR on embedding n_comps: {n_comps}")
        return pc_regression(adata.obsm[embed], batch, n_comps=n_comps)
    
    # use existing PCA computation
    elif (recompute_pca == False) and ('X_pca' in adata.obsm) and ('pca' in adata.uns):
        if verbose:
            print("using existing PCA")
        return pc_regression(adata.obsm['X_pca'], batch, pca_sd=adata.uns['pca']['variance'])
    
    # recompute PCA
    else:
        if verbose:
            print(f"compute PCA n_comps: {n_comps}")
        return pc_regression(adata.X, batch, n_comps=n_comps)

def pc_regression(data, covariate, pca_sd=None, n_comps=50, svd_solver='arpack', verbose=False):
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
            print("compute PCA")
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
        print("fit regression on PCs")
    
    # one-hot encode categorical values
    covariate = pd.get_dummies(covariate).to_numpy()
    
    # fit linear model for n_comps PCs
    from sklearn.linear_model import LinearRegression
    from sklearn import metrics as scm
    r2 = []
    for i in range(n_comps):
        pc = X_pca[:, [i]]
        lm = LinearRegression()
        lm.fit(pc, covariate)
        r2_score = lm.score(pc, covariate)
        #pred = lm.predict(pc)
        #r2_score = scm.r2_score(pred, covariate, multioutput='uniform_average')
        #print(r2_score)
        #print(pred)
        #print(covariate)
        r2.append(r2_score)
    
    Var = pca_sd**2 / sum(pca_sd**2) * 100
    R2Var = sum(r2*Var)/100
    
    return R2Var

### lisi score
def get_hvg_indices(adata, verbose=True):
    if "highly_variable" not in adata.var.columns:
        if verbose:
            print(f"No highly variable genes computed, continuing with full matrix {adata.shape}")
        return np.array(range(adata.n_vars))
    return np.where((adata.var["highly_variable"] == True))[0]

def select_hvg(adata, select=True):
    if select and 'highly_variable' in adata.var:
        return adata[:, adata.var['highly_variable']].copy()
    else:
        return adata

### diffusion for connectivites matrix extension
def diffusion_conn(adata, min_k=50, copy=True, max_iterations=16):
    '''
    This function performs graph diffusion on the connectivities matrix until a
    minimum number `min_k` of entries per row are non-zero.
    
    Note:
    Due to self-loops min_k-1 non-zero connectivies entries is actually the stopping
    criterion. This is equivalent to `sc.pp.neighbors`.
    
    Returns:
       The diffusion-enhanced connectivities matrix of a copy of the AnnData object
       with the diffusion-enhanced connectivities matrix is in 
       `adata.uns["neighbors"]["conectivities"]`
    '''
    if 'neighbors' not in adata.uns:
        raise ValueError('`neighbors` not in adata object. '
                         'Please compute a neighbourhood graph!')

    if 'connectivities' not in adata.uns['neighbors']:
        raise ValueError('`connectivities` not in `adata.uns["neighbors"]`. '
                         'Please pass an object with connectivities computed!')


    T = adata.uns['neighbors']['connectivities']

    #Normalize T with max row sum
    # Note: This keeps the matrix symmetric and ensures |M| doesn't keep growing
    T = sparse.diags(1/np.array([T.sum(1).max()]*T.shape[0]))*T
    
    M = T

    # Check for disconnected component
    n_comp, labs = connected_components(adata.uns['neighbors']['connectivities'],
                                                       connection='strong')
    
    if n_comp > 1:
        tab = pd.value_counts(labs)
        small_comps = tab.index[tab<min_k]
        large_comp_mask = np.array(~pd.Series(labs).isin(small_comps))
    else:
        large_comp_mask = np.array([True]*M.shape[0])

    T_agg = T
    i = 2
    while ((M[large_comp_mask,:][:,large_comp_mask]>0).sum(1).min() < min_k) and (i < max_iterations):
        print(f'Adding diffusion to step {i}')
        T_agg *= T
        M += T_agg
        i+=1

    if (M[large_comp_mask,:][:,large_comp_mask]>0).sum(1).min() < min_k:
        raise ValueError('could not create diffusion connectivities matrix'
                         f'with at least {min_k} non-zero entries in'
                         f'{max_iterations} iterations.\n Please increase the'
                         'value of max_iterations or reduce k_min.\n')

    M.setdiag(0)

    if copy:
        adata_tmp = adata.copy()
        adata_tmp.uns['neighbors'].update({'diffusion_connectivities': M})
        return adata_tmp

    else:
        return M

    
### diffusion neighbourhood score
def diffusion_nn(adata, k, max_iterations=16):
    '''
    This function generates a nearest neighbour list from a connectivities matrix
    as supplied by BBKNN or Conos. This allows us to select a consistent number
    of nearest neighbours across all methods.

    Return:
       `k_indices` a numpy.ndarray of the indices of the k-nearest neighbors.
    '''
    if 'neighbors' not in adata.uns:
        raise ValueError('`neighbors` not in adata object. '
                         'Please compute a neighbourhood graph!')
    
    if 'connectivities' not in adata.uns['neighbors']:
        raise ValueError('`connectivities` not in `adata.uns["neighbors"]`. '
                         'Please pass an object with connectivities computed!')
        
    T = adata.uns['neighbors']['connectivities']

    # Row-normalize T
    T = sparse.diags(1/T.sum(1).A.ravel())*T
    
    T_agg = T**3
    M = T+T**2+T_agg
    i = 4
    
    while ((M>0).sum(1).min() < (k+1)) and (i < max_iterations): 
        #note: k+1 is used as diag is non-zero (self-loops)
        print(f'Adding diffusion to step {i}')
        T_agg *= T
        M += T_agg
        i+=1

    if (M>0).sum(1).min() < (k+1):
        raise ValueError(f'could not find {k} nearest neighbors in {max_iterations}'
                         'diffusion steps.\n Please increase max_iterations or reduce'
                         ' k.\n')
    
    M.setdiag(0)
    k_indices = np.argpartition(M.A, -k, axis=1)[:, -k:]
    
    return k_indices


def lisi_knn(adata, batch_key, label_key, perplexity=None, verbose=False):
    """
    Deprecated
    Compute LISI score on kNN graph provided in the adata object. By default, perplexity
    is chosen as 1/3 * number of nearest neighbours in the knn-graph.
    """
    
    if 'neighbors' not in adata.uns:
        raise AttributeError(f"key 'neighbors' not found. Please make sure that a " +
                             "kNN graph has been computed")
    elif verbose:
        print("using precomputed kNN graph")
    
    #get knn index matrix
    if verbose:
        print("Convert nearest neighbor matrix and distances for LISI.")
    dist_mat = sparse.find(adata.uns['neighbors']['distances'])
    #get number of nearest neighbours parameter
    if 'params' not in adata.uns['neighbors']:
        #estimate the number of nearest neighbors as the median 
        #of the distance matrix
        _, e = np.unique(dist_mat[0], return_counts=True)
        n_nn = np.nanmin(e)
        n_nn = n_nn.astype('int')
    else:
        n_nn = adata.uns['neighbors']['params']['n_neighbors']-1
    nn_index = np.empty(shape=(adata.uns['neighbors']['distances'].shape[0],
                               n_nn))
    nn_dists = np.empty(shape=(adata.uns['neighbors']['distances'].shape[0],
                               n_nn))
    index_out = []
    for cell_id in np.arange(np.min(dist_mat[0]), np.max(dist_mat[0])+1):
        get_idx = dist_mat[0] == cell_id
        num_idx = get_idx.sum()
        #in case that get_idx contains more than n_nn neighbours, cut away the outlying ones
        #potential enhancement: handle case where less than n_nn neighbours are reported
        if num_idx >= n_nn:
            nn_index[cell_id,:] = dist_mat[1][get_idx][np.argsort(dist_mat[2][get_idx])][:n_nn]
            nn_dists[cell_id,:] = np.sort(dist_mat[2][get_idx])[:n_nn]
        else:
            index_out.append(cell_id)
    
    out_cells = len(index_out)
    
    if out_cells > 0:
        #remove all indexes in nn_index and nn_dists, which are 0
        #COMMENT: Terrible idea and commented out
        #nn_dists = np.delete(nn_dists, index_out, 0)
        #nn_index = np.delete(nn_index, index_out, 0)
        if verbose:
            print(f"{out_cells} had less than {n_nn} neighbors and were omitted in LISI score.")
    
    if perplexity is None:
        # use LISI default
        perplexity = np.floor(nn_index.shape[1]/3)
    
    # run LISI in R
    anndata2ri.activate()
    ro.r("library(lisi)")
    
    if verbose:
        print("importing knn-graph")  
    ro.globalenv['nn_indx'] = nn_index.astype('int').T
    ro.globalenv['nn_dst'] = nn_dists.T
    ro.globalenv['perplexity'] = perplexity
    ro.globalenv['batch'] = adata.obs[batch_key].cat.codes.values
    ro.globalenv['n_batches'] = len(np.unique(adata.obs[batch_key]))
    ro.globalenv['label'] = adata.obs[label_key].cat.codes.values
    ro.globalenv['n_labels'] = len(np.unique(adata.obs[label_key]))
    
    if verbose:
        print("LISI score estimation")
    simpson_estimate_batch = ro.r(f"simpson.estimate_batch <- compute_simpson_index(nn_dst, nn_indx, batch, n_batches, perplexity)") #batch_label_keys)")
    simpson_estimate_label = ro.r(f"simpson.estimate_label <- compute_simpson_index(nn_dst, nn_indx, label, n_labels, perplexity)") #batch_label_keys)")
    simpson_est_batch = 1/np.squeeze(ro.r("simpson.estimate_batch"))
    simpson_est_label = 1/np.squeeze(ro.r("simpson.estimate_label"))
    
    anndata2ri.deactivate()
    
    # extract results
    d = {batch_key : simpson_est_batch, label_key : simpson_est_label}
    lisi_estimate = pd.DataFrame(data=d, index=np.arange(0,len(simpson_est_label)))
    
    return lisi_estimate

#LISI core functions (which we want to implement in cython for speed 
def Hbeta(D_row, beta):
    """
    Helper function for simpson index computation
    """
    P = np.exp(- D_row * beta)
    sumP = np.nansum(P)
    if (sumP == 0):
        H = 0    
        P = np.zeros(len(D_row))
    else:
        H = np.log(sumP) + beta * np.nansum(D_row*P) / sumP
        P /= sumP
    return H, P

#helper function for LISI
def convertToOneHot(vector, num_classes=None):
    """
    Converts an input 1-D vector of integers into an output
    2-D array of one-hot vectors, where an i'th input value
    of j will set a '1' in the i'th row, j'th column of the
    output array.

    Example:
        v = np.array((1, 0, 4))
        one_hot_v = convertToOneHot(v)
        print one_hot_v

        [[0 1 0 0 0]
         [1 0 0 0 0]
         [0 0 0 0 1]]
    """

    #assert isinstance(vector, np.ndarray)
    #assert len(vector) > 0

    if num_classes is None:
        num_classes = np.max(vector)+1
    #else:
    #    assert num_classes > 0
    #    assert num_classes >= np.max(vector)

    result = np.zeros(shape=(len(vector), num_classes))
    result[np.arange(len(vector)), vector] = 1
    return result.astype(int)

#LISI core functions (which we want to implement in cython for speed 
def compute_simpson_index(D = None, knn_idx = None, batch_labels = None, n_batches = None,
                          perplexity = 15, tol = 1e-5): 
    """
    Simpson index of batch labels subsetted for each group.
    params:
        D: distance matrix n_cells x n_nearest_neighbors
        knn_idx: index of n_nearest_neighbors of each cell
        batch_labels: a vector of length n_cells with batch info
        n_batches: number of unique batch labels 
        perplexity: effective neighborhood size
        tol: a tolerance for testing effective neighborhood size
    returns:
        simpson: the simpson index for the neighborhood of each cell
    """
    n = D.shape[0]
    P = np.zeros(D.shape[1])
    simpson = np.zeros(n)
    logU = np.log(perplexity)
    
    #loop over all cells
    for i in np.arange(0, n, 1):
        beta = 1
        # negative infinity
        betamin = -np.inf
        # positive infinity
        betamax = np.inf
        #get active row of D
        D_act = D[i,:]
        H, P = Hbeta(D_act, beta)
        Hdiff = H - logU
        tries = 0
        #first get neighbor probabilities
        while (np.logical_and(np.abs(Hdiff) > tol, tries < 50)):
            if (Hdiff > 0):
                betamin = beta
                if (betamax == np.inf): 
                    beta *= 2
                else:
                    beta = (beta + betamax) / 2
            else:
                betamax = beta
                if (betamin== -np.inf):
                    beta /= 2
                else:
                    beta = (beta + betamin) / 2
    
            H, P = Hbeta(D_act, beta)
            Hdiff = H - logU
            tries += 1 
        
        if (H == 0):
            simpson[i] = -1
            continue        
    
        #then compute Simpson's Index
        non_nan_knn = knn_idx[i][np.invert(np.isnan(knn_idx[i]))].astype('int')
        batch = batch_labels[non_nan_knn] 
        #convertToOneHot omits all nan entries. 
        #Therefore, we run into errors in np.matmul.
        if len(batch) == len(P):
            B = convertToOneHot(batch, n_batches)
            sumP = np.matmul(P,B) #sum P per batch
            simpson[i] = np.dot(sumP, sumP) #sum squares
        else: #assign worst possible score 
            simpson[i] = 1      
  
    return simpson


def lisi_knn_py(adata, batch_key, label_key, perplexity=None, verbose=False):
    """
    Compute LISI score on kNN graph provided in the adata object. By default, perplexity
    is chosen as 1/3 * number of nearest neighbours in the knn-graph.
    """
    
    if 'neighbors' not in adata.uns:
        raise AttributeError(f"key 'neighbors' not found. Please make sure that a " +
                             "kNN graph has been computed")
    elif verbose:
        print("using precomputed kNN graph")
    
    #get knn index matrix
    if verbose:
        print("Convert nearest neighbor matrix and distances for LISI.")
    dist_mat = sparse.find(adata.uns['neighbors']['distances'])
    #get number of nearest neighbours parameter
    if 'params' not in adata.uns['neighbors']:
        #estimate the number of nearest neighbors as the median 
        #of the distance matrix
        _, e = np.unique(dist_mat[0], return_counts=True)
        n_nn = np.nanmedian(e)
        n_nn = n_nn.astype('int')
    else:
        n_nn = adata.uns['neighbors']['params']['n_neighbors']-1
    #initialise index and fill it with NaN values
    nn_index = np.empty(shape=(adata.uns['neighbors']['distances'].shape[0],
                               n_nn))
    nn_index[:] = np.NaN
    nn_dists = np.empty(shape=(adata.uns['neighbors']['distances'].shape[0],
                               n_nn))
    nn_dists[:] = np.NaN
    index_out = []
    for cell_id in np.arange(np.min(dist_mat[0]), np.max(dist_mat[0])+1):
        get_idx = dist_mat[0] == cell_id
        num_idx = get_idx.sum()
        #in case that get_idx contains more than n_nn neighbours, cut away the outlying ones
        fin_idx = np.min([num_idx, n_nn])
        nn_index[cell_id,:fin_idx] = dist_mat[1][get_idx][np.argsort(dist_mat[2][get_idx])][:fin_idx]
        nn_dists[cell_id,:fin_idx] = np.sort(dist_mat[2][get_idx])[:fin_idx]
        if num_idx < n_nn:
            index_out.append(cell_id)
    
    out_cells = len(index_out)
    
    if out_cells > 0:
        if verbose:
            print(f"{out_cells} had less than {n_nn} neighbors.")
    
    if perplexity is None:
        # use LISI default
        perplexity = np.floor(nn_index.shape[1]/3)
    
    # run LISI in python
    if verbose:
        print("importing knn-graph")  
        
    batch = adata.obs[batch_key].cat.codes.values
    n_batches = len(np.unique(adata.obs[batch_key]))
    label = adata.obs[label_key].cat.codes.values
    n_labels = len(np.unique(adata.obs[label_key]))  
    
    if verbose:
        print("LISI score estimation")
    
    simpson_estimate_batch = compute_simpson_index(D = nn_dists, 
                                                   knn_idx = nn_index,
                                                   batch_labels = batch,                           
                                                   n_batches = n_batches,
                                                   perplexity = perplexity, 
                                                   )
    simpson_estimate_label = compute_simpson_index(D = nn_dists, 
                                                   knn_idx = nn_index,
                                                   batch_labels = label,
                                                   n_batches = n_labels,
                                                   perplexity = perplexity
                                                   )
    simpson_est_batch = 1/simpson_estimate_batch
    simpson_est_label = 1/simpson_estimate_label
    # extract results
    d = {batch_key : simpson_est_batch, label_key : simpson_est_label}
    lisi_estimate = pd.DataFrame(data=d, index=np.arange(0,len(simpson_est_label)))
    
    return lisi_estimate

    
def lisi_matrix(adata, batch_key, label_key, matrix=None, verbose=False):
    """
    deprecated
    Computes the LISI scores for a given data matrix in adata.X. The scoring function of the 
    LISI R package is called with default parameters. This function takes a data matrix and
    recomputes nearest neighbours.
    """
    
    if matrix is None:
        matrix = adata.X
    
    #lisi score runs only on dense matrices (knn search)
    if sparse.issparse(matrix):
        matrix = matrix.todense()
    
    # run LISI in R
    anndata2ri.activate()
    ro.r("library(lisi)")
    
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
    anndata2ri.deactivate()
    
    return lisi_estimate

def lisi(adata, batch_key, label_key, k0=90, type_= None, scale=True, verbose=False):
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
    
    
    #if type_ != 'knn':
    #    if verbose: 
    #        print("recompute kNN graph with {k0} nearest neighbors.")
    #recompute neighbours
    if (type_ == 'embed'):
        adata_tmp = sc.pp.neighbors(adata,n_neighbors=k0, use_rep = 'X_emb', copy=True)
    elif (type_ == 'full'):
        if 'X_pca' not in adata.obsm.keys():
            sc.pp.pca(adata, svd_solver = 'arpack')
        adata_tmp = sc.pp.neighbors(adata, n_neighbors=k0, copy=True)
    else:
        adata_tmp = adata.copy()
    #if knn - do not compute a new neighbourhood graph (it exists already)
    
    #lisi_score = lisi_knn(adata=adata, batch_key=batch_key, label_key=label_key, verbose=verbose)
    lisi_score = lisi_knn_py(adata=adata_tmp, batch_key=batch_key, label_key=label_key, verbose=verbose)
    
    # iLISI: 2 good, 1 bad
    ilisi_score = np.nanmedian(lisi_score[batch_key])
    # cLISI: 1 good, 2 bad
    clisi_score = np.nanmedian(lisi_score[label_key])
    
    if scale:
        #Comment: Scaling should be applied at the end when all scenarios are rated 
        ilisi_score = ilisi_score - 1
        #scale clisi score to 0 bad 1 good
        clisi_score = 2 - clisi_score
    
    return ilisi_score, clisi_score

#LISI core function for shortest paths 
def compute_simpson_index_graph(D = None, batch_labels = None, n_batches = None, n_neighbors = 90,
                                  perplexity = 30, subsample = None, n_chunks = 10, chunk_no = 1,tol = 1e-5, 
                                verbose = False):
    """
    Simpson index of batch labels subsetted for each group.
    params:
        D: graph object
        batch_labels: a vector of length n_cells with batch info
        n_batches: number of unique batch labels
        n_neighbors: number of nearest neighbors
        perplexity: effective neighborhood size
        tol: a tolerance for testing effective neighborhood size
    returns:
        simpson: the simpson index for the neighborhood of each cell
    """
    #compute shortest paths of everyone to everyone
    #Update: We don't actually need that, because we can compute 
    #the distance from one to all others when we actually need it
    #dist = nx.all_pairs_dijkstra_path_length(D)
    
    n = len(batch_labels)
    P = np.zeros(n_neighbors)
    logU = np.log(perplexity)
    
    #prepare chunk
    if n_chunks is not None:
        n_ch = n_chunks #number of chunks
        #get start and endpoint of chunk
        bounds = np.arange(0,n, np.ceil(n/n_ch).astype('int'))
        if chunk_no < n_ch - 1:
            chunk_ids = np.arange(bounds[chunk_no], bounds[chunk_no+1])
            if verbose:
                print(f"Entering chunk {chunk_no}.")
        else: #last chunk
            chunk_ids = np.arange(bounds[chunk_no], n)
            if verbose:
                print("Entering last chunk.")
    else:
        chunk_ids = np.arange(0, n)
        
    #remove chunk_ids, which are not in subsample
    if subsample is not None:
        chunk_ids = chunk_ids[np.in1d(chunk_ids, subsample)]
        
    simpson = np.zeros(len(chunk_ids))
    #chunk has a start and an end
    #if chunk_ids[0] != 0:
    #    consume(dist, chunk_ids[0]) #fast forward to first element of chunk
    
    #loop over all cells in chunk number
    for i in enumerate(chunk_ids): 
        #get neighbors and distances
        res = nx.single_source_dijkstra_path_length(D, i[1])
        if len(res)<n_neighbors:
            #not enough neighbors
            simpson[i[0]] = 1 # np.nan #set nan for testing
            continue
        #get sorted list of neighbours (keys) and distances (values)
        keys = np.array(list(res.keys()))
        values = np.array(list(res.values()))
        
        #start lisi estimation
        beta = 1
        # negative infinity
        betamin = -np.inf
        # positive infinity
        betamax = np.inf
        #set distances
        D_act = values[1:][:n_neighbors]
        H, P = Hbeta(D_act, beta)
        Hdiff = H - logU
        tries = 0
        #first get neighbor probabilities
        while (np.logical_and(np.abs(Hdiff) > tol, tries < 50)):
            if (Hdiff > 0):
                betamin = beta
                if (betamax == np.inf): 
                    beta *= 2
                else:
                    beta = (beta + betamax) / 2
            else:
                betamax = beta
                if (betamin == -np.inf):
                    beta /= 2
                else:
                    beta = (beta + betamin) / 2
                
            H, P = Hbeta(D_act, beta)
            Hdiff = H - logU
            tries += 1
        
        if (H == 0):
            simpson[i[0]] = -1
            continue        
        #then compute Simpson's Index
        knn_idx = keys[1:][:n_neighbors]
        batch = batch_labels[knn_idx] 
        B = convertToOneHot(batch, n_batches)
        sumP = np.matmul(P,B) #sum P per batch
        simpson[i[0]] = np.dot(sumP, sumP) #sum squares
        
    return simpson

#function to prepare call of compute_simpson_index
def lisi_graph_py(adata, batch_key, n_neighbors = 90, perplexity=None, subsample = None, 
                  multiprocessing = None, nodes = None, verbose=False):
    """
    Compute LISI score on shortes path based on kNN graph provided in the adata object. 
    By default, perplexity is chosen as 1/3 * number of nearest neighbours in the knn-graph.
    """
    
    if 'neighbors' not in adata.uns:
        raise AttributeError(f"key 'neighbors' not found. Please make sure that a " +
                              "kNN graph has been computed")    
    elif verbose:                                                    
        print("using precomputed kNN graph")
                                                                        
    #get knn index matrix
    if verbose:
        print("Convert nearest neighbor matrix and distances for LISI.")
               
    batch = adata.obs[batch_key].cat.codes.values        
    n_batches = len(np.unique(adata.obs[batch_key])) 
                                        
    if perplexity is None or perplexity >=n_neighbors:
        # use LISI default
        perplexity = np.floor(n_neighbors/3)
                                                                                                                                
    # run LISI in python
    if verbose:
        print("Compute shortest paths") 
                                                                                                
    #turn connectivities matrix into graph
    G = nx.from_scipy_sparse_matrix(adata.uns['neighbors']['connectivities'])  
          
    if verbose:
        print("LISI score estimation")
    
    #do the simpson call 
    if multiprocessing is not None:
        #import tools needed for multiprocessing
        import itertools
        from multiprocessing import Pool
        import multiprocessing
        
        #set up multiprocessing
        if nodes is None:
            #take all but one CPU and 1 CPU, if there's only 1 CPU.
            n_cpu = multiprocessing.cpu_count()
            n_processes = np.max([ n_cpu, 
                               np.ceil(n_cpu/2)]).astype('int')
        else:
            n_processes = nodes

        if verbose:
            print(f"{n_processes} processes started.")
        pool = Pool(processes=n_processes)
        count = np.arange(0, n_processes)
        
        #create argument list for each worker
        results = pool.starmap(compute_simpson_index_graph, zip(itertools.repeat(G),
                                                                itertools.repeat(batch),
                                                                itertools.repeat(n_batches),
                                                                itertools.repeat(n_neighbors),
                                                                itertools.repeat(perplexity),
                                                                itertools.repeat(subsample),
                                                                itertools.repeat(n_processes),
                                                                count))
        pool.close()
        pool.join()
        
        simpson_est_batch = 1/np.concatenate(results)    
     
    else: 
        simpson_estimate_batch = compute_simpson_index_graph(D = G, 
                                                  batch_labels = batch,                           
                                                  n_batches = n_batches,
                                                  perplexity = perplexity, 
                                                  subsample = subsample,
                                                  n_neighbors = n_neighbors,
                                                  n_chunks = None,
                                                  chunk_no = None,
                                                  verbose = verbose
                                                 )
        simpson_est_batch = 1/simpson_estimate_batch
    # extract results
    d = {batch_key : simpson_est_batch}
    if subsample is None:
        lisi_estimate = pd.DataFrame(data=d, index=np.arange(0,len(simpson_est_batch)))
    else:
        lisi_estimate = pd.DataFrame(data=d, index=np.sort(subsample))
    
    return lisi_estimate

#function to prepare call of compute_simpson_index_graph_cy (cythonized version of lisi_graph_py)
def lisi_graph_cy(adata, batch_key, n_neighbors = 90, perplexity=None, subsample = None, 
                  multiprocessing = None, nodes = None, verbose=False):
    """
    Compute LISI score on shortes path based on kNN graph provided in the adata object. 
    By default, perplexity is chosen as 1/3 * number of nearest neighbours in the knn-graph.
    """
    
    if 'neighbors' not in adata.uns:
        raise AttributeError(f"key 'neighbors' not found. Please make sure that a " +
                              "kNN graph has been computed")    
    elif verbose:                                                    
        print("using precomputed kNN graph")
                                                                        
    #get knn index matrix
    if verbose:
        print("Convert nearest neighbor matrix and distances for LISI.")
               
    batch = adata.obs[batch_key].cat.codes.values.astype('int')         
    n_batches = len(np.unique(adata.obs[batch_key])) 
                                        
    if perplexity is None or perplexity >=n_neighbors:
        # use LISI default
        perplexity = np.floor(n_neighbors/3)
    
    if subsample is None:
        subsample = np.ndarray(shape=0, dtype=np.long)

    # run LISI in python
    #if verbose:
    #    print("Compute shortest paths") 
                                                                                                
    #turn connectivities matrix into graph
    #G = nx.from_scipy_sparse_matrix(adata.uns['neighbors']['connectivities'])  
    G_data = adata.uns['neighbors']['connectivities'].data
    G_indices = adata.uns['neighbors']['connectivities'].indices
    G_indptr = adata.uns['neighbors']['connectivities'].indptr
          
    if verbose:
        print("LISI score estimation")
    
    #do the simpson call 
    if multiprocessing is not None:
        #import tools needed for multiprocessing
        import itertools
        from multiprocessing import Pool
        import multiprocessing
        
        #set up multiprocessing
        if nodes is None:
            #take all but one CPU and 1 CPU, if there's only 1 CPU.
            n_cpu = multiprocessing.cpu_count()
            n_processes = np.max([ n_cpu, 
                               np.ceil(n_cpu/2)]).astype('int')
        else:
            n_processes = nodes

        if verbose:
            print(f"{n_processes} processes started.")
        pool = Pool(processes=n_processes)
        count = np.arange(0, n_processes)
        
        #create argument list for each worker
        results = pool.starmap(compute_simpson_index_graph_cy, zip(itertools.repeat(G_data),
                                                                   itertools.repeat(G_indices),
                                                                   itertools.repeat(G_indptr),
                                                                itertools.repeat(batch),
                                                                itertools.repeat(n_batches),
                                                                itertools.repeat(n_neighbors),
                                                                itertools.repeat(perplexity),
                                                                itertools.repeat(subsample),
                                                                itertools.repeat(n_processes),
                                                                count, 
                                                                   itertools.repeat(1e-5)))
        pool.close()
        pool.join()
        
        simpson_est_batch = 1/np.concatenate(results)    
     
    else: 
        simpson_estimate_batch = compute_simpson_index_graph_cy(data= G_data, 
                                                                indices = G_indices,
                                                                indptr = G_indptr,
                                                  batch_labels = batch,                           
                                                  n_batches = n_batches,
                                                  perplexity = perplexity, 
                                                  subsample = subsample,
                                                  n_neighbors = n_neighbors,
                                                  n_chunks = 1,
                                                  chunk_no = 0,
                                                  tol= 1e-5
                                                 )
        simpson_est_batch = 1/simpson_estimate_batch
    # extract results
    d = {batch_key : simpson_est_batch}
    if len(subsample)==0:
        lisi_estimate = pd.DataFrame(data=d, index=np.arange(0,len(simpson_est_batch)))
    else:
        lisi_estimate = pd.DataFrame(data=d, index=np.sort(subsample))
    
    return lisi_estimate

#LISI graph function (analoguous to lisi function) 
def lisi_graph(adata, batch_key=None, label_key=None, k0=90, type_= None, 
               subsample = None, scale=True, 
               multiprocessing = None, nodes = None, verbose=False):
    """
    Compute lisi score (after integration)
    params:
        adata: adata object to calculate on
        batch_key: variable to compute iLISI on
        label_key: variable to compute cLISI on
        k0: number of nearest neighbors to compute lisi score
            Please note that the initial neighborhood size that is
            used to compute shortest paths is 15.
        type_: type of data integration, either knn, full or embed
        subsample: Fraction of observations (between 0 and 1) 
                   to which lisi scoring should be subsampled
        scale: scale output values (True/False)
        multiprocessing: parallel computation of LISI scores, if None, no parallisation 
                         via multiprocessing is performed
        nodes: number of nodes (i.e. CPUs to use for multiprocessing); ignored, if
               multiprocessing is set to None
    return:
        pd.DataFrame with median cLISI and median iLISI scores 
        (following the harmony paper)
    """
    
    checkAdata(adata)
    checkBatch(batch_key, adata.obs)
    checkBatch(label_key, adata.obs)
        
    #recompute neighbours
    if (type_ == 'embed'):
        adata_tmp = sc.pp.neighbors(adata,n_neighbors=15, use_rep = 'X_emb', copy=True)
    if (type_ == 'full'):
        if 'X_pca' not in adata.obsm.keys():
            sc.pp.pca(adata, svd_solver = 'arpack')
        adata_tmp = sc.pp.neighbors(adata, n_neighbors=15, copy=True)
    else:
        adata_tmp = adata.copy()
    #if knn - do not compute a new neighbourhood graph (it exists already)
    
    if subsample is not None:
        #check if subsample is indeed 1 float number between 0 and 1
        if (not isinstance(subsample,float)) or np.logical_or(subsample<0, subsample>1):
            raise ValueError('`subsample` not a fraction between 0 and 1 or has wrong size.')
        subset = np.random.choice(np.arange(0,adata_tmp.n_obs), 
                     np.floor(subsample*adata_tmp.n_obs).astype('int'),
                     replace=False
                     )
    else:
        subset = None
    
    #compute LISI score
    ilisi_score = lisi_graph_cy(adata = adata, batch_key = batch_key, 
                  n_neighbors = k0, perplexity=None, subsample = subset, 
                  multiprocessing = multiprocessing, nodes = nodes, verbose=verbose)
    
    clisi_score = lisi_graph_cy(adata = adata, batch_key = label_key, 
                  n_neighbors = k0, perplexity=None, subsample = subset, 
                  multiprocessing = multiprocessing, nodes = nodes, verbose=verbose)
    
    # iLISI: 2 good, 1 bad
    ilisi_score = np.nanmedian(ilisi_score)
    # cLISI: 1 good, 2 bad
    clisi_score = np.nanmedian(clisi_score)
    
    if scale:
        #Comment: Scaling should be applied at the end when all scenarios are rated 
        ilisi_score = ilisi_score - 1
        #scale clisi score to 0 bad 1 good
        clisi_score = 2 - clisi_score
    
    return ilisi_score, clisi_score


### kBET
def kBET_single(matrix, batch, type_ = None, k0 = 10, knn=None, subsample=0.5, heuristic=True, verbose=False):
    """
    params:
        matrix: expression matrix (at the moment: a PCA matrix, so do.pca is set to FALSE
        batch: series or list of batch assignemnts
        subsample: fraction to be subsampled. No subsampling if `subsample=None`
    returns:
        kBET p-value
    """
        
    anndata2ri.activate()
    ro.r("library(kBET)")
    
    if verbose:
        print("importing expression matrix")
    ro.globalenv['data_mtrx'] = matrix
    ro.globalenv['batch'] = batch
    #print(matrix.shape)
    #print(len(batch))
    
    if verbose:
        print("kBET estimation")
    #k0 = len(batch) if len(batch) < 50 else 'NULL'
    
    ro.globalenv['knn_graph'] = knn
    ro.globalenv['k0'] = k0
    batch_estimate = ro.r(f"batch.estimate <- kBET(data_mtrx, batch, knn=knn_graph, k0=k0, plot=FALSE, do.pca=FALSE, heuristic=FALSE, adapt=FALSE, verbose={str(verbose).upper()})")
            
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
    #compute connectivities for non-knn type data integrations
    #and increase neighborhoods for knn type data integrations
    if type_ != 'knn':
        adata_tmp = sc.pp.neighbors(adata, n_neighbors = 50, use_rep=embed, copy=True)
    else:
        #check if pre-computed neighbours are stored in input file
        adata_tmp = adata.copy()
        if 'diffusion_connectivities' not in adata.uns['neighbors']:
            if verbose:
                print(f"Compute: Diffusion neighbours.")
            adata_tmp = diffusion_conn(adata, min_k = 50, copy = True)
        adata_tmp.uns['neighbors']['connectivities'] = adata_tmp.uns['neighbors']['diffusion_connectivities']
            
    if verbose:
        print(f"batch: {batch_key}")
        
    #set upper bound for k0
    size_max = 2**31 - 1
    
    kBET_scores = {'cluster': [], 'kBET': []}
    for clus in adata_tmp.obs[label_key].unique():
        
        adata_sub = adata_tmp[adata_tmp.obs[label_key] == clus,:].copy()
        #check if neighborhood size too small or only one batch in subset
        if np.logical_or(adata_sub.n_obs < 10, 
                         len(adata_sub.obs[batch_key].cat.categories)==1):
            print(f"{clus} consists of a single batch or is too small. Skip.")
            score = np.nan
        else:
            quarter_mean = np.floor(np.mean(adata_sub.obs[batch_key].value_counts())/4).astype('int')
            k0 = np.min([70, np.max([10, quarter_mean])])
            #check k0 for reasonability
            if (k0*adata_sub.n_obs) >=size_max:
                k0 = np.floor(size_max/adata_sub.n_obs).astype('int')
           
            matrix = np.zeros(shape=(adata_sub.n_obs, k0+1))
                
            if verbose:
                print(f"Use {k0} nearest neighbors.")
            n_comp, labs = connected_components(adata_sub.uns['neighbors']['connectivities'], 
                                                              connection='strong')
            if n_comp > 1:
                #check the number of components where kBET can be computed upon
                comp_size = pd.value_counts(labs)
                #check which components are small
                comp_size_thresh = 3*k0
                idx_nonan = np.flatnonzero(np.in1d(labs, 
                                                   comp_size[comp_size>=comp_size_thresh].index))
                #check if 75% of all cells can be used for kBET run
                if len(idx_nonan)/len(labs) >= 0.75:
                    #create another subset of components, assume they are not visited in a diffusion process
                    adata_sub_sub = adata_sub[idx_nonan,:].copy()
                    nn_index_tmp = np.empty(shape=(adata_sub.n_obs, k0))
                    nn_index_tmp[:] = np.nan
                    nn_index_tmp[idx_nonan] = diffusion_nn(adata_sub_sub, k=k0).astype('float') 
                    #need to check neighbors (k0 or k0-1) as input?   
                    score = kBET_single(
                            matrix=matrix,
                            batch=adata_sub.obs[batch_key],
                            knn = nn_index_tmp+1, #nn_index in python is 0-based and 1-based in R
                            subsample=subsample,
                            verbose=verbose,
                            heuristic=False,
                            k0 = k0,
                            type_ = type_
                            )
                else:
                    #if there are too many too small connected components, set kBET score to 1 
                    #(i.e. 100% rejection)
                    score = 1
                
            else: #a single component to compute kBET on 
                #need to check neighbors (k0 or k0-1) as input?  
                nn_index_tmp = diffusion_nn(adata_sub, k=k0).astype('float')
                score = kBET_single(
                            matrix=matrix,
                            batch=adata_sub.obs[batch_key],
                            knn = nn_index_tmp+1, #nn_index in python is 0-based and 1-based in R
                            subsample=subsample,
                            verbose=verbose,
                            heuristic=False,
                            k0 = k0,
                            type_ = type_
                            )
        
        kBET_scores['cluster'].append(clus)
        kBET_scores['kBET'].append(score)
    
    kBET_scores = pd.DataFrame.from_dict(kBET_scores)
    kBET_scores = kBET_scores.reset_index(drop=True)
    
    return kBET_scores

# determine root cell for trajectory conservation metric
def get_root(adata_pre, adata_post, ct_key, dpt_dim=3):
    n_components, adata_post.obs['neighborhood'] = connected_components(csgraph=adata_post.uns['neighbors']['connectivities'], directed=False, return_labels=True)
    
    start_clust = adata_pre.obs.groupby([ct_key]).mean()['dpt_pseudotime'].idxmin()
    min_dpt = np.flatnonzero(adata_pre.obs[ct_key] == start_clust)
    max_neigh = np.flatnonzero(adata_post.obs['neighborhood']== adata_post.obs['neighborhood'].value_counts().argmax())
    min_dpt = [value for value in min_dpt if value in max_neigh]
    
    # compute Diffmap for adata_post
    sc.tl.diffmap(adata_post)
    
    # determine most extreme cell in adata_post Diffmap
    min_dpt_cell = np.zeros(len(min_dpt))
    for dim in np.arange(dpt_dim):
        
        diffmap_mean = adata_post.obsm["X_diffmap"][:, dim].mean()
        diffmap_min_dpt = adata_post.obsm["X_diffmap"][min_dpt, dim]
        
        # choose optimum function
        if diffmap_min_dpt.mean() < diffmap_mean:
            opt = np.argmin
        else:
            opt = np.argmax
        # count opt cell
        min_dpt_cell[opt(diffmap_min_dpt)] += 1
    
    # root cell is cell with max vote
    return min_dpt[np.argmax(min_dpt_cell)]

def trajectory_conservation(adata_pre, adata_post, label_key):
    cell_subset = adata_pre.obs.index[adata_pre.obs["dpt_pseudotime"].notnull()]
    adata_pre_sub = adata_pre[cell_subset]
    adata_post_sub = adata_post[cell_subset]
    
    adata_post_sub.uns['iroot'] = get_root(adata_pre_sub, adata_post_sub, label_key)
    
    sc.tl.dpt(adata_post_sub)
    adata_post_sub.obs['dpt_pseudotime'][adata_post_sub.obs['dpt_pseudotime']>1]=0
    return (adata_post_sub.obs['dpt_pseudotime'].corr(adata_pre_sub.obs['dpt_pseudotime'], 'spearman')+1)/2


def graph_connectivity(adata_post, label_key):
    """"
    Metric that quantifies how connected the subgraph corresponding to each batch cluster is.
    """
    if 'neighbors' not in adata_post.uns:
        raise KeyError('Please compute the neighborhood graph before running this '
                       'function!')

    clust_res = [] 

    for ct in adata_post.obs[label_key].cat.categories:
        adata_post_sub = adata_post[adata_post.obs[label_key].isin([ct]),]
        _,labs = connected_components(adata_post_sub.uns['neighbors']['connectivities'], connection='strong')
        tab = pd.value_counts(labs)
        clust_res.append(tab[0]/sum(tab))

    return np.mean(clust_res)


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
            hvg_score_=True, cluster_nmi=None,
            nmi_=False, ari_=False, nmi_method='arithmetic', nmi_dir=None, 
            silhouette_=False,  embed='X_pca', si_metric='euclidean',
            pcr_=False, cell_cycle_=False, organism='mouse', verbose=False,

            isolated_labels_=False, n_isolated=None, graph_conn_=False,
            kBET_=False, kBET_sub=0.5, lisi_graph_=False,
            trajectory_= False, type_ = None
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
        res_max, nmi_max, nmi_all = opt_louvain(adata_int,
                label_key=label_key, cluster_key=cluster_key, function=nmi,
                plot=False, verbose=verbose, inplace=True, force=True)
        if cluster_nmi is not None:
            nmi_all.to_csv(cluster_nmi, header=False)
            print(f'saved clustering NMI values to {cluster_nmi}')

    results = {}
    
    if nmi_:
        print('NMI...')
        nmi_score = nmi(adata_int, group1=cluster_key, group2=label_key,
                    method=nmi_method, nmi_dir=nmi_dir)
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
        sil_clus = np.nan
    results['ASW_label'] = sil_global
    results['ASW_label/batch'] = sil_clus

    if pcr_:
        print('PC regression...')
        pcr_score = pcr_comparison(adata, adata_int, embed=embed, covariate=batch_key, verbose=verbose)
    else:
        pcr_score = np.nan
    results['PCR_batch'] = pcr_score
    
    if cell_cycle_:
        print('cell cycle effect...')
        cc_score = cell_cycle(adata, adata_int, batch_key=batch_key, embed=embed,
                           agg_func=np.mean, organism=organism)
    else:
        cc_score = np.nan
    results['cell_cycle_conservation'] = cc_score
    
    if isolated_labels_:
        print("isolated labels...")
        il_score_clus = isolated_labels(adata_int, label_key=label_key, batch_key=batch_key,
                                cluster=True, n=n_isolated, verbose=False)
        if silhouette_:
            il_score_sil = isolated_labels(adata_int, label_key=label_key, batch_key=batch_key,
                                           cluster=False, n=n_isolated, verbose=False)
        else:
            il_score_sil = np.nan
    else:
        il_score_clus = np.nan
        il_score_sil  = np.nan
    results['isolated_label_F1'] = il_score_clus
    results['isolated_label_silhouette'] = il_score_sil

    if graph_conn_:
        print('Graph connectivity...')
        graph_conn_score = graph_connectivity(adata_int, label_key=label_key)
    else:
        graph_conn_scores = np.nan
    results['graph_conn'] = graph_conn_score
    
    if kBET_:
        print('kBET...')
        kbet_score = 1-np.nanmean(kBET(adata_int, batch_key=batch_key, label_key=label_key, type_=type_,
                                       embed = embed, subsample=kBET_sub, 
                                       heuristic=True, verbose=verbose)['kBET'])
    else: 
        kbet_score = np.nan
    results['kBET'] = kbet_score

    #if lisi_:
    #    print('LISI score...')
    #    ilisi_score, clisi_score = lisi(adata_int, batch_key=batch_key, label_key=label_key,
    #                                    type_ = type_, verbose=verbose)
    #else:
    #    ilisi_score = np.nan
    #    clisi_score = np.nan
    #results['iLISI'] = ilisi_score
    #results['cLISI'] = clisi_score
    
    if lisi_graph_:
        print('LISI graph score...')
        ilisi_g_score, clisi_g_score = lisi_graph(adata_int, batch_key=batch_key, label_key=label_key,
                                        type_ = type_, subsample = kBET_sub, 
                                        multiprocessing = True, verbose=verbose)
    else:
        ilisi_g_score = np.nan
        clisi_g_score = np.nan
    results['iLISI'] = ilisi_g_score
    results['cLISI'] = clisi_g_score
        
    if hvg_score_:
        hvg_score = hvg_overlap(adata, adata_int, batch_key)
    else:
        hvg_score = np.nan
    results['hvg_overlap'] = hvg_score
    
    if trajectory_:
        print('Trajectory conservation score...')
        trajectory_score = trajectory_conservation(adata, adata_int, label_key=label_key)
    else:
        trajectory_score = np.nan
    results['trajectory'] = trajectory_score
    
    return pd.DataFrame.from_dict(results, orient='index')
