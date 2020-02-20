import numpy as np
from scipy import sparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata
from scIB.utils import *
from scIB.preprocessing import score_cell_cycle
from scIB.clustering import opt_louvain

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

def lisi_knn(adata, batch_key, label_key, perplexity=None, verbose=False):
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
        nn_dists = np.delete(nn_dists, index_out, 0)
        nn_index = np.delete(nn_index, index_out, 0)
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
    if out_cells > 0:  
        batch_adapt = np.delete(adata.obs[batch_key].cat.codes.values, index_out)
        label_adapt = np.delete(adata.obs[label_key].cat.codes.values, index_out)
        ro.globalenv['batch'] = batch_adapt
        ro.globalenv['n_batches'] = len(np.unique(batch_adapt))
        ro.globalenv['label'] = label_adapt
        ro.globalenv['n_labels'] = len(np.unique(label_adapt))
    else:
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

def lisi(adata, batch_key, label_key, scale=True, verbose=False):
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
    
    lisi_score = lisi_knn(adata=adata, batch_key=batch_key, label_key=label_key, verbose=verbose)
    
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
    
    if verbose:
        print("kBET estimation")
    #k0 = len(batch) if len(batch) < 50 else 'NULL'
    if type_ == 'knn':
        ro.globalenv['knn_graph'] = knn
        ro.globalenv['k0'] = np.min([knn.shape[1], matrix.shape[0]])
        batch_estimate = ro.r(f"batch.estimate <- kBET(data_mtrx, batch, knn=knn_graph, k0=k0, plot=FALSE, do.pca=FALSE, heuristic=FALSE, adapt=FALSE, verbose={str(verbose).upper()})")
    else:
        #in this case, we do a knn search in R with FNN package
        #FNN has an upper limit for the data size it can handle
        size_max = 2**31 - 1 #limit before R uses long vector format
        #if the input matrix is potentially too large, we set an upper limit for k0
        if (matrix.shape[0]*matrix.shape[1]) >= size_max:
            ro.globalenv['k0'] = np.floor(size_max/matrix.shape[0])
            batch_estimate = ro.r(f"batch.estimate <- kBET(data_mtrx, batch, k0=k0, plot=FALSE, do.pca=FALSE, heuristic={str(heuristic).upper()}, verbose={str(verbose).upper()})")
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
        #get number of nearest neighbours parameter
        if 'params' not in adata.uns['neighbors']:
            #estimate the number of nearest neighbors as the median 
            #of the distance matrix
            _, e = np.unique(dist_mat[0], return_counts=True)
            n_nn = np.nanmedian(e)
            #set type of n_nn to int to avoid type errors downstream
            n_nn = n_nn.astype('int')
        else:
            n_nn = adata.uns['neighbors']['params']['n_neighbors']-1
        nn_index = np.empty(shape=(adata.uns['neighbors']['distances'].shape[0],
                                   n_nn))
        index_out = []
        for cell_id in np.arange(np.min(dist_mat[0]), np.max(dist_mat[0])+1):
            get_idx = dist_mat[0] == cell_id
            num_idx = get_idx.sum()
            if num_idx >= n_nn:
                nn_index[cell_id,:] = dist_mat[1][get_idx][np.argsort(dist_mat[2][get_idx])][:n_nn]
            else:
                index_out.append(cell_id)
        
        out_cells = len(index_out)
        
        if out_cells > 0:
        #remove all indexes in nn_index and nn_dists, which are 0 
            nn_index = np.delete(nn_index, index_out, 0)
            #adapt adata for the time being
            adata_tmp = adata[np.invert(np.in1d(np.arange(0, adata.shape[0]), index_out))].copy()
            if verbose:
                print(f"{out_cells} had less than {n_nn} neighbors and were omitted in kBET.")
        else:
            adata_tmp = adata.copy()
    
    matrix = adata_tmp.obsm[embed]
    
    if verbose:
        print(f"batch: {batch_key}")
    batch = adata_tmp.obs[batch_key]
    
    kBET_scores = {'cluster': [], 'kBET': []}
    for clus in adata_tmp.obs[label_key].unique():
        idx = np.where((adata_tmp.obs[label_key] == clus))[0]
        if type_ == 'knn':
            nn_index_tmp = nn_index[idx,:] #reduce nearest neighbor matrix to the desired indices
            nn_index_tmp[np.invert(np.isin(nn_index_tmp, idx))] = np.nan #set the rest nan
            score = kBET_single(
                matrix[idx, :],
                batch[idx],
                knn = nn_index_tmp+1, #nn_index in python is 0-based and 1-based in R
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

# determine root cell for trajectory conservation metric
def get_root(adata_pre, adata_post, dpt_dim=3):
    
    # minimum DPT candidate cell indices
    min_dpt = np.flatnonzero(adata_pre.obs["dpt_pseudotime"] == 0)
    #min_dpt = adata_pre.obs.index[adata_pre.obs.dpt_pseudotime == 0]
    
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

def trajectory_conservation(adata_pre, adata_post):
    cell_subset = adata_pre.obs.index[adata_pre.obs["dpt_pseudotime"].notnull()]
    adata_pre_sub = adata_pre[cell_subset]
    adata_post_sub = adata_post[cell_subset]
    
    adata_post_sub.uns['iroot'] = get_root(adata_pre_sub, adata_post_sub)
    
    sc.tl.dpt(adata_post_sub)
    adata_post_sub.obs['dpt_pseudotime'][adata_post_sub.obs['dpt_pseudotime']>1]=0
    return (adata_post_sub.obs['dpt_pseudotime'].corr(adata_pre_sub.obs['dpt_pseudotime'], 'spearman')+1)/2

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
            isolated_labels_=False, n_isolated=None,
            kBET_=False, kBET_sub=0.5, lisi_=False, trajectory_= False, type_ = None
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
    
    if kBET_:
        print('kBET...')
        kbet_score = 1-np.nanmean(kBET(adata_int, batch_key=batch_key, label_key=label_key, type_=type_,
                           subsample=kBET_sub, heuristic=True, verbose=verbose)['kBET'])
    else: 
        kbet_score = np.nan
    results['kBET'] = kbet_score

    if lisi_:
        print('LISI score...')
        ilisi_score, clisi_score = lisi(adata_int, batch_key=batch_key, label_key=label_key,
                                        verbose=verbose)
    else:
        ilisi_score = np.nan
        clisi_score = np.nan
    results['iLISI'] = ilisi_score
    results['cLISI'] = clisi_score
    
    if hvg_score_:
        hvg_score = hvg_overlap(adata, adata_int, batch_key)
    else:
        hvg_score = np.nan
    results['hvg_overlap'] = hvg_score
    
    if trajectory_:
        print('Trajectory conservation score...')
        trajectory_score = trajectory_conservation(adata, adata_int)
    else:
        trajectory_score = np.nan
    results['trajectory'] = trajectory_score
    
    return pd.DataFrame.from_dict(results, orient='index')
