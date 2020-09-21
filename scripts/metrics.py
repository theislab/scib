#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scIB
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# types of integration output
RESULT_TYPES = [
    "full", # reconstructed expression data
    "embed", # embedded/latent space
    "knn" # only corrected neighbourhood graph as output
]
ASSAYS = ["expression", "atac", "simulation"]

if __name__=='__main__':
    """
    read adata object, compute all metrics and output csv.
    """

    import argparse
    import os

    parser = argparse.ArgumentParser(description='Compute all metrics')

    parser.add_argument('-u', '--uncorrected', required=True)
    parser.add_argument('-i', '--integrated', required=True)
    parser.add_argument('-o', '--output', required=True, help='Output file')
    parser.add_argument('-m', '--method', required=True, help='Name of method')

    parser.add_argument('-b', '--batch_key', required=True, help='Key of batch')
    parser.add_argument('-l', '--label_key', required=True, help='Key of annotated labels e.g. "cell_type"')

    parser.add_argument('--organism', required=True)
    parser.add_argument('--type', required=True, choices=RESULT_TYPES, help='Type of result: full, embed, knn\n full: scanorama, seurat, MNN\n embed: scanorama, Harmony\n knn: BBKNN')
    parser.add_argument('--assay', default='expression', choices=ASSAYS, help='Experimental assay')
    parser.add_argument('--hvgs', default=0, help='Number of highly variable genes. Use 0 to specify that no feature selection had been used.', type=int)
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()

    verbose = args.verbose
    type_ = args.type
    batch_key = args.batch_key
    label_key = args.label_key
    assay = args.assay
    organism = args.organism
    n_hvgs = args.hvgs if args.hvgs > 0 else None

    # encode setup for column name
    setup = f'{args.method}_{args.type}'

    # create cluster NMI output file
    file_stump = os.path.splitext(args.output)[0]
    cluster_nmi = f'{file_stump}_nmi.txt'

    if verbose:
        print('Options')
        print(f'    type:\t{type_}')
        print(f'    batch_key:\t{batch_key}')
        print(f'    label_key:\t{label_key}')
        print(f'    assay:\t{assay}')
        print(f'    organism:\t{organism}')
        print(f'    n_hvgs:\t{n_hvgs}')
        print(f'    setup:\t{setup}')
        print(f'    optimised clustering results:\t{cluster_nmi}')

    empty_file = False

    print("reading adata before integration")
    adata = sc.read(args.uncorrected, cache=True)
    print(adata)
    print("reading adata after integration")
    if os.stat(args.integrated).st_size == 0:
        print(f'{args.integrated} is empty, setting all metrics to NA.')
        adata_int = adata
        empty_file = True
    else:
        adata_int = sc.read(args.integrated, cache=True)
        print(adata_int)

    if (n_hvgs is not None):
        if (adata_int.n_vars < n_hvgs):
            raise ValueError("There are less genes in the corrected adata than specified for HVG selection")

    # check input files
    if adata.n_obs != adata_int.n_obs:
        message = "The datasets have different numbers of cells before and after integration."
        message += "Please make sure that both datasets match."
        raise ValueError(message)

    #check if the obsnames were changed and rename them in that case
    if len(set(adata.obs_names).difference(set(adata_int.obs_names))) > 0:
        #rename adata_int.obs[batch_key] labels by overwriting them with the pre-integration labels
        new_obs_names = ['-'.join(idx.split('-')[:-1]) for idx in adata_int.obs_names]

        if len(set(adata.obs_names).difference(set(new_obs_names))) == 0:
            adata_int.obs_names = new_obs_names
        else:
            raise ValueError('obs_names changed after integration!')

    #batch_key might be overwritten, so we match it to the pre-integrated labels
    adata_int.obs[batch_key] = adata_int.obs[batch_key].astype('category')
    if not np.array_equal(adata.obs[batch_key].cat.categories,adata_int.obs[batch_key].cat.categories):
        #pandas uses the table index to match the correct labels
        adata_int.obs[batch_key] = adata.obs[batch_key]

    if (n_hvgs is not None) and (adata_int.n_vars < n_hvgs):
        # check number of HVGs to be computed
        message = "There are fewer genes in the uncorrected adata "
        message += "than specified for HVG selection."
        raise ValueError(message)

    # DATA REDUCTION
    # select options according to type

    # case 1: full expression matrix, default settings
    precompute_pca = True
    recompute_neighbors = True
    embed = 'X_pca'

    # distinguish between subsetted and full expression matrix
    # compute HVGs only if output is not already subsetted
    if adata.n_vars > adata_int.n_vars:
        n_hvgs = None

    # case 2: embedding output
    if (type_ == "embed"):
        n_hvgs = None
        embed = "X_emb"
        # legacy check
        if ('emb' in adata_int.uns) and (adata_int.uns['emb']):
            adata_int.obsm["X_emb"] = adata_int.obsm["X_pca"].copy()

    # case3: kNN graph output
    elif (type_ == "knn"):
        n_hvgs = None
        precompute_pca = False
        recompute_neighbors = False

    if verbose:
        print('reduce integrated data:')
        print(f'    HVG selection:\t{n_hvgs}')
        message = f'    compute neighbourhood graph:\t{recompute_neighbors}'
        if recompute_neighbors:
            message += f' on {embed}'
        print(message)
        print(f'    precompute PCA:\t{precompute_pca}')

    if not empty_file:
        scIB.preprocessing.reduce_data(adata_int,
                                       n_top_genes=n_hvgs,
                                       neighbors=recompute_neighbors, use_rep=embed,
                                       pca=precompute_pca, umap=False)

    print("computing metrics")
    # DEFAULT
    silhouette_ = True
    nmi_ = True
    ari_ = True
    pcr_ = True
    cell_cycle_ = True
    isolated_labels_ = True
    hvg_score_ = True
    graph_conn_ = True
    kBET_ = True
    #lisi_ = True
    lisi_graph_ = True


    # by output type
    if (type_ == "embed"):
        hvg_score_ = False
    elif (type_ == "knn"):
        silhouette_ = False
        pcr_ = False
        cell_cycle_ = False
        hvg_score_ = False
        #lisi_ = False

    # by assay
    if args.assay == 'atac':
        cell_cycle_ = False
        hvg_score_ = False
    elif args.assay == 'simulation':
        cell_cycle_ = False

    # check if pseudotime data exists in original data
    if 'dpt_pseudotime' in adata.obs:
        trajectory_ = True
    else:
        trajectory_ = False

    if empty_file:
        silhouette_=False
        nmi_=False
        ari_=False
        pcr_=False
        cell_cycle_=False
        isolated_labels_=False
        hvg_score_=False
        graph_conn_=False
        kBET_=False
        #lisi_=False
        lisi_graph_=False
        trajectory_=False

    if adata.n_obs > 300000:
        kBET_=False

    if verbose:
        print(f'type:\t{type_}')
        print(f'    ASW:\t{silhouette_}')
        print(f'    NMI:\t{nmi_}')
        print(f'    ARI:\t{ari_}')
        print(f'    PCR:\t{pcr_}')
        print(f'    cell cycle:\t{cell_cycle_}')
        print(f'    iso lab F1:\t{isolated_labels_}')
        print(f'    iso lab ASW:\t{isolated_labels_ & silhouette_}')
        print(f'    HVGs:\t{hvg_score_}')
        print(f'    kBET:\t{kBET_}')
        #print(f'    LISI:\t{lisi_}')
        print(f'    LISI:\t{lisi_graph_}')
        print(f'    Trajectory:\t{trajectory_}')

    results = scIB.me.metrics(adata, adata_int, verbose=verbose,
                              hvg_score_=hvg_score_, cluster_nmi=cluster_nmi,
                              batch_key=batch_key, label_key=label_key,
                              silhouette_=silhouette_, embed=embed,
                              type_ = type_,
                              nmi_=nmi_, nmi_method='arithmetic', nmi_dir=None,
                              ari_=ari_,
                              pcr_=pcr_,
                              cell_cycle_=cell_cycle_, organism=organism,
                              isolated_labels_=isolated_labels_, n_isolated=None,
                              graph_conn_=graph_conn_,
                              kBET_=kBET_,
                              #lisi_=lisi_,
                              lisi_graph_= lisi_graph_,
                              trajectory_=trajectory_
                             )
    results.rename(columns={results.columns[0]:setup}, inplace=True)
    if verbose:
        print(results)
    # save metrics' results
    results.to_csv(args.output)

    # Save UMAPs
    outdir = os.path.dirname(args.output)
    print(f'Calculating UMAP...')
    sc.tl.umap(adata_int)
    print(f'Saving UMAP for labels "{label_key}"...')
    fig = sc.pl.umap(adata_int, color=label_key, return_fig=True)
    fig.set_size_inches(10, 8)
    fig.savefig(os.path.join(outdir, f'{args.method}_{args.type}_labels.png'),
                bbox_inches='tight')
    print(f'Saving UMAP for batches "{batch_key}"...')
    fig = sc.pl.umap(adata_int, color=batch_key, return_fig=True)
    fig.set_size_inches(10, 8)
    fig.savefig(os.path.join(outdir, f'{args.method}_{args.type}_batch.png'),
                bbox_inches='tight')

    print("done")
