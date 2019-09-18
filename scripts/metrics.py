#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scIB
import warnings
warnings.filterwarnings('ignore')

# types of integration output
RESULT_TYPES = [
    "full", # reconstructed expression data
    "embed", # embedded/latent space
    "knn" # only corrected neighbourhood graph as output
]

if __name__=='__main__':
    """
    read adata object, compute all metrics and output csv.
    """
    
    import argparse
    import os

    parser = argparse.ArgumentParser(description='Compute all metrics')
    
    parser.add_argument('-u', '--uncorrected', required=True)
    parser.add_argument('-i', '--integrated', required=True)
    parser.add_argument('-o', '--output', required=True, help='output directory')
    
    parser.add_argument('-b', '--batch_key', required=True, help='Key of batch')
    parser.add_argument('-l', '--label_key', required=True, help='Key of annotated labels e.g. "cell_type"')
    
    parser.add_argument('--organism', required=True)
    parser.add_argument('--type', required=True, choices=RESULT_TYPES, help='Type of result: full, embed, knn\n full: scanorama, seurat, MNN\n embed: scanorama, Harmony\n knn: BBKNN')
    parser.add_argument('--hvgs', default=None, help='Number of highly variable genes', type=int)
    parser.add_argument('-v', '--verbose', action='store_true')
    
    args = parser.parse_args()
    
    verbose = args.verbose
    type_ = args.type
    batch_key = args.batch_key
    label_key = args.label_key
    organism = args.organism
    n_hvgs = args.hvgs
    
    # set prefix for output and results column name
    base = os.path.basename(args.integrated)
    out_prefix = f'{os.path.splitext(base)[0]}_{args.type}'
    cluster_nmi = os.path.join(args.output, f'{out_prefix}_int_nmi.txt')

    if verbose:
        print('Options')
        print(f'    type:\t{type_}')
        print(f'    batch_key:\t{batch_key}')
        print(f'    label_key:\t{label_key}')
        print(f'    organism:\t{organism}')
        print(f'    n_hvgs:\t{n_hvgs}')
        print(f'    out_prefix:\t{out_prefix}')
        print(f'    optimised clustering results:\t{cluster_nmi}')
    
    ###
    
    print("reading adata before integration")
    adata = sc.read(args.uncorrected, cache=True)
    print(adata)
    print("reading adata after integration")
    adata_int = sc.read(args.integrated, cache=True)
    print(adata_int)
    if (n_hvgs is not None):
        if (adata_int.n_vars < n_hvgs):
            raise ValueError("There are less genes in the uncorrected adata than specified for HVG selection")

    # DATA REDUCTION
    # select options according to type
    if adata.n_vars > adata_int.n_vars: # no HVG selection if output is not already subsetted
        n_hvgs = None
    pca = True
    neighbors = True
    embed = 'X_pca'
    
    if (type_ == "embed"):
        n_hvgs = None
        embed = "X_emb"
        # legacy check
        if ('emb' in adata_int.uns) and (adata_int.uns['emb']):
            adata_int.obsm["X_emb"] = adata_int.obsm["X_pca"].copy()
    elif (type_ == "knn"):
        n_hvgs = None
        pca = False
    
    if verbose:
        print('reduce integrated data:')
        print(f'    HVG selection:\t{n_hvgs}')
        message = f'    neighbourhood graph:\t{neighbors}'
        if neighbors:
            message += f' on {embed}'
        print(message)
        print(f'    PCA:\t{pca}')
    scIB.preprocessing.reduce_data(adata_int,
                                   n_top_genes=n_hvgs,
                                   neighbors=neighbors, use_rep=embed,
                                   pca=pca, umap=False)
    
    # METRICS
    print("computing metrics")
    silhouette_ = True
    nmi_ = True
    ari_ = True
    pcr_ = True
    cell_cycle_ = True
    kBET_ = True
    lisi_ = True
    if (type_ == "knn"):
        silhouette_ = False
        pcr_ = False
        cell_cycle_ = False
    
    if verbose:
        print(f'type:\t{type_}')
        print(f'    ASW:\t{silhouette_}')
        print(f'    NMI:\t{nmi_}')
        print(f'    ARI:\t{ari_}')
        print(f'    cell cycle:\t{cell_cycle_}')
        print(f'    kBET:\t{kBET_}')
        print(f'    LISI:\t{lisi_}')
        
    results = scIB.me.metrics(adata, adata_int, verbose=verbose,
                              hvgs=n_hvgs, cluster_nmi=cluster_nmi,
                              batch_key=batch_key, label_key=label_key,
                              silhouette_=silhouette_, embed=embed,
                              type_ = type_, 
                              nmi_=nmi_, nmi_method='arithmetic', nmi_dir=None,
                              ari_=ari_,
                              pcr_=pcr_,
                              cell_cycle_=cell_cycle_, organism=organism,
                              kBET_=kBET_,
                              lisi_=lisi_
                             )
    results.rename(columns={results.columns[0]:out_prefix}, inplace=True)
    if verbose:
        print(results)
    # save metrics' results
    results.to_csv(os.path.join(args.output, f'{out_prefix}_metrics.csv'))

    print("done")

