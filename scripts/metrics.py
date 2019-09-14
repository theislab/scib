#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scIB
import warnings
warnings.filterwarnings('ignore')

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
    parser.add_argument('-g', '--organism', required=True)
    parser.add_argument('-t', '--type', required=True, help='Type of result: full, embed, knn\n full: scanorama, seurat, MNN\n embed: scanorama, Harmony\n knn: BBKNN')
    parser.add_argument('-v', '--hvgs', default=None, help='Number of highly variable genes', type=int)
    args = parser.parse_args()
    
    result_types = [
        "full", # reconstructed expression data
        "embed", # embedded/latent space
        "knn" # only corrected neighbourhood graph as output
    ]
    if args.type not in result_types:
        raise ValueError(f'{args.type} is not a valid result type flag')
    
    batch_key = args.batch_key
    label_key = args.label_key
    organism = args.organism
    n_hvgs = args.hvgs
    
    base = os.path.basename(args.integrated)
    out_prefix = f'{os.path.splitext(base)[0]}_{args.type}'
    cluster_nmi = os.path.join(args.output, f'{out_prefix}_int_nmi.txt')

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

    # preprocessing
    if adata.n_vars > adata_int.n_vars: # no HVG selection if output is not already subsetted
        n_hvgs = None
    pca = True
    neighbors = True
    embed = 'X_pca'
    
    # default metric flags
    silhouette_ = True
    nmi_ = True
    ari_ = True
    pcr_ = True
    kBET_ = False
    cell_cycle_ = True
    
    if (args.type == "embed"):
        n_hvgs = None
        embed = "X_emb"
        # legacy check
        if ('emb' in adata_int.uns) and (adata_int.uns['emb']):
            adata_int.obsm["X_emb"] = adata_int.obsm["X_pca"].copy()
    elif (args.type == "knn"):
        n_hvgs = None
        pca = False
        neighbors = False
        silhouette_ = False
        pcr_ = False
        cell_cycle_ = False
    
    print("reducing integrated and uncorrected data")
    scIB.preprocessing.reduce_data(adata_int,
                                   n_top_genes=n_hvgs,
                                   neighbors=neighbors, use_rep=embed,
                                   pca=pca, umap=False)
    # select HVGs according to corrected data
    adata = adata[:,adata_int.var_names].copy()
    scIB.preprocessing.reduce_data(adata,
                                   n_top_genes=None,
                                   neighbors=True,
                                   pca=True, umap=False)
    
    print("computing metrics")
    results = scIB.me.metrics(adata, adata_int,
                              hvg=n_hvgs is not None, cluster_nmi=cluster_nmi,
                              batch_key=batch_key, label_key=label_key,
                              silhouette_=silhouette_, embed=embed,
                              nmi_=nmi_, nmi_method='arithmetic', nmi_dir=None,
                              ari_=ari_,
                              pcr_=pcr_,
                              kBET_=kBET_,
                              cell_cycle_=cell_cycle_, organism=organism,
                              verbose=False
                    )
    results.rename(columns={results.columns[0]:out_prefix}, inplace=True)
    # save metrics' results
    results.to_csv(os.path.join(args.output, f'{out_prefix}_metrics.csv'))

    print("done")

