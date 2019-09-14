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
    hvg = args.hvgs is not None
    cc = False
    
    base = os.path.basename(args.integrated)
    out_prefix = f'{os.path.splitext(base)[0]}_{args.type}'
    cluster_nmi = os.path.join(args.output, f'{out_prefix}_int_nmi.txt')

    ###
    
    print("reading adata before integration")
    adata = sc.read(args.uncorrected, cache=True)
    print(adata)
    if (args.hvgs is not None):
        if (adata.n_vars < args.hvgs):
            raise ValueError("There are less genes in the uncorrected adata than specified for HVG selection")

    print("reading adata after integration")
    adata_int = sc.read(args.integrated, cache=True)
    print(adata_int)

    # metric flags
    si_embed = 'X_pca'
    neighbors = True
    pca = True
    pcr_ = True
    kBET_ = False #ready to use for embedded and full matrix
    lisi_ = False
    hvg = adata.n_vars == adata_int.n_vars
    silhouette_ = True
    
    if (args.type == "embed"):
        si_embed = "X_emb"
        if ('emb' in adata_int.uns) and (adata_int.uns['emb']): # legacy check
            adata_int.obsm["X_emb"] = adata_int.obsm["X_pca"].copy()
        hvg = False
    elif (args.type == "knn"):
        hvg = False
        neighbors = False
        pcr_ = False
        silhouette_ = False
        cc = False
        kBET_ = False #until we have determined how to convert the bbknn knn-graph to FNN format, which kBET uses
        lisi_ = False

    print("reducing integrated and uncorrected data")
    scIB.preprocessing.reduce_data(adata_int, umap=False,
                                   neighbors=neighbors, pca=pca,
                                   n_top_genes=args.hvgs if hvg else None,
                                   use_rep='X_emb' if 'X_emb' in adata_int.obsm else 'X_pca')

    adata=adata[:,adata_int.var_names].copy()
    scIB.preprocessing.reduce_data(adata, umap=False,
                                   neighbors=True, pca=True, n_top_genes=None)
    
    print("computing metrics")
    results = scIB.me.metrics(adata, adata_int, hvg=hvg, cluster_nmi=cluster_nmi,
                    batch_key=batch_key, label_key=label_key,
                    silhouette_=silhouette_, si_embed=si_embed,
                    nmi_=True, ari_=True, nmi_method='max', nmi_dir=None,
                    pcr_=pcr_, kBET_=kBET_, lisi_ = lisi_, cell_cycle_=cc, verbose=False, organism=organism
                    )
    results.rename(columns={results.columns[0]:out_prefix}, inplace=True)
    # save metrics' results
    results.to_csv(os.path.join(args.output, f'{out_prefix}_metrics.csv'))
    print(results)

    print("done")

