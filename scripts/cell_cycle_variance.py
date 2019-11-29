import pandas as pd
import scanpy as sc
import scIB

RESULT_TYPES = ["full", "embed"]

if __name__=='__main__':
    
    import argparse
    import os

    parser = argparse.ArgumentParser(description='Compute all PC regression scores of cell cycle metric')
    
    parser.add_argument('-u', '--uncorrected', required=True)
    parser.add_argument('-i', '--integrated', required=True)
    parser.add_argument('-o', '--output', required=True, help='output file')
    parser.add_argument('-b', '--batch_key', required=True, help='Key of batch')
    parser.add_argument('--organism', required=True)
    parser.add_argument('--type', required=True, choices=RESULT_TYPES, help='Type of result: full, embed, knn\n full: scanorama, seurat, MNN\n embed: scanorama, Harmony')
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()
    
    embed = 'X_emb' if args.type == 'embed' else 'X_pca'
    
    print("reading adata before integration")
    adata_pre = sc.read(args.uncorrected, cache=True)
    print("reading adata after integration")
    adata_post = sc.read(args.integrated, cache=True)
    
    scores = scIB.me.cell_cycle(adata_pre, adata_post,
                                batch_key=args.batch_key,
                                organism=args.organism,
                                embed=embed,
                                agg_func=None, 
                                verbose=args.verbose)
    print(scores)
    scores.to_csv(args.output, index=False)
