#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scIB
import warnings
warnings.filterwarnings('ignore')


def runPP(inPath, outPath, hvg, batch, rout, scale, seurat):
    """
    params:
        method: name of method
        batch: name of `adata.obs` column of the batch
        max_genes_hvg: maximum number of HVG
    """

    adata = sc.read(inPath)
    hvgs=adata.var.index

    # remove HVG if already precomputed
    if 'highly_variable' in adata.var:
        del adata.var['highly_variable']
    
    if hvg > 500:
        if seurat:
            hvgs= scIB.preprocessing.hvg_batch(adata,batch_key=batch, target_genes=hvg, adataOut=False)
        else:
            adata = scIB.preprocessing.hvg_batch(adata,
                                                batch_key=batch,
                                                target_genes=hvg,
                                                adataOut=True)
    if scale:
        adata = scIB.preprocessing.scale_batch(adata, batch)


    if rout:
        import rpy2.robjects as ro
        import anndata2ri
        from scipy.sparse import issparse

        ro.r('library(Seurat)')
        ro.r('library(scater)')
        anndata2ri.activate()

        if issparse(adata.X):
            if not adata.X.has_sorted_indices:
                adata.X.sort_indices()

        for key in adata.layers:
            if issparse(adata.layers[key]):
                if not adata.layers[key].has_sorted_indices:
                    adata.layers[key].sort_indices()
        ro.globalenv['adata']=adata
        ro.r('adata <- as.Seurat(adata,counts=NULL, data="X")')
        ro.r(f'saveRDS(adata, file="{outPath}")')
        
        if hvgs is not None:
            hvg_out = outPath+'_hvg.rds'
            ro.globalenv['hvgs']=hvgs
            ro.r('unlist(hvgs)')
            ro.r(f'saveRDS(hvgs, file="{hvg_out}")')


    else:
        sc.write(outPath, adata)

if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run the integration methods')

    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)
    parser.add_argument('-b', '--batch', required=True, help='Batch variable')
    parser.add_argument('-v', '--hvgs', help='Number of highly variable genes', default=2000)
    parser.add_argument('-r', '--rout', help='Save output for R methods', action='store_true')
    parser.add_argument('-s', '--scale', action='store_true', help='Scale the data per batch')
    parser.add_argument('-l', '--seurat', help='Generate output for seurat including hvg list', action='store_true')

    args = parser.parse_args()
    file = args.input_file
    out = args.output_file
    batch = args.batch
    hvg = int(args.hvgs)
    rout = args.rout
    seurat = args.seurat
    scale = args.scale
    
    runPP(file, out, hvg, batch, rout, scale, seurat)
