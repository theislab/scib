#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scIB
import warnings
warnings.filterwarnings('ignore')


def runIntegration(inPath, outPath, method, hvg, batch, scale):
    """
    params:
        method: name of method
        batch: name of `adata.obs` column of the batch
        max_genes_hvg: maximum number of HVG
    """

    adata = sc.read(inPath)

    # Prevent HVG selection for Seurat
    if (method == scIB.integration.runSeurat) and (hvg > 500):
        raise ValueError('Should not run Seurat with a priori HVG selection. '
                         'Seurat will run its own feature selection.')
    
    if hvg > 500:
        adata = scIB.preprocessing.hvg_batch(adata,
                                             batch_key=batch,
                                             target_genes=hvg,
                                             adataOut=True)
    
    if scale:
        scIB.preprocessing.scale_batch(adata, batch)

    integrated_tmp = scIB.metrics.measureTM(method, adata, batch)

    integrated = integrated_tmp[2][0]


    integrated.uns['mem'] = integrated_tmp[0]
    integrated.uns['runtime'] = integrated_tmp[1]

    sc.write(outPath, integrated)

if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run the integration methods')

    parser.add_argument('-m', '--method', required=True)
    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)
    parser.add_argument('-b', '--batch', required=True, help='Batch variable')
    parser.add_argument('-v', '--hvgs', help='Number of highly variable genes', default=2000)
    parser.add_argument('-s', '--scale', action='store_true', help='Scale the data per batch')

    args = parser.parse_args()
    file = args.input_file
    out = args.output_file
    batch = args.batch
    hvg = int(args.hvgs)
    scale = args.scale
    method = args.method
    methods = {
        'scanorama': scIB.integration.runScanorama,
        'trvae': scIB.integration.runTrVae,
        'seurat': scIB.integration.runSeurat,
        'harmony': scIB.integration.runHarmony,
        'mnn': scIB.integration.runMNN,
        'bbknn': scIB.integration.runBBKNN,
        'conos': scIB.integration.runConos,
        'scvi': scIB.integration.runScvi
    }
    
    if method not in methods.keys():
        raise ValueError('Method does not exist. Please use one of the following:\n'+str(list(methods.keys())))
    
    run= methods[method]
    runIntegration(file, out, run, hvg, batch, scale)
