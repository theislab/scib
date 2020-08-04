#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scIB
import warnings
warnings.filterwarnings('ignore')


def runIntegration(inPath, outPath, method, hvg, batch):
    """
    params:
        method: name of method
        batch: name of `adata.obs` column of the batch
        max_genes_hvg: maximum number of HVG
    """

    adata = sc.read(inPath)
    
    if timing:
        integrated_tmp = scIB.metrics.measureTM(method, adata, batch)

        integrated = integrated_tmp[2][0]


        integrated.uns['mem'] = integrated_tmp[0]
        integrated.uns['runtime'] = integrated_tmp[1]

    else:
        integrated = method(adata, batch)
    
    sc.write(outPath, integrated)

if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run the integration methods')

    parser.add_argument('-m', '--method', required=True)
    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)
    parser.add_argument('-b', '--batch', required=True, help='Batch variable')
    parser.add_argument('-v', '--hvgs', help='Number of highly variable genes', default=2000)
    parser.add_argument("-t", '--timing', help='Activate runtime and memory profiling', action='store_true')

    args = parser.parse_args()
    file = args.input_file
    out = args.output_file
    batch = args.batch
    hvg = int(args.hvgs)
    timing = args.timing
    method = args.method
    methods = {
        'scanorama': scIB.integration.runScanorama,
        'trvae': scIB.integration.runTrVae,
        'trvaep': scIB.integration.runTrVaep,
        'mnn': scIB.integration.runMNN,
        'bbknn': scIB.integration.runBBKNN,
        'scvi': scIB.integration.runScvi,
        'combat': scIB.integration.runCombat,
        'saucie': scIB.integration.runSaucie
        'desc': scIB.integration.runDESC
    }
    
    if method not in methods.keys():
        raise ValueError('Method does not exist. Please use one of the following:\n'+str(list(methods.keys())))
    
    run= methods[method]
    runIntegration(file, out, run, hvg, batch)
