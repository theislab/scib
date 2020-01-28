#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import scIB
import warnings
warnings.filterwarnings('ignore')


def runPost(inPath, outPath, conos):
    """
    params:
        inPath: path of the R object
        outPath: path of the processed anndata object
        conos: set if input is conos obect
    """
    if conos:
        adata = scIB.pp.readConos(inPath)
    else:
        adata = scIB.pp.readSeurat(inPath)

    adata.write(outPath)


if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run the integration methods')

    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)
    parser.add_argument('-c', '--conos', help='set for conos input', action='store_true')

    args = parser.parse_args()
    file = args.input_file
    out = args.output_file
    conos = args.conos
    
    runPost(file, out, conos)
