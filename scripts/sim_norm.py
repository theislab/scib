#!/usr/bin/env python


import scanpy as sc
import scIB

def normSim(adata):
    adata.layers['counts'] = adata.X
    scIB.pp.normalize(adata)



if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run the integration methods')

    parser.add_argument('-i', '--input_file', required=True)
    parser.add_argument('-o', '--output_file', required=True)

    args = parser.parse_args()
    infile = args.input_file
    out = args.output_file

    adata = sc.read(infile)
    normSim(adata)
    adata.write(out)
