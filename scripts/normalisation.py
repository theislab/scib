#!/usr/bin/env python
# coding: utf-8

import argparse

parser = argparse.ArgumentParser(description='Normalisation with scran')
parser.add_argument('-i', '--input_file', required=True)
parser.add_argument('-o', '--output_file', required=True)
args = parser.parse_args()

input_file = args.input_file
output_file = args.output_file

# # Normalisation

import pandas as pd
import scanpy as sc
import scIB

adata = sc.read(input_file, cache=True)
scIB.preprocessing.normalize(adata)
adata.write(output_file)
