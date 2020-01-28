#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import scIB
import warnings
warnings.filterwarnings('ignore')
import argparse
from functools import reduce

if __name__=='__main__':
    """
    Scale LISI scores in the output of merge_metrics.py
    """
    
    parser = argparse.ArgumentParser(description='Scale LISI scores')

    parser.add_argument('-i', '--input', required=True, help='Merged metrics file in csv file format')
    parser.add_argument('-o', '--output', required=True, help='Merged and LISI-scaled metrics file in csv file format')
    
    args = parser.parse_args()
    
    #get input file
    file = args.input
    res = pd.read_csv(file, index_col=0)
    
    #select columns with LISI score and scale
   
    
    results = res.T
    #write output file
    results.to_csv(args.output)
