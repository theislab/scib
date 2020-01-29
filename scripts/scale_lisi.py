#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
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
    #original cLISI in [1,max] with 1 good and max bad
    #scale to 1 good and 0 bad
    max_cl = np.max([2,res['cLISI'].max()]) #take the max of 2 and observed score
    min_cl = np.min([1,res['cLISI'].min()]) #take the min of 1 and observed score
    res['cLISI'] = (res['cLISI']-min_cl)/(max_cl-min_cl) 
    
    #original iLISI in [1,max] with 1 bad and max good
    max_il = np.max([2,res['iLISI'].max()])
    min_il = np.min([1,res['iLISI'].min()])
    #scale to 1 good and 0 bad
    res['iLISI']=(res['iLISI']-min_il)/(max_il-min_il)
    
    results = res #do not transpose the file
    #write output file
    results.to_csv(args.output)
