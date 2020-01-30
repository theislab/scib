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
    #get scenario info
    scenario= pd.Series([idx.split('/')[0] for idx in res.index], dtype='category')
    
    #scale by scenario
    for scen in scenario.cat.categories:
        scen_idx = res.index[scenario == scen]
        res_tmp = res.loc[scen_idx]
        #select columns with LISI score and scale (assuming scale=True in metrics call)
        #original cLISI in [0,1] with 1 good and 0 bad
        #scale to 1 good and 0 bad
        max_cl = np.max([1,res_tmp['cLISI'].max()]) #take the max of 2 and observed score
        min_cl = np.min([0,res_tmp['cLISI'].min()]) #take the min of 1 and observed score
        res.loc[scen_idx,'cLISI'] = (res.loc[scen_idx,'cLISI']-min_cl)/(max_cl-min_cl) 
        #original iLISI in [0,max] with 0 bad and max good
        max_il = np.max([1,res_tmp['iLISI'].max()])
        min_il = np.min([0,res_tmp['iLISI'].min()])
        #scale to 1 good and 0 bad
        res.loc[scen_idx,'iLISI']=(res.loc[scen_idx,'iLISI']-min_il)/(max_il-min_il)
    
    results = res #do not transpose the file
    #write output file
    results.to_csv(args.output)
