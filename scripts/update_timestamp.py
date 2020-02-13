#!/usr/bin/env python

from snakemake.io import load_configfile
from pathlib import Path
from os.path import isfile

def touch_if_exists(file_path):
    if isfile(file_path):
        Path(file_path).touch()
    else:
        print(f'{file_path} does not exist.')
        

def update_timestamp_task(config, task, update_metrics=True):
    """
    This function updates the timestamps of all the files generated for a scenario.
    Updating timestamps is done in order of file generation by the snakemake pipeline
    such that snakemake doesn't unnecessarily rerun rules that generate files that
    already exist..

    Note that this function does not update the timestamp of the aggregated metrics 
    files.
    """
    
    base_folder = config['ROOT']
    scaling = config['SCALING']
    hvgs = list(config['FEATURE_SELECTION'].keys())
    methods = list(config['METHODS'].keys())

    # Prep files
    prep_endings = ['.h5ad', '.RDS', '_hvg.RDS']

    for scal in scaling:
        for feat in hvgs:
            folder_path = '/'.join([folder,task,'prepare',scal,feat])+'/'
            file_base = 'adata_pre'
            for end in prep_endings:
                filename = file_base+end
                full_path = folder_path+filename
                touch_if_exists(file_path)

            full_path = folder_path+'prep_h5ad.benchmark'
            touch_if_exists(file_path)

            full_path = folder_path+'prep_RDS.benchmark'
            touch_if_exists(file_path)

    # Integration & convert files
    for scal in scaling:
        for feat in hvgs:
            folder_base = '/'.join([folder,task,'integration',scal,feat])+'/'
            for method in methods:
                if 'R' in config['METHODS']:
                    r_folder = 'R/'
                    method_endings = ['.RDS', '.RDS.benchmark', '.h5ad']
                else:
                    r_folder = ''
                    method_endings = ['.h5ad', '.h5ad.benchmark']

                for end in method_endings:
                    folder_path = folder_base+r_folder
                    full_path = folder_path+method+end
                    touch_if_exists(file_path)

    # Metrics files
    metric_endings = ['_int_nmi.txt', '.csv']

    if update_metrics:
        for scal in scaling:
            for feat in hvgs:
                folder_base = '/'.join([folder,task,'metrics',scal,feat])+'/'
                for method in methods:
                    for out_type in config['METHODS'][method]['output_type']:
                        file_base = '_'.join([method,out_type])
                        for end in metric_endings:
                            file_path = folder_base+file_base+end
                            touch_if_exists(file_path)
        


if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Update timestamp on all output files'
                                     ' for an integration task')

    parser.add_argument('-c', '--config', help='Snakemake config file', required=True)
    parser.add_argument('-t', '--task', help='Integration task to update', 
                        required=True)
    parser.add_argument('-m', '--include-metrics', action='store_true', 
                        help='Also update timestamp of metrics files')

    args = parser.parse_args()
    config = args.config
    task = args.task
    metrics_flag = args.include_metrics

    # Load config file
    params = load_configfile(config)

    if task not in params['DATA_SCENARIOS']:
        raise ValueError(f'{task} is not a valid integration task.\n' 
                         'Please choose one of:\n'
                         f'{list(params["DATA_SCENARIOS"].keys())}')

    update_timestamp_task(params, task, metrics_flag)
