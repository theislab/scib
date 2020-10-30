#!/usr/bin/env python

from snakemake.io import load_configfile
from pathlib import Path
from os.path import isfile, dirname

def rewrite(file_path, file_path_recomp):
    if isfile(file_path):
        with open(file_path, 'r') as file:
            data = file.readlines()
        try:    
            with open(file_path_recomp, 'r') as file:
                data_recomp = file.readlines()
                data[1] = data_recomp[1]
                data[2] = data_recomp[2]
                data[7] = data_recomp[7]
        except:
            print('recomp file not present')
        new_path = file_path.replace("metrics", "metrics_combined", 1)
        print(dirname(new_path))
        Path(dirname(new_path)).mkdir(parents=True, exist_ok=True)
        with open(new_path, 'w') as file:
            file.writelines(data)
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
    # Metrics files
    metric_endings = ['.csv']

    if True:
        for scal in scaling:
            for feat in hvgs:
                folder_base = '/'.join([base_folder,task,'metrics',scal,feat])+'/'
                folder_base_recomp = '/'.join([base_folder,task,'metrics_recomp',scal,feat])+'/'
                for method in methods:
                    out_types = config['METHODS'][method]['output_type']
                    if isinstance(out_types, str):
                        out_types = [out_types]

                    for out_type in out_types:
                        file_base = '_'.join([method,out_type])
                        for end in metric_endings:
                            full_path = folder_base+file_base+end
                            full_path_recomp = folder_base_recomp+file_base+end
                            rewrite(full_path, full_path_recomp)
        


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
