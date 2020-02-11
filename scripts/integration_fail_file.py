
from snakemake.io import load_configfile
from pathlib import Path

if __name__=='__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Create an empty output file for failed integration runs')

    parser.add_argument('-c', '--config', help='Snakemake config file', required=True)
    parser.add_argument('-t', '--task', required=True)
    parser.add_argument('-m', '--method', required=True)
    parser.add_argument("-v", '--hvgs', help='pre-processed by HVG filtering', action='store_true')
    parser.add_argument('-s', '--scale', action='store_true', help='pre-processed by scaling')

    args = parser.parse_args()
    config = args.config
    task = args.task
    hvgs = args.hvgs
    scale = args.scale
    method = args.method

    # Load config file
    params = load_configfile(config)

    # Check inputs
    if method not in params['METHODS']:
        raise ValueError(f'{method} is not a valid method.\n' 
                         f'Please choose one of: {list(params["METHODS"].keys())}')

    if task not in params['DATA_SCENARIOS']:
        raise ValueError(f'{task} is not a valid integration task.\n' 
                         f'Please choose one of: {list(params["DATA_SCENARIOS"].keys())}')
        
    # Get path values
    folder = params['ROOT']
    t_folder = task
    s_folder = 'scaled' if scale else 'unscaled'
    h_folder = 'hvg' if hvgs else 'full_feature'
    r_folder = 'R/' if 'R' in params['METHODS'][method] else ''
    filename = method+'.h5ad'

    folder_path = '/'.join([folder,task,'integration',s_folder,h_folder])+'/'+r_folder
    full_path = folder_path+filename

    #print(full_path)
    Path(full_path).touch()
