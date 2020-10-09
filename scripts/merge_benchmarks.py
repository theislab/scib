import pandas as pd
import argparse
import os

if __name__=='__main__':
    """
    Merge benchmark output for all scenarios, methods and settings
    """

    parser = argparse.ArgumentParser(description='Collect all benchmarks')

    parser.add_argument('-o', '--output', required=True, help='output file')
    parser.add_argument('-r', '--root', required=True,
                        help='root directory for scIB output')
    args = parser.parse_args()


    print("Searching for .benchmark files...")
    bench_files = []
    for path, dirs, files in os.walk(args.root):
        for file in files:
            if 'integration' in path and file.endswith('.benchmark'):
                bench_files.append(os.path.join(path, file))
    print(f"Found {len(bench_files)} .benchmark files")

    print("Collecting benchmarks...")
    res_list = []
    for file in bench_files:

        if os.stat(file).st_size == 0:
            print(f'{file} is empty and will be skipped')
            continue

        clean_name = file.replace(args.root, "").replace(".benchmark", "")
        res = pd.read_csv(file, sep='\t')
        res.rename(columns={res.columns[1]: 'h_m_s'}, inplace=True)
        res['scenario'] = clean_name
        res.set_index('scenario', inplace=True)
        res_list.append(res)

    print(f"Saving merged benchmarks to {args.output}")
    results = pd.concat(res_list)
    results.to_csv(args.output, index_label='scenario')

    print("Done!")

