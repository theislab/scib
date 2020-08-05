import pathlib
from snakemake.io import expand, load_configfile


def as_list(x):
        return x if isinstance(x, list) else [x]


def join_path(*args):
    path = pathlib.Path(args[0])
    for d in args[1:]:
        path = path / d
    return str(path)


class ParsedConfig:

    OUTPUT_FILE_TYPES = ['prepare', 'integration', 'metrics', 'metrics_unintegrated', 'cc_variance']
    OUTPUT_LEVELS     = ['single', 'final', 'scaled_final', 'by_method', 'by_method_scaling', 
                         'directory_by_setting']
    OUTPUT_TYPES      = ['full', 'embed', 'knn']
        
    def __init__(self, config):

        # TODO: define and check schema of config

        self.ROOT              = config["ROOT"]
        self.DATA_SCENARIOS    = config["DATA_SCENARIOS"]
        self.SCALING           = config["SCALING"]
        self.FEATURE_SELECTION = config["FEATURE_SELECTION"]
        self.METHODS           = config["METHODS"]
        self.timing            = config["timing"]
        self.r_env             = config["r_env"]
        self.py_env            = config["py_env"]
        self.conv_env          = config["conv_env"]
        try:
            self.unintegrated_m = config["unintegrated_metrics"]
        except:
            self.unintegrated_m = False


    def get_all_scalings(self):
        return self.SCALING


    def get_all_feature_selections(self):
        return self.FEATURE_SELECTION.keys()

    # --------------------------------------------------------------------------
    # Gets all available methods. filter for framework (R/python) if needed.
    #
    # @param framework Only methods based on the framework will be retrieved.
    # one of ["python", "R", "both"], default: both
    # --------------------------------------------------------------------------
    def get_all_methods(self, framework="both"):
        all_methods = []
        for method in self.METHODS:
            is_r = self.get_from_method(method, "R")
            if framework == "both":
                all_methods.append(method)
            elif (framework == "python") and (not is_r):
                all_methods.append(method)
            elif (framework == "R") and (is_r):
                all_methods.append(method)

        return all_methods


    def get_all_scenarios(self):
        return self.DATA_SCENARIOS.keys()


    def get_feature_selection(self, key):

        if key not in self.FEATURE_SELECTION:
            raise ValueError(f"{key} not a valid key for scaling")

        return self.FEATURE_SELECTION[key]

    def get_from_method(self, method, key):

        if method not in self.METHODS:
            raise ValueError(f"{method} not defined as method")

        if key not in self.METHODS[method]:
            return False
            #raise ValueError(f"{key} not a valid attribute of scenario {scenario}")

        return self.METHODS[method][key]

    def get_hvg(self, wildcards, adata_path):

        n_hvgs = self.get_feature_selection(wildcards.hvg)

        if n_hvgs == 0:
                return ""

        if self.get_from_method(wildcards.method, "R"):
            path_parts = adata_path.split('.')
            path_parts[-2] += '_hvg'
            hvg_path = '.'.join(path_parts)
            return f'-v "{hvg_path}"'

        return f"-v {n_hvgs}"


    def get_from_scenario(self, scenario, key):

        if scenario not in self.DATA_SCENARIOS:
            raise ValueError(f"{scenario} not defined as scenario")

        if key not in self.DATA_SCENARIOS[scenario]:
            return False

        return self.DATA_SCENARIOS[scenario][key]


    def get_filename_pattern(self, file_type, level, file_suffix=None):
        """
        file_type: ParsedConfig.OUTPUT_FILE_TYPES
        level: one of ParsedConfig.OUTPUT_LEVELS
        """

        if file_type not in ParsedConfig.OUTPUT_FILE_TYPES:
            raise ValueError(f"{file_type} not a valid output file type")

        if level not in ParsedConfig.OUTPUT_LEVELS:
            raise ValueError(f"{level} not a valid output level")

        if file_suffix not in ["rds", "rds_to_h5ad", "h5ad", None]:
            raise ValueError(f"{file_suffix} not a valid output file suffix")

        file_suffixes = {
            "prepare"     : "{method}.h5ad",
            "integration" : "{method}.h5ad",
            "metrics"     : "{method}_{o_type}.csv",
            "cc_variance" : "{method}_{o_type}.csv"
        }

        # in case of R, we need a different suffix for the integration part
        if file_suffix == "rds":
            file_suffixes["integration"] = "R/{method}.RDS"
        elif file_suffix == "rds_to_h5ad":
            file_suffixes["integration"] = "R/{method}.h5ad"

        suffix = file_suffixes[file_type]
        
        if level == "single":
            return join_path(self.ROOT, "{scenario}", file_type, "{scaling}", "{hvg}", suffix)
        elif level == "directory_by_setting":
            return join_path(self.ROOT, "{scenario}", file_type, "{scaling}", "{hvg}")
        elif level == "by_method":
            return join_path(self.ROOT, "{{scenario}}", file_type, "{{scaling}}", "{{hvg}}", suffix)
        elif level == "by_method_scaling":
            return join_path(self.ROOT, "{{scenario}}", file_type, "{scaling}", "{{hvg}}", suffix)
        elif level == "final":
            return join_path(self.ROOT, f"{file_type}.csv")
        elif level == "scaled_final":
            return join_path(self.ROOT, f"{file_type}_scaled.csv")

        


    def get_all_file_patterns(self, file_type, output_types=None):
        """
        Collect all expanded file patterns
        file_type: one of ['integration', 'metrics', 'cc_variance']
        output_types: output type or list of output types to be considered.
            If output_types==None, all output types are included.
            Useful if a certain metric is examined on a specific output type.
            Output types are ['full', 'embed', 'knn']
        """

        if file_type not in ParsedConfig.OUTPUT_FILE_TYPES:
            raise ValueError(f"{file_type} not a valid output file type")

        if output_types is None:
            output_types = ParsedConfig.OUTPUT_TYPES
        else:
            for ot in output_types:
                if ot not in ParsedConfig.OUTPUT_TYPES:
                    raise ValueError(f"{output_types} not a valid output type")

        all_files = []
        
        if file_type == 'metrics_unintegrated':
                # add unintegrated
                file_pattern = self.get_filename_pattern("metrics", "single")
                all_files = expand(file_pattern,
                                   scenario=self.get_all_scenarios(),
                                   hvg="full_feature",
                                   scaling="unscaled",
                                   method="unintegrated", o_type="full")
        else:
            for method in self.METHODS:

                ot = self.get_from_method(method, "output_type")
                # keep only types to be considered
                ot = [x for x in as_list(ot) if x in as_list(output_types)]
                if not ot:
                    continue # skip if method does not have any

                # skip scaling for certain methods
                scaling = self.SCALING.copy()
                skip_scaling = self.get_from_method(method, "no_scale")
                if skip_scaling and "scaled" in scaling:
                    scaling.remove("scaled")

                file_level = "by_method_scaling"

                if file_type == 'integration':
                    is_r = self.get_from_method(method, "R")
                    if is_r:
                        file_pattern = self.get_filename_pattern(file_type, file_level, "rds_to_h5ad")
                    else:
                        file_pattern = self.get_filename_pattern(file_type, file_level, "h5ad")

                    expanded = expand(file_pattern, method=method, scaling=scaling)

                elif file_type in ['metrics', 'cc_variance']:
                    file_pattern = self.get_filename_pattern(file_type, file_level)
                    expanded = expand(file_pattern, method=method, o_type=ot, scaling=scaling)
                else:
                    raise ValueError(f"{file_type} is not a valid file type")

                for p in expanded:
                    f = expand(p, scenario=self.get_all_scenarios(),
                               hvg=self.get_all_feature_selections())
                    all_files.extend(f)

        return all_files

