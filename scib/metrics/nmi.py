import os
import subprocess

from sklearn.metrics.cluster import normalized_mutual_info_score

from ..utils import check_adata, check_batch


def nmi(adata, group1, group2, method="arithmetic", nmi_dir=None):
    """Normalized mutual information

    Wrapper for normalized mutual information NMI between two different cluster assignments

    :param adata: Anndata object
    :param group1: column name of ``adata.obs``
    :param group2: column name of ``adata.obs``
    :param method: NMI implementation.
        'max': scikit method with ``average_method='max'``;
        'min': scikit method with ``average_method='min'``;
        'geometric': scikit method with ``average_method='geometric'``;
        'arithmetic': scikit method with ``average_method='arithmetic'``;
        'Lancichinetti': implementation by A. Lancichinetti 2009 et al. https://sites.google.com/site/andrealancichinetti/mutual;
        'ONMI': implementation by Aaron F. McDaid et al. https://github.com/aaronmcdaid/Overlapping-NMI
    :param nmi_dir: directory of compiled C code if 'Lancichinetti' or 'ONMI' are specified as ``method``.
        These packages need to be compiled as specified in the corresponding READMEs.
    :return: Normalized mutual information NMI value
    """

    check_adata(adata)
    check_batch(group1, adata.obs)
    check_batch(group2, adata.obs)

    group1 = adata.obs[group1].tolist()
    group2 = adata.obs[group2].tolist()

    if len(group1) != len(group2):
        raise ValueError(
            f"different lengths in group1 ({len(group1)}) and group2 ({len(group2)})"
        )

    # choose method
    if method in ["max", "min", "geometric", "arithmetic"]:
        nmi_value = normalized_mutual_info_score(group1, group2, average_method=method)
    elif method == "Lancichinetti":
        nmi_value = nmi_Lanc(group1, group2, nmi_dir=nmi_dir)
    elif method == "ONMI":
        nmi_value = onmi(group1, group2, nmi_dir=nmi_dir)
    else:
        raise ValueError(f"Method {method} not valid")

    return nmi_value


def onmi(group1, group2, nmi_dir=None, verbose=True):
    """
    Based on implementation https://github.com/aaronmcdaid/Overlapping-NMI
    publication: Aaron F. McDaid, Derek Greene, Neil Hurley 2011

    :param group1: list or series of cell assignments
    :param group2: list or series of cell assignments
    :param nmi_dir: directory of compiled C code
    """

    if nmi_dir is None:
        raise FileNotFoundError(
            "Please provide the directory of the compiled C code from "
            "https://sites.google.com/site/andrealancichinetti/mutual3.tar.gz"
        )

    group1_file = write_tmp_labels(group1, to_int=False)
    group2_file = write_tmp_labels(group2, to_int=False)

    nmi_call = subprocess.Popen(
        [nmi_dir + "onmi", group1_file, group2_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    stdout, stderr = nmi_call.communicate()
    if stderr:
        print(stderr)

    nmi_out = stdout.decode()
    if verbose:
        print(nmi_out)

    nmi_split = [x.strip().split("\t") for x in nmi_out.split("\n")]
    nmi_max = float(nmi_split[0][1])

    # remove temporary files
    os.remove(group1_file)
    os.remove(group2_file)

    return nmi_max


def nmi_Lanc(group1, group2, nmi_dir="external/mutual3/", verbose=True):
    """
    paper by A. Lancichinetti 2009
    https://sites.google.com/site/andrealancichinetti/mutual
    recommended by Malte
    """

    if nmi_dir is None:
        raise FileNotFoundError(
            "Please provide the directory of the compiled C code from https://sites.google.com/site/andrealancichinetti/mutual3.tar.gz"
        )

    group1_file = write_tmp_labels(group1, to_int=False)
    group2_file = write_tmp_labels(group2, to_int=False)

    nmi_call = subprocess.Popen(
        [nmi_dir + "mutual", group1_file, group2_file],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    stdout, stderr = nmi_call.communicate()
    if stderr:
        print(stderr)
    nmi_out = stdout.decode().strip()

    return float(nmi_out.split("\t")[1])


def write_tmp_labels(group_assignments, to_int=False, delim="\n"):
    """
    write the values of a specific obs column into a temporary file in text format
    needed for external C NMI implementations (onmi and nmi_Lanc functions), because they require files as input

    :param group_assignments: list or series of cell assignments
    :param to_int: rename the unique column entries by integers in range(1,len(group_assignments)+1)
    """
    import tempfile

    if to_int:
        label_map = {}
        i = 1
        for label in set(group_assignments):
            label_map[label] = i
            i += 1
        delim.join([str(label_map[name]) for name in group_assignments])
    else:
        delim.join([str(name) for name in group_assignments])

    clusters = {label: [] for label in set(group_assignments)}
    for i, label in enumerate(group_assignments):
        clusters[label].append(str(i))

    output = "\n".join([" ".join(c) for c in clusters.values()])
    output = str.encode(output)

    # write to file
    with tempfile.NamedTemporaryFile(delete=False) as f:
        f.write(output)
        filename = f.name

    return filename
