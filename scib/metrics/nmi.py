import os
import subprocess

from sklearn.metrics.cluster import normalized_mutual_info_score

try:
    from scanpy._utils import renamed_arg
except ImportError:
    from .._package_tools import renamed_arg

from ..utils import check_adata, check_batch


@renamed_arg("group1", "cluster_key")
@renamed_arg("group2", "label_key")
@renamed_arg("method", "implementation")
def nmi(adata, cluster_key, label_key, implementation="arithmetic", nmi_dir=None):
    """Normalized mutual information

    The normalized mutual information is a version of the mutual information corrected by the entropy of clustering and
    ground truth labels (e.g. cell type).
    The score ranges between 0 and 1, with 0 representing no sharing and 1 representing perfect sharing of information
    between clustering and annotated cell labels.

    :param adata: anndata object with cluster assignments in ``adata.obs[cluster_key]``
    :param cluster_key: string of column in adata.obs containing cluster assignments
    :param label_key: string of column in adata.obs containing labels
    :param implementation: NMI implementation.
        ``'max'``: scikit method with ``average_method='max'``;
        ``'min'``: scikit method with ``average_method='min'``;
        ``'geometric'``: scikit method with ``average_method='geometric'``;
        ``'arithmetic'``: scikit method with ``average_method='arithmetic'``;
        ``'Lancichinetti'``: implementation by A. Lancichinetti 2009 et al. https://sites.google.com/site/andrealancichinetti/mutual;
        ``'ONMI'``: implementation by Aaron F. McDaid et al. https://github.com/aaronmcdaid/Overlapping-NMI
    :param nmi_dir: directory of compiled C code if 'Lancichinetti' or 'ONMI' are specified as ``method``.
        These packages need to be compiled as specified in the corresponding READMEs.
    :return: Normalized mutual information NMI value

    This function can be applied to all integration output types.
    The ``adata`` must contain cluster assignments that are based off the knn graph given or derived from the integration
    method output.
    For this metric you need to include all steps that are needed for clustering.
    See :ref:`preprocessing` for more information on preprocessing.

    **Examples**

    .. code-block:: python

        # feature output
        scib.pp.reduce_data(
            adata, n_top_genes=2000, batch_key="batch", pca=True, neighbors=True
        )
        scib.me.cluster_optimal_resolution(adata, cluster_key="cluster", label_key="celltype")
        scib.me.nmi(adata, cluster_key="cluster", label_key="celltype")

        # embedding output
        sc.pp.neighbors(adata, use_rep="X_emb")
        scib.me.cluster_optimal_resolution(adata, cluster_key="cluster", label_key="celltype")
        scib.me.nmi(adata, cluster_key="cluster", label_key="celltype")

        # knn output
        scib.me.cluster_optimal_resolution(adata, cluster_key="cluster", label_key="celltype")
        scib.me.nmi(adata, cluster_key="cluster", label_key="celltype")

    """

    check_adata(adata)
    check_batch(cluster_key, adata.obs)
    check_batch(label_key, adata.obs)

    cluster_key = adata.obs[cluster_key].tolist()
    label_key = adata.obs[label_key].tolist()

    if len(cluster_key) != len(label_key):
        raise ValueError(
            f"different lengths in cluster_key ({len(cluster_key)}) and label_key ({len(label_key)})"
        )

    # choose method
    if implementation in ["max", "min", "geometric", "arithmetic"]:
        nmi_value = normalized_mutual_info_score(
            cluster_key, label_key, average_method=implementation
        )
    elif implementation == "Lancichinetti":
        nmi_value = nmi_Lanc(cluster_key, label_key, nmi_dir=nmi_dir)
    elif implementation == "ONMI":
        nmi_value = onmi(cluster_key, label_key, nmi_dir=nmi_dir)
    else:
        raise ValueError(f"Method {implementation} not valid")

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
