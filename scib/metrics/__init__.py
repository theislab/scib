# fmt: off
try:
    from sklearnex import patch_sklearn
    patch_sklearn()
except ImportError:
    print('INFO: sklearnex not installed.  No hardware acceleration for sklearn-based metrics available.')

from .ari import ari
from .cell_cycle import cell_cycle
from .clustering import (cluster_optimal_resolution, get_resolutions,
                         opt_louvain)
from .graph_connectivity import graph_connectivity
from .highly_variable_genes import hvg_overlap
from .isolated_labels import (isolated_labels, isolated_labels_asw,
                              isolated_labels_f1)
from .kbet import kBET
from .lisi import clisi_graph, ilisi_graph, lisi_graph
from .metrics import metrics, metrics_all, metrics_fast, metrics_slim
from .nmi import nmi
from .pcr import pc_regression, pcr, pcr_comparison
from .silhouette import silhouette, silhouette_batch
from .trajectory import trajectory_conservation
