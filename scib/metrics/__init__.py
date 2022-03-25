"""
This package contains all the metrics used for benchmarking scRNA-seq data integration performance.
The metrics can be classified into biological conservation and batch removal metrics:

Biological conservation
    + HVG overlap
    + Cell type ASW
    + Isolated label ASW
    + Isolated label F1
    + NMI cluster/label
    + ARI cluster/label
    + Cell cycle conservation
    + cLISI

Batch correction
    + Graph connectivity
    + Batch ASW
    + PC regression
    + kBET
    + iLISI

"""

from .metrics import *
from .lisi import lisi_graph
