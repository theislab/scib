"""
This package contains all the metrics used for benchmarking scRNA-seq data integration performance.
The metrics can be classified into biological conservation and batch removal metrics.

+----------------------------------------------------+------------------------------------------------------+
| Biological Conservation                            | Batch Correction                                     |
+====================================================+======================================================+
| + Cell type ASW                                    | + Batch ASW                                          |
|                                                    |                                                      |
| + Cell cycle conservation                          | + Principal component regression                     |
|                                                    |                                                      |
| + cLISI                                            | + iLISI                                              |
|                                                    |                                                      |
| + ARI cluster/label                                | + Graph connectivity                                 |
|                                                    |                                                      |
| + NMI cluster/label                                | + kBET                                               |
|                                                    |                                                      |
| + Highly variable gene overlap                     |                                                      |
|                                                    |                                                      |
| + Isolated label ASW                               |                                                      |
|                                                    |                                                      |
| + Isolated label F1                                |                                                      |
|                                                    |                                                      |
| + Trajectory conservation                          |                                                      |
+----------------------------------------------------+------------------------------------------------------+


"""

from .metrics import *
from .lisi import lisi_graph
from .pcr import pc_regression, pcr
