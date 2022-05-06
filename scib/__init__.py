try:
    from importlib import metadata
except ImportError:  # for Python<3.8
    import importlib_metadata as metadata

__version__ = metadata.version("scib")

from . import integration, metrics, preprocessing, utils
from ._package_tools import rename_func
from .metrics import clustering

alias_func_map = {
    "runScanorama": integration.scanorama,
    "runTrVae": integration.trvae,
    "runTrVaep": integration.trvaep,
    "runScGen": integration.scgen,
    "runScvi": integration.scvi,
    "runScanvi": integration.scanvi,
    "runMNN": integration.mnn,
    "runBBKNN": integration.bbknn,
    "runSaucie": integration.saucie,
    "runCombat": integration.combat,
    "runDESC": integration.desc,
    "readConos": preprocessing.read_conos,
    "readSeurat": preprocessing.read_seurat,
    "saveSeurat": preprocessing.save_seurat,
}

for alias, func in alias_func_map.items():
    rename_func(func, alias)

pp = preprocessing
ig = integration
me = metrics
cl = clustering
