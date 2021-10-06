from . import integration, metrics, preprocessing
from . import utils as utils
from ._tools import rename_func
from .metrics import clustering

alias_func_map = {
    'runScanorama': integration.scanorama,
    'runTrVae': integration.trvae,
    'runTrVaep': integration.trvaep,
    'runScGen': integration.scgen,
    'runScvi': integration.scvi,
    'runScanvi': integration.scanvi,
    'runMNN': integration.mnnpy,
    'runBBKNN': integration.bbknn,
    'runSaucie': integration.saucie,
    'runCombat': integration.combat,
    'runDESC': integration.desc,
    'readSeurat': preprocessing.read_seurat,
    'readConos': preprocessing.read_conos,
}

for alias, func in alias_func_map.items():
    rename_func(func, alias)

pp = preprocessing
ig = integration
me = metrics
cl = clustering
