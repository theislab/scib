from . import integration, metrics, preprocessing
from . import utils as utils
from ._tools import wrap_func_naming
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
}

for alias, func in alias_func_map.items():
    if callable(func):
        func = wrap_func_naming(func, alias)
    setattr(integration, alias, func)

pp = preprocessing
ig = integration
me = metrics
cl = clustering
