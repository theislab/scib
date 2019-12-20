args <- commandArgs(trailingOnly=T)
source('integration.R')
print(args)
sobj = loadSeuratObject(args[[1]])
out = func_profiler(runSeurat(sobj, args[[2]]))
saveSeuratObject(out$results, args[[3]])

