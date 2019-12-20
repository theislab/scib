args <- commandArgs(trailingOnly=T)
source('integration.R')
print(args)
out = func_profiler(runConos(args[[1]], args[[2]]))
saveConos(out$results, args[[3]])

