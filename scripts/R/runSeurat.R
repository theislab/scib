args <- commandArgs(trailingOnly=T)
source('integration.R')
print(args)
sobj = loadSeuratObject(args[[1]])
hvg<-2000
if(length(args)==4) {
	hvg<-readRDS(args[[4]])
}

out = runSeurat(sobj, args[[2]])
saveSeuratObject(out, args[[3]])


