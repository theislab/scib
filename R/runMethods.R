library('optparse')

option_list <- list(make_option(c("-m", "--method"), type="character", default=NA, help="integration method to use"),
		    make_option(c("-i", "--input"), type="character", default=NA, help="input data"),
		    make_option(c("-o", "--output"), type="character", default=NA, help="output file"),
		    make_option(c("-b", "--batch"), type="character", default=NA, help="batch variable"),
		    make_option(c("-v", "--hvg"), type="character", default=NA, help="hvg list for seurat"),
		    make_option(c("-t", "--timing"), action="store_true", default=FALSE, help="time the function run"))



opt = parse_args(OptionParser(option_list=option_list))

#args <- commandArgs(trailingOnly=T)
source('integration.R')
sobj = loadSeuratObject(opt$i)

if(opt$method=='seurat'){
	hvg<-2000
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
	}
	
	out = runSeurat(sobj, opt$b, hvg)
}

if(opt$method=='conos'){
	out = runConos(sobj, opt$b)
}

saveSeuratObject(out, opt$o)

