getScriptPath <- function(){
	    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
        script.dir <- dirname(regmatches(cmd.args, m))
        if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
	    if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
	    return(script.dir)
}
# Setting the script path would then be:
setwd(getScriptPath())

library('optparse')

option_list <- list(make_option(c("-m", "--method"), type="character", default=NA, help="integration method to use"),
		    make_option(c("-i", "--input"), type="character", default=NA, help="input data"),
		    make_option(c("-o", "--output"), type="character", default=NA, help="output file"),
		    make_option(c("-b", "--batch"), type="character", default=NA, help="batch variable"),
		    make_option(c("-v", "--hvg"), type="character", default=NA, help="hvg list for seurat"),
		    make_option(c("-t", "--timing"), action="store_true", default=FALSE, help="time the function run"))



opt = parse_args(OptionParser(option_list=option_list))

#args <- commandArgs(trailingOnly=T)
#source(here('R','integration.R'))
source('../R/integration.R', chdir=T)
sobj = loadSeuratObject(opt$i)

if(opt$method=='seurat'){
	if(!is.na(opt$hvg)) {
		hvg<-unlist(readRDS(opt$hvg), use.names=FALSE)
	}
	else {
		hvg <- rownames(sobj@assays$RNA)
	}
	
	out = runSeurat(sobj, opt$b, hvg)
}

if(opt$method=='conos'){
	out = runConos(sobj, opt$b)
}

if(opt$method=='harmony'){
	
	out=runHarm(sobj, opt$b)
}
saveSeuratObject(out, opt$o)

