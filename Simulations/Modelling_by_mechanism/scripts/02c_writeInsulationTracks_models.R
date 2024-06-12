#!/bin/env Rscript

options(warn=-1)

#require(pastecs)
suppressWarnings(library(misha))
mDBloc <-  '../mishaDB/trackdb/'
db <- 'mm10'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

options("scipen"=999)

#Extended Zfp608 regions
chrom <- "chr18"
start <- 53700000
end   <- 56700000

cellType <- "XXXcellTypeXXX"
tag      <- "_XXXdirXXX_"

# Make a plot of the insulation profile in the ROI
interv <- gintervals(chrom, start, end)

for(insResolution in c(5000))
{
    for(k in c(400))
    {
	insScale <- k * 1e+3

	nLines <- 0
	for(condition in c(cellType))	
	{
	    for(track in gtrack.ls(paste0("insulation.models_",cellType),tag))
	    {
		currentIns <- gextract(track, interv, colnames='ins')

		insMax <-  max(-currentIns$ins) 
		insMin <-  min(-currentIns$ins) 
		currentIns$ins <- (-currentIns$ins - insMin) / (insMax - insMin)

		d <- smooth.spline(currentIns$ins, tol=1)
		outFile    = paste0('./XXXdirXXX/',gsub(paste0('insulation.models_',cellType,'.'),'',track),'.bed')
		if(file.exists(outFile)){next;}
		print(outFile)

		peaksDF <- data.frame()

		peaksDF <- cbind("chr18",currentIns$start,currentIns$end,currentIns$ins)
		peaksDF <- as.data.frame(peaksDF)
		colnames(peaksDF) <- c("chrom","start","end","Insulation")

		write.table(peaksDF,file=outFile,row.names = F,quote = FALSE, col.names=F, sep="\t")		
	    } 
	}
    }
}
warnings()
