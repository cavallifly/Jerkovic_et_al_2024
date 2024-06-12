#!/bin/env Rscript

require(pastecs)
library(misha)
mDBloc <-  '/zdata/data/mishaDB/trackdb/'
db <- 'mm10'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()
source('/zdata/data/auxFunctions/auxFunctions.R')

#options(gmax.data.size = 1e+09)
options("scipen"=999)
conditions <- c("mESC_DAB","mESC_DA","mESC_DB","mESC_dCTCF","mESC_Dprom","mESC_polyA","mESC_UPPax6","mESC_WT","NPC_DAB","NPC_DA","NPC_DB","NPC_dCTCF","NPC_Dprom","NPC_polyA","NPC_WT")

for(insResolution in c(5e+3, 2e+3))
{
    #for(k in seq(150, 500, by=50))
    for(k in seq(300, 500, by=50))    
    {
	insScale <- k * 1e+3

	for(condition in conditions)
	{
	    print(paste0(insResolution," ",insScale))	
	    inFile=paste0("chic_",condition,"_merge_mm10_IJ_w",insScale/1000,"kb_r",insResolution/1000,"kb.tab")

	    insTrack <- read.table(inFile, header=T)
	    insTrackName <- paste0('insulation.',gsub('.tab','_cooler',inFile))

	    print(insTrackName)
	    #print(insTrack)

	    if(gtrack.exists(insTrackName)){print(paste0(insTrackName,' exists')); next}		

	    description=paste0("Loaded from ",inFile," from cooler")

	    gtrack.create_sparse(track=insTrackName, insTrack[,c(1,2,3)], insTrack[,c(4)], description=description)

	}
    }
}
warnings()
quit()
