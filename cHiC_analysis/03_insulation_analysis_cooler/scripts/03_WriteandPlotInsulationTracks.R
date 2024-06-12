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
celltypes <- c('NPC')
celltype <- 'NPC'
conditions <- c("mESC_WT_","mESC_UPPax6_","NPC_DAB_","NPC_DA_","NPC_DB_","NPC_Dprom_","NPC_WT_","NPC_dCTCF_","NPC_polyA_")



print(celltypes)
print(conditions)

#Extended Zfp608 regions
chrom <- "chr18"
#start <- 52598845
#end   <- 57599408
start <- 53859000
end   <- 56456400
#start <- 52775574
#end  <-  57266687
#start <- 53700000
#end  <-  56700000

# Make a plot of the insulation profile in the ROI
interv <- gintervals(chrom, start, end)


for(insResolution in c(5000,2000))
{
    for(k in c(300, 350, 400, 450, 500))
    {
	insScale <- k * 1e+3

	#nLines <- 0
	for(condition in conditions)	
	{
	    #chic_NPC_DAB_merge_mm10_IJ_w150kb_r2kb_cooler.track
 	    track <- paste0('insulation.chic_',condition,'merge_mm10_IJ_w',insScale/1000,'kb_r',insResolution/1000,'kb_cooler')
	    print(paste0("Plotting the insulation tracks for ",track))	    

	    if(gtrack.exists(track))
	    {
		#print(paste0(track," ",nLines))

		currentIns <- gextract(track, interv, colnames=c('ins'))
		currentIns$ins[is.na(currentIns$ins)] = 0
		print(head(currentIns))
		insMax <-  max(-currentIns$ins) 
		insMin <-  min(-currentIns$ins) 
		print(paste0(insMin,' ',insMax))
		currentIns$ins <- (-currentIns$ins - insMin) / (insMax - insMin)

		print(paste0("Find maxima on shoothed profile"))
		probThr <- 1e-3
		#probThr <- 1e-7

		#if(file.exists(outFile)){next}
		d <- smooth.spline(currentIns$ins, tol=1)
		tp <- turnpoints(d$y)
		insThr <-  mean(currentIns$ins) 
		tp.pos <- tp$tppos
		tp.proba <- tp$proba
		keep <- 1:(tp$nturns / 2) * 2
		if (tp$firstispeak) keep <- keep - 1
		tp.pos <- tp.pos[keep]
		tp.proba <- tp.proba[keep]    
		tp.pos <- tp.pos[tp.proba < probThr]
		tp.proba <- tp.proba[tp.proba < probThr]
		print(paste0("Take only the maxima whose insulation is above the mean"))
		tp.pos <- tp.pos[which(d$y[tp.pos] > insThr)]

		#outFile    = paste0('./insulation_profiles_cooler/Insulation_maxima_',gsub('insulation.','',track),'.bedGraph')
		#peaksDF <- data.frame()
		#peaksDF <- cbind("chr18",as.integer(d$x[tp.pos]*insResolution/insResolution)*insResolution+start,as.integer((d$x[tp.pos])*insResolution/insResolution+1)*insResolution+start,d$y[tp.pos])
		#peaksDF <- as.data.frame(peaksDF)
		#colnames(peaksDF) <- c("chrom","start","end","Insulation")
		#print(head(peaksDF))		   
		#write.table(peaksDF,file=outFile,row.names = F,quote = FALSE, col.names=F, sep="\t")

		outFile    = paste0('./insulation_profiles_cooler/Insulation_',gsub('insulation.','',track),'.bedGraph')
		peaksDF <- data.frame()
		peaksDF <- cbind("chr18",as.integer(d$x*insResolution/insResolution)*insResolution+start,as.integer((d$x)*insResolution/insResolution+1)*insResolution+start,d$y)
		peaksDF <- as.data.frame(peaksDF)
		colnames(peaksDF) <- c("chrom","start","end","Insulation")
		write.table(peaksDF,file=outFile,row.names = F,quote = FALSE, col.names=F, sep="\t")		
		#quit()
	    } 
	}
    }
}
warnings()
quit()
