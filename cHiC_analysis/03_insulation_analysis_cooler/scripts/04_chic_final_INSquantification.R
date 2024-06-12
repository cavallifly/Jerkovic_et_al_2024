require(misha)
require(shaman)
# require(plot3D)
# require(plotrix)
# require(pheatmap)

require(doParallel)

##########################################
### misha working DB
mDBloc <-  '/zdata/data/mishaDB/trackdb/'
db <- 'mm10'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

source('/zdata/data/auxFunctions/auxFunctions.R')
options(scipen=20,gmax.data.size=0.5e8,shaman.sge_support=1)

##############################################################################
tssCoordinates <- gintervals.load('intervals.ucscCanTSS')
rownames(tssCoordinates) <- tssCoordinates$geneName
geneCoordinates <- gintervals.load('intervals.ucscCanGenes')
rownames(geneCoordinates) <- geneCoordinates$geneName
##############################################################################

#ls -1 /zdata/data/mishaDB/trackdb/mm10/tracks/insulation/  | grep -v Nanog | grep -v Sox2 | grep 400kb

#refTrack <- "NPC_WTDStoNPCWT"

insTracks <- list(
	  mESC_DAB_merge    = "insulation.chic_mESC_DAB_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  mESC_DB_merge     = "insulation.chic_mESC_DB_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  mESC_DA_merge     = "insulation.chic_mESC_DA_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  mESC_dCTCF_merge  = "insulation.chic_mESC_dCTCF_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  mESC_Dprom_merge  = "insulation.chic_mESC_Dprom_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  mESC_polyA_merge  = "insulation.chic_mESC_polyA_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  mESC_UPPax6_merge = "insulation.chic_mESC_UPPax6_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  mESC_WT_merge     = "insulation.chic_mESC_WT_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  NPC_DAB_merge     = "insulation.chic_NPC_DAB_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  NPC_DA_merge      = "insulation.chic_NPC_DA_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  NPC_DB_merge      = "insulation.chic_NPC_DB_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  NPC_dCTCF_merge   = "insulation.chic_NPC_dCTCF_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  NPC_Dprom_merge   = "insulation.chic_NPC_Dprom_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  NPC_polyA_merge   = "insulation.chic_NPC_polyA_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler",
	  NPC_WT_merge      = "insulation.chic_NPC_WT_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp_cooler"
	)

#print(insTracks)

refTrack <- "NPC_WT_merge"

chr   <- 'chr18'
#start <- 52598845
#end   <- 57599408
start <- 53859000
end   <- 56456400
print(paste0("Considering insulation in ",chr,":",start,"-",end))

# Genomic locations of interest
C1 = gintervals("chr18",54986882,54991228)
C2 = gintervals("chr18",55130571,55137516)
C3 = gintervals("chr18",55194331,55199818)
C4 = gintervals("chr18",55544055,55552304)
C5 = gintervals("chr18",55758514,55765215)
C6 = gintervals("chr18",55811387,55818779)
targetPoints <- list(C1=C1,C2=C2,C3=C3,C4=C4,C5=C5,C6=C6)
#C1_offset = gintervals("chr18",54059846,54064192)
#C2_offset = gintervals("chr18",54203535,54210480)
#C3_offset = gintervals("chr18",54267295,54272782)
#C4_offset = gintervals("chr18",54617019,54625268)
#C5_offset = gintervals("chr18",54831478,54838179)
#C6_offset = gintervals("chr18",54884351,54891743)
#targetPoints <- list(C1=C1,C2=C2,C3=C3,C4=C4,C5=C5,C6=C6,C1_offset=C1_offset,C2_offset=C2_offset,C3_offset=C3_offset,C4_offset=C4_offset,C5_offset=C5_offset,C6_offset=C6_offset)

cols <- rainbow(length(names(insTracks)))

insIntrv <- gintervals(chr,start,end)

#####
### Insulation quantification
insQuantData <- list()
# 0. Compute insulation tracks using scripts/sandrine_computeInsulationTracks.R
# Retrieve the names of all insulation tracks 
allTracks <- gtrack.ls('insulation')
#print(allTracks)

for(set in names(insTracks))
{
    #set is the ID of the condition e.g. INS_hic_larvae_DWT_merge_dm6_BeS
    print(set)
    #Retrieve all tracks for the condition with different windows
    tracks     <- allTracks[grep(set,allTracks)]
    tracks     <- tracks[grep("_r2kb",tracks)]    
    #tracks     <- tracks[grep("_r5kb",tracks)]
    tracks     <- tracks[grep('cooler',tracks)]	    
    trackNames <- gsub('_Res','',gsub('INS_Scale_','',gsub('chic_','',gsub('_merge_mm10_IJ','',gsub('insulation.','',tracks)))))
    print(tracks)
    print(trackNames)    

    data <- gextract(tracks,insIntrv,iterator=2e3,colnames=trackNames)
    data[is.na(data)] <- 0
    #print(head(data))
		
    for(t in trackNames)
    {
	data[,t] <- clipQuantile(-data[,t],0.999)
	data[,t] <- scaleData(data[,t],0,1,1e-9)		
    }
    insQuantData[[set]] <- data
    #print(head(insQuantData))
}
#print(head(insQuantData))

#pdf(paste0("chic_Insulation_quantifications_at_Cpoints_5kb.pdf"),width=12,height=9)
pdf(paste0("chic_Insulation_quantifications_at_Cpoints_2kb.pdf"),width=12,height=9)
 
par(mfrow=c(2,2),mar=c(8,5,1.5,2))
print(names(insQuantData))

#print(targetPoints)
#quit()
for(target in names(targetPoints))
{
    targetPos <- targetPoints[[target]]
    print(paste0("Analysing ",target))
    maxneighbors = as.integer((targetPos$end - targetPos$start) / 2e3) + 1
    print(maxneighbors)

    print(targetPos)
    #print(head(insQuantData[[refTrack]]))
    #print(names(insQuantData[[refTrack]]))
    wtData <- gintervals.neighbors(targetPos,insQuantData[[refTrack]],maxneighbors=maxneighbors)[,grep('kb',colnames(insQuantData[[refTrack]]))+ncol(targetPos)]    
    print(gintervals.neighbors(targetPos,insQuantData[[refTrack]],maxneighbors=maxneighbors))
    npoints <- nrow(wtData)*ncol(wtData)
    print(paste0("Number of points ",nrow(wtData)*ncol(wtData)))


    #print(head(wtData))
    pvs <- c()
    pvalues <- c()
    #signs <- c()    
    for(set in names(insQuantData))
    {	
	mutData <- gintervals.neighbors(targetPos,insQuantData[[set]],maxneighbors=maxneighbors)[,grep('kb',colnames(insQuantData[[set]]))+ncol(targetPos)]
	print(paste0("Number of points ",nrow(mutData)*ncol(mutData)))
	print(set)

	print(t.test(wtData,mutData))
	pv <- -log(t.test(wtData,mutData)$p.value)
	pvalue <- t.test(wtData,mutData)$p.value
	#print(pv)
	pvs <- append(pvs,pv)
	pvalues <- append(pvalues,pvalue)
    }

    names(pvs) <- names(insQuantData)
    names(pvalues) <- names(insQuantData)

    cols <- c("blue","red","#b22222","violet","orange","gray50","black")
    barplot(pvs, main=paste0(target," ",npoints),las=2, ylim=c(0,50),col=cols, ylab='-log(Pval)')
    text('P=0.01',x=1,y=-log(0.01),pos=3)

    abline(h=-log(0.01))
    print(pvalues)
}
dev.off() 
