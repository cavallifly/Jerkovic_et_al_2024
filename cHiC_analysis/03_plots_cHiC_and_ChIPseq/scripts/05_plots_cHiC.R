options("scipen"=999, max.print=999999)

library(misha)
library(shaman)
library(plyr)
source('./scripts/auxFunctions.R')

mDBloc <- '../mishaDB/trackdb/'
db <- 'mm10'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

###Region to visualize
#Visualization: chr18:53,859,000-56,456,400
#Modelling region -> chr18:53700000-56700000
chrom <- "chr18"
start <- 53859000
end  <-  56456400

Capture    = 0 # 1 is FALSE 0 is TRUE
Genes      = 0
ViewPoints = 1 # 1 is FALSE 0 is TRUE
Insulation = 1 # 1 is FALSE 0 is TRUE
ChIPseqs   = 0 # 1 is FALSE 0 is TRUE

###Resolutions
chicRes    = 5000
virt4CRes  = 5000
insRes     = 5000

###Colors
insCols     =  c("#b22222")
#ChIPseq tracks      R45     G45      B45    A255
colChip = rgb(45 /255, 45/255, 45/255, alpha=255/255)
#RNAseq tracks       R194   G37      B92    A255
colRNA  = rgb(194/255, 37/255, 92/255, alpha=255/255)

###Initialization of variables
allViewPoints = list()

# Plots per condition
#celltypes  <- list('NPC','mESC')
celltypes  <- list('mESC')

for(ChIPres in c(5000))
{
    for(celltype in celltypes)
    {
        if(celltype == "NPC")
	{
	    conditions <- list('WT_') #,'DpromDStoNPCWT_','polyADStoNPCWT_','DADStoNPCWT_','DBDStoNPCWT_','DABDStoNPCWT_','dCTCFDStoNPCWT_','Dprom_','polyA_','DA_','DB_','DAB_','dCTCF_','WTDStoNPCWT_')

	    #chipDatasets <- c('Rad21','CTCF','H3K27ac','H3K4me3','H3K36me3','Pol2','Pax6','RNAseq')
    	    #plotOrder    <- c('Rad21','CTCF','H3K27ac','H3K4me3','H3K36me3','Pol2','Pax6','RNAseq')
    	    #plotTitles   <- c('Rad21','CTCF','H3K27ac','H3K4me3','H3K36me3','Pol2','Pax6','RNAseq')
	    chipDatasets <- c('Rad21','CTCF','H3K27ac','H3K4me3','H3K36me3','Pol2','RNAseq')
    	    plotOrder    <- c('CTCF','Rad21','H3K4me3','H3K27ac','H3K36me3','Pol2','RNAseq')
    	    plotTitles   <- c('CTCF','Rad21','H3K4me3','H3K27ac','H3K36me3','Pol2','RNAseq')	    

	    #Ivana's suggestion
	    plotMin      <- list(0,0,0,0,0,0,2,0)
	    if(ChIPres ==  500)
	    {
	    	plotMax      <- list(70,50,70,70,70,70,100)	
		#plotMax      <- list(70,100,100,100,70,100,50,100)
	    }
	    if(ChIPres == 1000)
	    {
		plotMax      <- list(70,100,100,100,70,100,50,100)
	    }	
	    if(ChIPres == 2500)
	    {
		plotMax      <- list(70,100,100,100,70,100,50,100)
	    }
	    if(ChIPres == 5000)
	    {
		#Ivanas plotMax      <- list(70,100,100,100,70,100,50,100)
		plotMax      <- list(40,40,70,70,30,50,50)
	    }
	    plotCols     = c(rep(colChip,length(plotTitles)-1),colRNA)
	}
	if(celltype == "mESC")
	{
	    #conditions <- list('WTDStomESCWT_')
	    ###Ivanas conditions <- list('UPPax6DStomESCWT_','UPPax6_','DpromDStomESCWT_','polyADStomESCWT_','DADStomESCWT_','DBDStomESCWT_','DABDStomESCWT_','dCTCFDStomESCWT_','Dprom_','polyA_','DA_','DB_','DAB_','dCTCF_','WTDStomESCWT_','Pax6recr_','Pax6recrDStoNPCWT_')
	    conditions <- list('UPPax6DStomESCWT_','UPPax6_','DpromDStomESCWT_','polyADStomESCWT_','DADStomESCWT_','DBDStomESCWT_','DABDStomESCWT_','dCTCFDStomESCWT_','Dprom_','polyA_','DA_','DB_','DAB_','dCTCF_','WTDStomESCWT_','Pax6recr_','Pax6recrDStoNPCWT_')	    
	
	    #chipDatasets <- c('CTCF','Rad21','H3K4me3','H3K27ac','H3K36me3','Pol2','Pax6','RNAseq')
    	    #plotOrder    <- c('CTCF','Rad21','H3K4me3','H3K27ac','H3K36me3','Pol2','Pax6','RNAseq')
    	    #plotTitles   <- c('CTCF','Rad21','H3K4me3','H3K27ac','H3K36me3','Pol2','Pax6','RNAseq')
	    chipDatasets <- c('CTCF','Rad21','H3K4me3','H3K27ac','H3K36me3','Pol2','RNAseq')
    	    plotOrder    <- c('CTCF','Rad21','H3K4me3','H3K27ac','H3K36me3','Pol2','RNAseq')
    	    plotTitles   <- c('CTCF','Rad21','H3K4me3','H3K27ac','H3K36me3','Pol2','RNAseq')	    

	    #Ivana's suggestion
	    plotMin      <- list(0,2,0,0,0,0,2,0)
	    if(ChIPres ==  500)
	    {
		plotMax      <- list(70,50,70,70,70,70,50,100)
	    }
	    if(ChIPres == 1000)
	    {
		plotMax      <- list(70,50,70,70,70,70,50,100)
	    }	
	    if(ChIPres == 2500)
	    {
		plotMax      <- list(70,50,70,70,70,70,50,100)
	    }
	    if(ChIPres == 5000)
	    {
		plotMax      <- list(40,40,70,70,30,50,50)	    
		###Ivanas plotMax      <- list(70,50,70,70,70,70,50,100)
	    }
	    plotCols     = c(rep(colChip,length(plotTitles)-1),colRNA)
	}

        for(condition in conditions)
        {
	    ChIPseqCondition = gsub('DStomESCWT','',gsub('DStoNPCWT','',condition))	
	    #print(ChIPseqCondition)

	    if(celltype == 'NPCNoPax6' && condition != 'WT_')
	    {
	        next
	    }

            if(ViewPoints == 0)
            {
	        #### List of all viepoints of the virtual 4C tracks
	        #viewpoint at the following points
	        #proximal promoter (P)     - deletion coordinates -    chr18:54,986,882-54,991,228
	        #alternative promoter (aP)                             chr18:54,991,244-54,993,455
	        #enhancer A (A)            - deletion coordinates -    chr18:55,544,055-55,552,304
	        #enhancer B (B)            - deletion coordinates -    chr18:55,811,387-55,818,779
	        allViewPoints <- list(
	        proximalPromoter    = list('chr18', 54986882, 54991228),
	        alternativePromoter = list('chr18', 54991244, 54993455),
	        enhancerA		  = list('chr18', 55544055, 55552304),
      	        enhancerB		  = list('chr18', 55811387, 55818779))
	        plotOrdervirt4C  <- c("proximalPromoter","alternativePromoter","enhancerA","enhancerB")	
	        plotTitlesvirt4C <- c("Proximal promoter","Alternative promoter","A","B")
	        #v4C tracks               R57      G106   B156  A255
	        plotColsvirt4C   <- c(rep(rgb(0.2235294, 0.4156863, 0.6117647, alpha=1.0),length(allViewPoints)))
	        plotMaxvirt4C    <- c(rep(100,length(allViewPoints)))
            }
            if(Insulation == 0)
            {
                ### Insulation ###
	        #insulation/chic_mESC_WT_merge_mm10_IJ_INS_Scale_400kb_Res_5000bp.track
	        insDatasets <- c('400kb')
    	        insOrder    <- c('400kb')
    	        insTitles   <- c(paste0('INS Scale 400kb Res 5000bp'))
    	        insMax      <- list(1.0)
		insTracks <- c(gtrack.ls('insulation',gsub('NoPax6','',celltype),paste0(condition,'_'),'Scale_400kb','_5000bp'))
		print(insTracks)
	    }
	    if(ChIPseqs == 0)
	    {
    	        ### ChIPseq ###
      	        print(chipDatasets)
	        #chipTracks <- c(gtrack.ls('chipseq',paste0(gsub('NoPax6','',celltype),"_",gsub('DStomESCWT','',gsub('DStoNPCWT','',condition))),'IJ'))	
	        chipTracks <- c(gtrack.ls('chipseq',paste0(gsub('NoPax6','',celltype)),"_",ChIPseqCondition,'IJ'))	
	        print(chipTracks)
	    }



       	    hicDatasets <- c(celltype)[1]	
    	    tag <- paste0(gsub('NoPax6','',celltype),'_',gsub("_","",condition))

	    print(tag)

            if(Capture == 0)
            {
 
 		interval2D <- gintervals.2d(chrom, start, end, chrom, start, end)
		binnedIterator <- giterator.intervals(intervals=interval2D, iterator=c(chicRes,chicRes))
		interval1D = gintervals(chrom, start, end)

		for(chicTrack in gtrack.ls('chic.',paste0(gsub('NoPax6','',celltype),"_",condition),'Scores'))
		{
        	    print(chicTrack)
	    	    if(length(chicTrack) == 0)
	    	    {
		        print(chicTrack)
		    	next
		    }

		    outfile <- paste0("./plots_shamanScores_Zfp608_",celltype,"_",gsub("_","",condition),"_",chicRes,"bp.png")
		    print(outfile)
		    if (file.exists(outfile)){next}

		    ##############################################################################
		    #plotRatios <- list(hic=1.5,genes=0.3,ins=0.5,virt4C=0.5,chip=0.2,clusters=0.075,scale=0.1)
		    plotRatios <- list(hic=0.5,genes=0.05,ins=1,virt4C=0.5,chip=0.05,clusters=0.075,scale=0.1)	    
		    band <- 1
		    ##############################################################################

		    nPlots <- length(hicDatasets) + length(chipDatasets) + 1 ### + 1 is for genes!

		    # PNG
		    width=850
		    png(outfile, width=width,height=(length(hicDatasets)*band*plotRatios[['hic']]+plotRatios[['genes']]+length(chipDatasets)*plotRatios[['chip']])*width, type = "cairo")
		    lM <- 0.4 # 3.5
		    bM <- (0.2*band)/2
		    tM <- (0.2*band)/2
		    rM <- lM
		    #mai		    A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
		    par(mai=c(bM,lM,tM,rM), xaxs="i", yaxs="i", family='sans', cex=1)	    

		    layout(matrix(1:nPlots, ncol=1), height=c(rep(plotRatios[['hic']],length(hicDatasets))*band, rep(plotRatios[['genes']],1), rep(plotRatios[['virt4C']],length(allViewPoints)), rep(plotRatios[['ins']],0), rep(plotRatios[['chip']],length(chipDatasets)), width=2, respect=T))


		    if(!is.null(chicTrack))
		    {
	
			print("##### Plot triangular matrix #####")
			data <- gextract(chicTrack, binnedIterator, iterator=binnedIterator,colnames = c('score'))

			# Plot triangular matrix		        
			#pointCEX <- 0.05 # Res 5kb
			pointCEX <- 0.45 # Res 10kb	
			if(nrow(data) < 50000){pointCEX <- 0.5;}

			#data <- read.table(infile, header=T)

			plotData <- data
			print(head(plotData))
			print(tail(plotData))		
			plotData <- plotData[plotData$start2 >= plotData$start1,]
	    		plotData <- plotData[order(plotData$score),]
			plotData$score[is.na(plotData$score)] = 0

			max <- 100 #cHi-C balanced and scores
			#max <- round_any(max(max(plotData$score),-min(plotData$score))*0.8,accuracy=1) #cHi-C balanced and scores
			print(max)
		
			scoreCol <- colorRampPalette(c('#000080','#ffffff','#8b0000'))(200)
			pal <- colorRampPalette(c('darkblue','white','darkred'))
			#print(scoreCol)

    			xAxis <- c(interval2D$start1,interval2D$end2)
	    		winSize <- xAxis[2]-xAxis[1]
    			yAxis <- c(0,winSize)

    			plot(0,xlim=xAxis,ylim=yAxis, col='white',cex.main=2,xlab='', ylab='', yaxt='n', xaxt='n', frame=F)
			maxDist <- 0
			apply(plotData, 1, function(x){
	        	s <- as.numeric(as.character(x['score']))+101
		    
		        col <- scoreCol[s]
		        #print(paste0(col," ",s))

		        start <- as.numeric(as.character(x['start1']))
		    	end   <- as.numeric(as.character(x['start2']))
	    	        dist  <- end-start
		        if(dist>maxDist){maxDist = dist}

		    	mid   <- min(c(start,end))+abs(dist)/2
		    
	   	        xx    <- mid
		        yy    <- abs(dist)

		        points(x=xx, y=yy, type='p', cex=pointCEX, col=col, pch=15)
		        })
			#rng <- paste0("Range ", round_any(min(plotData$score),accuracy=1)," - ",round_any(max(plotData$score),accuracy=1))
			#legend("topleft",legend=rng, cex=1.6)
			legend("topleft",legend=paste0(celltype,' ',gsub('_|DStomESCWT_','',condition)), cex=1.6, box.col = "white")
			#text(x=interval2D$end1-interval2D$start1, y=maxDist/2, labels=paste0(celltype,'_',condition),pos=2, cex=1.4, col='black')
			#text(x=interval2D$end1-interval2D$start1, y=maxDist/2, labels=expression(min(plotData$score),"-",max(plotData$score)),pos=2, cex=1.4, col='black')
		    }
	    
		    if(Genes == 0)
		    {
		        print("# Plot gene annotation")
		        ### Genes
	   	        ### Load TSS and GENE coordinates...
	  	        tssCoordinates <- gintervals.load('intervals.ucscCanTSS') # Ritorna un dataframe
		        rownames(tssCoordinates) <- tssCoordinates$geneName
	  	        geneCoordinates <- gintervals.load('intervals.ucscCanGenes')
	  	        rownames(geneCoordinates) <- geneCoordinates$geneName
	  	        genes <- gintervals.neighbors(geneCoordinates,interval1D) # Non e' 100% sicuro che funzioni bene. E' un intersect tra intervalli e geni presenti.
		        TSSs  <- gintervals.neighbors(tssCoordinates,interval1D) # Non e' 100% sicuro che funzioni bene. E' un intersect tra intervalli e geni presenti.
                                                        # Puoi impostare il numero di primi vicini...o usare massimo distance minimum distance: centre+/-2kb
                                                        # Usa un altro exact_neighbors che e' basato su stringhe!

		        genesPlot(genes,plotLim=c(interval1D$start,interval1D$end),cex=1.0,rHeight=90) # in auxiliary function!
		    }

		    if(Insulation == 0)
		    {
			if(length(insTracks) > 0)
		        {
		            print("# Plot insulation tracks")
		            insData <- gextract(insTracks,intervals=interval1D,iterator=insRes, colnames = gsub('insulation.|merge_','',insTracks))
	  		    print(head(insData))

			    plotIns(chipData=insData,plotOrder=insOrder,plotCols=insCols,plotTitles=insTitles,plotStart=interval1D$start,plotEnd=interval1D$end,chipRes=insRes,plotMax=insMax,chrom=chrom,main=FALSE,condition=tag, chipStart=interval1D$start)
 	        	}
		    }

	            if(ViewPoints == 0)
                    {
		        ##### Plot virtual4C track
		        virt4CData <- c()
		        for(viewPoint in names(allViewPoints))
		        {
	  	            print(viewPoint)
		    	    locus <- allViewPoints[[viewPoint]]
	 	 	    locus <- unlist(locus)
			    print(locus)

		  	    interval2Dvirt4C <- gintervals.2d(chrom, start, end, locus[1], as.numeric(locus[2]), as.numeric(locus[3]))
			    print(interval2Dvirt4C)
			    binnedIterator <- data.frame(matrix(ncol = 6, nrow = as.integer((end-start)/virt4CRes)+2))
			    colnames(binnedIterator) <- c("chrom1","start1","end1","chrom2","start2","end2")
			    binnedIterator$chrom1	 <- chrom
			    binnedIterator$chrom2	 <- chrom
			    print(colnames(binnedIterator))
			    nbin = 0
			    s1 = as.integer(start/virt4CRes+1)*virt4CRes
			    s2 = as.integer(end/virt4CRes)*virt4CRes
			    resolution1 = virt4CRes
			    resolution2 = as.numeric(locus[3])-as.numeric(locus[2])
			    print(s1)
			    print(s2)		
			    print(resolution1)
			    print(resolution2)		
		
			if(s1>start)
			{
			    nbin = nbin+1
			    binnedIterator$start1[nbin] <- start
			    binnedIterator$end1[nbin]   <- s1
			    binnedIterator$start2[nbin] <- as.numeric(locus[2])
			    binnedIterator$end2[nbin] <- as.numeric(locus[3])			    		    
			}

			for(binStart1 in seq(from=s1, to=s2-resolution1, by=resolution1))
			{
			    for(binStart2 in seq(from=as.numeric(locus[2]), to=as.numeric(locus[3])-1, by=resolution2))
			    {
				nbin = nbin+1
				#print(paste(nbin,binStart1,binStart1+resolution1,binStart2,binStart2+resolution2))
				binnedIterator$start1[nbin] <- binStart1
				binnedIterator$end1[nbin]   <- binStart1+resolution1
				binnedIterator$start2[nbin] <- binStart2
				binnedIterator$end2[nbin]   <- binStart2+resolution2-1
			    }
			}
			if(s2<end)
			{
			    nbin = nbin+1
			    binnedIterator$start1[nbin] <- s2
			    binnedIterator$end1[nbin]   <- end
			    binnedIterator$start2[nbin] <- as.numeric(locus[2])
			    binnedIterator$end2[nbin]   <- as.numeric(locus[3])			    		    
			}
			#print(binnedIterator)
			#binnedIterator <- giterator.intervals(intervals=interval2Dvirt4C, iterator=c(virt4CRes,as.numeric(locus[3])-as.numeric(locus[2])))
			#quit()

			virt4CTrack <- gextract(chicTrack, binnedIterator, iterator=binnedIterator, colnames= c(viewPoint))
			#print(virt4CTrack)
			print(length(virt4CTrack))
			print(nrow(virt4CTrack))

			if(is.null(virt4CData))
			    {
			        virt4CData <- virt4CTrack[,c("chrom1","start1","end1",viewPoint)]
			    } else {
			        virt4CData <- cbind(virt4CData,virt4CTrack[c(viewPoint)])
		   	    }

		  	    virt4CData[is.na(virt4CData)] <- 0
		        }

		        outfile <- paste0("./virtual4C_shamanScores_Zfp608_",celltype,"_",gsub("_","",condition),"_chic_",chicRes,"bp.tsv")	    	    
		        virt4CData[is.na(virt4CData)] <- 0
		        print(virt4CData)	    

		        print(head(virt4CData))	

		        plotAnn=allViewPoints
		        print(plotAnn)

		        plotVirt4C(virt4CData,plotOrder=plotOrdervirt4C,plotCols=plotColsvirt4C,plotTitles=plotTitlesvirt4C,plotStart=interval1D$start,plotEnd=interval1D$end,chipRes=virt4CRes,plotMax=plotMaxvirt4C,chrom=chrom,main=FALSE,plotAnn=plotAnn)
	            }

                    if(ChIPseqs == 0)
                    {

                        if(length(chipTracks) > 0)
                        {
                            chipData <- gextract(chipTracks,intervals=interval1D,iterator=ChIPres, colnames = gsub('chipseq.|_merge_|mm10_IJ|NPC_dCTCF|_mESC_UPPax6|_mESC_dCTCF','',chipTracks))       
			    
                            print("# Plot ChIPseq tracks")
                            print(chipTracks)
			    print(head(chipData))

			    #plotAnn=read.table("./scripts/deletions.bed", header=T)
                            #plotAnn=read.table("./scripts/contacts.bed", header=T)             
			    #plotAnn=read.table("./scripts/Bonev_TADs_ES.bed", header=T)             

                            CTCFDir <- "/work/cavalli/mdistefano/2021_04_13_Project_with_Ivana/CTCF_peaks_directionality/pwmscan_mm10_19159_19671_all_possible_CTCF_motifs.bed"
		  	    plotChIP(chipData,plotOrder=plotOrder,plotCols=plotCols,plotTitles=plotTitles,plotStart=interval1D$start,plotEnd=interval1D$end,chipRes=ChIPres,plotMax=plotMax,chrom=chrom,main=FALSE,plotMin=plotMin)
                        }
		    }
		    dev.off()
	        }     
            }
	    #quit()

	}

    }

    pdf(paste0('scoreMaps_colorbar.pdf'),width=2,height=8)
    labels <- c("-100","-50","0","50","100")
    scoreCol <- colorRampPalette(c('#000080','#ffffff','#8b0000'))(200)
    x <- length(scoreCol)
    plot(y=1:x,x=rep(0,x),col=scoreCol, pch=15, xlab='', ylab='', yaxt='n', xaxt='n',frame=F,ylim=c(-x,x),xlim=c(-5,5), main='');
    sapply(1:length(labels), function(y){text(y=seq(0,1,length.out=length(labels))[y]*x,x=0,labels=labels[y],pos=2)});
    dev.off()

}
warnings()
quit()
