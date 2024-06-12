library("plyr")
#library("dplyr")
library("ggplot2")
library("reshape2")
library("seqplots")
#library(BSgenome)
source('auxFunctions.R')

chromSizes <- read.table("mm10.chrom.sizes")

### Get the list of bigWig files ###
bwDir <- './'
### here define if you want it to take all the bigwigs there (c("")) or specific bigwigs i.e. c("H3K4me1","H3K36me3")
targets <- c("ES") #c("H3K4me1","H3K36me3","H3K9ac","H3K27me3","H3K27ac","H3K4me3","H3K79me2","H2AFZ","H3K4me2","H3K9me1","H3K9me3","H4K20me1")
assembly <- 'mm10'
cellline <- 'Ivana_Pax6_ATAC_NPC_q0.01'

### Get the list of bigWig files ###
testBW <- list.files(bwDir, pattern = '.*bigWig') # tutti i bigWig
testBW <- testBW[ grep(paste(targets,collapse='|'), testBW) ] # Select by target and then by cell line
print(testBW)
### Get the list of bigWig files ###

### Set the parameters of the heatmap-plots and the clustering ###
# this is +- 5 kb of the smmit
range <- 5e3
binSize <- range/100
raster <- TRUE
meanBins <- 5 # 10 25
kMax <- 2 #3 # Number of clusters. I can change
col <- colorRampPalette(c('grey95','grey25'))(100)
### Set the parameters of the heatmap-plots and the clustering ###

### Load the peaks ###
### To change!!!load merged files that you made before
peakFile='merged_and_filtered_Pax6_peaks_q0.01_summits_bedtools.bed'
print(peakFile)
sourcePeaks <- read.table(peakFile)
colnames(sourcePeaks) <- c('chrom','start','end')
head(sourcePeaks)
### Load the peaks ###

### name of the pdf
name <- 'Ivana_Pax6_ATAC'

### Analysis of the peak statistics: length and chromosome ###
#print(paste0("Make the histogram of the silencers lenghts"))
#silencers=length1D(sourcePeaks)[,1:3]
#pdf(file=paste('Silencers_length_distribution','pdf',sep='.'))
#Length <- silencers$start
#par(xaxs='i',yaxs='i')
#head(Length)
#hist(Length, col="darkgreen", freq=TRUE, breaks=seq(from=0, to=15000, by=50), main="", xlab="Length (nt)", ylab="Number of silencers", cex.lab = 1.5)
#dev.off()
### Make the histogram of the silencers lenghts ###

# Make the histogram of the chromosomes involved
#pdf(file=paste(paste('Silencers','chromosomes','distribution',sep='_'),'pdf',sep="."))
#Chromosome <- sort(as.numeric(gsub('chr','',silencers$chrom)))
#par(xaxs='i',yaxs='i')
#barplot(table(Chromosome), col="lightgreen", horiz = TRUE, las=1, beside=TRUE, main="", xlab="Number of silencers", ylab="Chromosome", cex.lab = 1.5)
#dev.off()
### Analysis of the peak statistics: length and chromosome ###

datasets       <- list(silencers=center1D(sourcePeaks)[,1:3])
cols     <- list(silencers=colorRampPalette(c('grey90','darkblue'))(200))

cluster  <- c(TRUE,FALSE)[1]

# In this case we have just one peak set 'silencers'
for(peakSet in names(datasets))
{
    peaks <- datasets[[peakSet]]
    col <- cols[[peakSet]]
    #print(nrow(peaks))
    #print(ncol(peaks))
    singlePeakDataFile <- 'spData.bed'
    # Initialize the table with the peaks
    write.table(center1D(peaks[,1:3],1)[,1:3], file=singlePeakDataFile, quote=F, row.names=F, col.names=F, sep="\t")
    # Load the peaks. PlotSetArray is a function from the SeqPlots package which allows to load the bigWig files with the Chip-seq tracks (first argument: tracks) 
    # and associate them to the peakset (second argument: features), which are the silencers in this case, on the reference (refgenome), which is 'hg38' in this case.
    # It returns the heatmaps in sqData
    sqData <- getPlotSetArray(paste0(bwDir,testBW), singlePeakDataFile, refgenome=assembly, bin = binSize, rm0 = TRUE, ignore_strand = TRUE, xmin = range, xmax = range, type = "mf", add_heatmap = TRUE)

    print("Read the BW one-by one and prepare the data in allMeanData for the clustering")
    # allMeanData will have the targets as columns and the peaks as rows
    allData <- list()
    allMeanData <- c()
    for (set in testBW)
    {
        print(paste0("Analysing bigWig ",set,sep=""))

        target   <- unlist(strsplit(set,'_'))[1:3]
        sName <- paste(target,collapse='_')
        sName <- gsub('.bigWig','',sName)
        print(paste0("Plot column ID: ",sName,sep=""))
		
        # Get the heatmap
        # The following is a table for the bigWig set. It contains:
        # $anno    -> Genomic location of the peak : chr centre
        # $heatmap -> heatmap per peak of the target bigWig set
	hData <- sqData$data$spData[[which(testBW == set)]]$heatmap # hData is the heatmap for target set. It has Nrows=Npeaks and Ncolumns=Nbins
	hData[is.na(hData)] <- 0         # Substitute nans with zeros
	hData <- clipMatrix(hData,0.999) # Don't consider outliers

        # cData contains per each peak 1->Npeaks the elements the contribution of the target (in the bigWig) around the peak centre (ncol(hData)/2)) from -meanBins to +minBins
        # summed (rowSums) and log2-transformed. So, allMeanData will contain for each peak one value per target! 
        cData <- log2(rowSums(hData[,seq(round(ncol(hData)/2)-meanBins,round(ncol(hData)/2)+meanBins)])+1)
	if(is.null(allMeanData))
        {
            #if allMeanData is empty initilise it
            allMeanData <- matrix(cData, ncol=1)
	}else{
            #else add a new column (cbind) to allMeanData
	    allMeanData <- cbind(allMeanData,matrix(cData, ncol=1))
	}
	colnames(allMeanData)[ncol(allMeanData)] <- sName

        print(dim(hData))
	allData[[sName]] <- hData        # Store the heatmap in allData    
    }
    print(paste0("Dimension of allData ",dim(allData)))
    allMeanData[is.na(allMeanData)] <- 0 # Substitute nans with zeros in allMeanData.
	
    # layout divides the device (R 'terminal' for a figure) into as many rows and columns 
    # as there are in matrix mat, with the column-widths and the row-heights specified in the respective arguments.
    layout(matrix(1:length(testBW),ncol=length(testBW)))
		
    if(peakSet == names(datasets)[1])
    {
        # Since dist() computes the distance between rows, so transposing allMeanData, we cluster the targets (rows in the transposed) depending on the contribution
        # per each peak (columns in the transposed) using a hierarchical clustering
        clOrder <- hclust(dist(t(allMeanData)))$order
        print(paste0("Cluster order ",clOrder))
		
        # Using allMeanData as it is we now cluster the peaks (rows) depending on the contribution per target (columns) using a kmeans clustering.
	if(cluster == TRUE)
        {
	    set.seed(32)
            # kmeans cluster the peaks in kMax clusters using a kmeans algorithm depending on the value of each peak for the different targets
	    kMeans <- kmeans(allMeanData, kMax)
            # kOrder is the order of the peaks after the kmeans clustering
            # kMeans$cluster is a vector of integers (from 1:kMax) indicating the cluster to which each point is allocated.
            # order() orders the peaks of allMeanData depending on the cluster they belong.
            # That is: first all the peaks of cluster 1, then the ones of cluster 2 and so on...
            kOrder <- order(kMeans$cluster,apply(allMeanData,1,max))
            print(paste0("Order ",kOrder))
        }else{
	    print('NOT CLUSTERING')
	    kMeans <- km[[peakSet]]
            # 	    kOrder <- 1:nrow(peaks)
	    kOrder <- order(kMeans$cluster,apply(allMeanData,1,max))
	}
        peaks$k           <- kMeans$cluster
        
        # Number of heatmaps-panels (columns) in the final figure
        nPanels <- length(testBW)*length(names(datasets)) + 1 + 1
        # Width of the output figure
        outWidth <- nPanels*1.32
 	
        # Open the R device on which to save the figure as a .pdf
	pdf(file=paste(paste(name,paste('meanBins',meanBins,sep=""),paste('Nclusters',kMax,sep=""),cellline,sep='_'),'pdf',sep="."),width=outWidth,height=9)
		
	layout(matrix(1:(nPanels*2),ncol=nPanels), heights=c(outWidth/3.75,1))
        #par(mai=rep(0.125,4))
	#par sets the margins (mai=c(bottom, left, top, right)) and the stiles of y and x axis (Style "i" (internal) just finds an axis with pretty labels that fits within the original data range)
	par(mai=rep(c(0.3,0.125),2), yaxs="i", xaxs="i")
    }

    # Plot the leftmost column with the labels
    nameCex <- 1.25
    plot(0,xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1), col='white', frame=F)
    clLabels <- sapply(max(kMeans$cluster):1,function(k){
	      y <- (max(which(sort(kMeans$cluster, decreasing=FALSE) == k)) - length(which(kMeans$cluster == k))/2)/length(kMeans$cluster)
       	      txt <- paste0(paste0('CL',k),'\n',length(which(kMeans$cluster == k)),' peaks')
              text(0.5,y,txt,cex=nameCex)
	      return(txt)})
    addClusterLines(kMeans,FALSE,FALSE,col=1,max=length(kMeans$cluster))
    plot(0,xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,1), ylim=c(0,1), col='white', frame=F)
	
    print("Plot the columns of the figure one per target in the order given by hclust (clOrder)")
    print(dim(allData))
    for(set in names(allData)[clOrder])
    {
        # m is the heatmap for the target (e.g., histone mark H4K4me3) named set
        m <- allData[[set]]
        #print(dim(m))
        # Use clipMatrix to remove outliers    
        m <- clipMatrix(m,0.99)
        # Order the rows of m (one per peak) depending on the kmeans-cluster they belong
        # And take the transpose matrix (XXXWhy???XXX)
        m <- t(m[kOrder,])
        print(paste("Dimension of the heatmap ",dim(m)[[1]],dim(m)[[2]],sep=" "))
        
        # Set the title of the column
	title <- gsub('12787','12878',set)
        # Plot matrix m using color-scale col
	image(m, col=col, useRaster=raster, axes=F, main=gsub('H3','',title), cex.main=1.5)

        # Add black lines to mark the cluster borders
	addClusterLines(kMeans,FALSE,FALSE,col=1,max=length(kMeans$cluster))
        # Add yellow lines to mark the bins around the centre of the peak used for clustering
 	abline(v=(nrow(m)/2-meanBins)/nrow(m),lty=2,col='yellow',lwd=0.5)
 	abline(v=(nrow(m)/2+meanBins)/nrow(m),lty=2,col='yellow',lwd=0.5)

	if(cluster == TRUE)
        {
	    clOrd <- sort(kMeans$cluster)
	}else{
	    kMeans <- km[[peakSet]]
	    clOrd <- kMeans$cluster[kOrder]
	}

        print(title)	
        yMax <- max(unlist(sapply(1:max(kMeans$cluster),function(x){max(smooth(rowMeans(m[,which(clOrd == x)],na.rm=T)))})))
        print(title)
        yMin <- min(unlist(sapply(1:max(kMeans$cluster),function(x){min(smooth(rowMeans(m[,which(clOrd == x)],na.rm=T)))})))
	#yMin <- min(unlist(sapply(1:max(kMeans$cluster),function(x){min(smooth(rowMeans(m[,which(clOrd == x)],na.rm=T)))})))
        print(title)

        for (i in 1:max(kMeans$cluster))
        {
	    print(paste0("Cluster ",i))
	    clMtx <- m[,which(clOrd == i)]
	    if(i==1)
            {
		plot(smooth(rowMeans(clMtx,na.rm=T)), type='l', lwd=2, xlab='', ylab='', xaxt='n', yaxt='n',col=i,ylim=c(yMin,yMax),main=gsub('H3','',set))
		axis(1,labels=c(paste0('-',range/1e3,'kb'),0,paste0(range/1e3,'kb')), at=c(1,nrow(clMtx)/2,nrow(clMtx)))
	    }else{
                lines(smooth(rowMeans(clMtx,na.rm=T)), lwd=2,col=i)
		#lines(smooth(rowMeans(clMtx,na.rm=T)), lwd=2,col=i)
	    }
	}
    }
}

# Write the peaks per cluster
 dataframe <- as.data.frame(cbind(sourcePeaks$chrom,as.numeric(peaks$start[seq(1,length(sourcePeaks$chrom))]),peaks$end[seq(1,length(sourcePeaks$chrom))],kMeans$cluster))
 colnames(dataframe) <- c('chrom','start','end','k') 
 print(head(dataframe))
 for (i in 1:(max(kMeans$cluster)))
 {
    filename <- paste(paste(name,i,'MeanBins',meanBins,'Nclusters',kMax,cellline,sep='_'),'bed',sep=".")
    print(filename)
    write.table(subset(dplyr::filter(peaks,peaks$k == i),select=c("chrom", "start", "end", "k")), file=filename, quote=FALSE, sep='\t',row.names=FALSE)
}
dev.off()

