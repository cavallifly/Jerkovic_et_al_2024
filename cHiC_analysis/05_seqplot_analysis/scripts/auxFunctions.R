# require(misha)
# require(shaman)
# require(plyr)
# require(svMisc)
# require(reshape2)
# require(pheatmap)
# require(plotrix)
# require(pastecs)


##########################
snapToGrid <- function(source,target,d=10e3){
	
# 	t <- gintervals(target$chrom1,target$start1,target$start1+1)
	t <- g2anchors(target)
	
	ret <- source
	for(i in 1:nrow(source)){
		
		p1 <- gintervals(source[i,'chrom1'],source[i,'start1']-1,source[i,'start1']+1)
		p2 <- gintervals(source[i,'chrom1'],source[i,'end1']-1,source[i,'end1']+1)
		pair1 <- gintervals.neighbors(p1,t,maxneighbors=1,mindist=-d,maxdist=d)
		pair2 <- gintervals.neighbors(p2,t,maxneighbors=1,mindist=-d,maxdist=d)
		
		if(!is.null(pair1)){
			ret[i,'start1'] <- pair1[1,5]
			ret[i,'start2'] <- pair1[1,5]
		}
		
		if(!is.null(pair2)){
			ret[i,'end1'] <- pair2[1,5]
			ret[i,'end2'] <- pair2[1,5]
		}
	}
# 	print(ret)
	return(ret)
	
	next
	s1 <- gintervals(source$chrom1,source$start1,source$start1+1)
	s2 <- gintervals(source$chrom1,source$end1,source$end1+1)
	
	
	pair1 <- gintervals.neighbors(s1,t1,maxneighbors=1,mindist=-d,maxdist=d)
	pair2 <- gintervals.neighbors(s2,t2,maxneighbors=1,mindist=-d,maxdist=d)
	
	chrs <- pair1[,1]
	starts <- pair1[,5]
	ends <- pair2[,5]
	
	retData <- gintervals.2d(chrs,starts,ends,chrs,starts,ends)
	
}

shufflePositions <- function(points, minDist=50e3, maxDist=500e3, binSize=10e3){
	
	shifts <- apply(points,1,function(x){sample(seq(minDist,maxDist,by=binSize),1)})
	signs <- apply(points,1,function(x){sample(c(-1,1),1)})
	points$start <- points$start + shifts*signs
	
# 	chrEnds <- gintervals.all()
	
	for(i in 1:nrow(points)){
		if(points[i,'start'] > gintervals.all()[which(as.character(gintervals.all()$chrom) == as.character(points[i,'chrom'])),'end']){
			points[i,'start'] <- gintervals.all()[which(as.character(gintervals.all()$chrom) == as.character(points[i,'chrom'])),'end'] - sample(2:100,1)
		}
	}
	points[which(points$start < 0),'start'] <- 1
	points$end <- points$start + 1
	
	return(points)
}


pileUpBed <- function(intrv,data,bin=1e3,od=0){
	outBed <- c()
	iter <- giterator.intervals(intervals=intrv,iterator=bin)
	for(i in 1:nrow(iter)){
		binScore <- gintervals.neighbors(iter[i,],data,mindist=-od,maxdist=od,maxneighbors=5e3)
		if(is.null(binScore)){
			outBed <- append(outBed,0)
		}else{
			outBed <- append(outBed,nrow(binScore))
		}
	}
	iter$score <- outBed
	return(iter)
}

##########################
create.marginal.intervals=function(width, chroms=NULL)
{
  fc = gintervals.all()

  if (!is.null(chroms)) {
    if (sum(grep("chr", chroms[1])) == 0) {
      chroms = paste("chr", chroms, sep="")
    }
    fc = fc[is.element(fc$chrom, chroms),]
  }

  ints = data.frame()
  for (i in 1:nrow(fc)) {
    starts = seq(fc[i, 'start'], fc[i, 'end'], by=width)
    ends = pmin(starts + width, fc[i, 'end'])
    ints = rbind(ints, data.frame(chrom=fc[i, 'chrom'], s1=starts, e1=ends, s2=fc[i, 'start'], e2=fc[i, 'end']))
  }
  gintervals.2d(chroms1=ints$chrom, starts1=ints$s1, ends1=ints$e1, starts2=ints$s2, ends2=ints$e2)
}

########################
get.marginal=function(tn, width=10^6, chroms=c(1:19, 'X'))
{
  tns = paste(tn, "area", sep="_")
  gvtrack.create(tns, tn, "area")
  ints = create.marginal.intervals(width, chroms=chroms)
  gextract(tns, intervals=ints, iterator=ints)
}


addSomCluster <- function(m,col=1,lwd=2){
	for(i in 1:nrow(m)){
		for(j in 1:ncol(m)){
			center <- c(i,j)
# 			if(m[center[1],center[2]] != 0){centers[[paste0('CL',m[center[1],center[2]])]] <- rbind(centers[[paste0('CL',m[center[1],center[2]])]],matrix(center,ncol=2))}
			r <- c(i+1,j)
			tr <- c(i+1,j+1)
			if(center[2] %% 2 == TRUE){tr <- c(i,j+1)}
			tl <- c(i,j+1)
			if(center[2] %% 2 == TRUE){tl <- c(i-1,j+1)}
			br <- c(i+1,j-1)
			if(center[2] %% 2 == TRUE){br <- c(i,j-1)}
			bl <- c(i,j-1)
			if(center[2] %% 2 == TRUE){bl <- c(i-1,j-1)}
			
			l <- c(i-1,j)
			
			c <- 0.5
			xR <- (1/(nrow(m)+c))/2
			yR <- (1/(nrow(m)+c))/3
			
			if(r[1] > nrow(m)){next}
			
			if(m[center[1],center[2]] != m[r[1],r[2]]){
				d <- ifelse(center[2] %% 2 == TRUE,-0.5,0)
				x <- (center[1]+d)/(nrow(m)+c)
				y <- (center[2]-c)/ncol(m)
				segments(x+xR,y-yR,x+xR,y+yR,lwd=lwd,col=col)
				
				if(tr[1] <= nrow(m) & tr[2] <= ncol(m) & tr[1] > 0 & tr[2] > 0){
					if(m[center[1],center[2]] != m[tr[1],tr[2]]){segments(x,y+yR*2,x+xR,y+yR,lwd=lwd,col=col)}
				}
				if(tl[1] <= nrow(m) & tl[2] <= ncol(m) & tl[1] > 0 & tl[2] > 0){
					if(m[center[1],center[2]] != m[tl[1],tl[2]]){segments(x,y+yR*2,x-xR,y+yR,lwd=lwd,col=col)}
				}
				if(br[1] <= nrow(m) & br[2] <= ncol(m) & br[1] > 0 & br[2] > 0){
					if(m[center[1],center[2]] != m[br[1],br[2]]){segments(x,y-yR*2,x+xR,y-yR,lwd=lwd,col=col)}
				}
				if(bl[1] <= nrow(m) & bl[2] <= ncol(m) & bl[1] > 0 & bl[2] > 0){
					if(m[center[1],center[2]] != m[bl[1],bl[2]]){segments(x,y-yR*2,x-xR,y-yR,lwd=lwd,col=col)}
				}
			}
			
			if(l[1] < 1){next}
			if(m[center[1],center[2]] != m[l[1],l[2]]){
				d <- ifelse(center[2] %% 2 == TRUE,-0.5,0)
				x <- (center[1]+d)/(nrow(m)+c)
				y <- (center[2]-c)/ncol(m)
				segments(x-xR,y-yR,x-xR,y+yR,lwd=lwd,col=col)
				if(tr[1] <= nrow(m) & tr[2] <= ncol(m) & tr[1] > 0 & tr[2] > 0){
					if(m[center[1],center[2]] != m[tr[1],tr[2]]){segments(x,y+yR*2,x+xR,y+yR,lwd=lwd,col=col)}
				}
				if(tl[1] <= nrow(m) & tl[2] <= ncol(m) & tl[1] > 0 & tl[2] > 0){
					if(m[center[1],center[2]] != m[tl[1],tl[2]]){segments(x,y+yR*2,x-xR,y+yR,lwd=lwd,col=col)}
				}
				if(br[1] <= nrow(m) & br[2] <= ncol(m) & br[1] > 0 & br[2] > 0){
					if(m[center[1],center[2]] != m[br[1],br[2]]){segments(x,y-yR*2,x+xR,y-yR,lwd=lwd,col=col)}
				}
				if(bl[1] <= nrow(m) & bl[2] <= ncol(m) & bl[1] > 0 & bl[2] > 0){
					if(m[center[1],center[2]] != m[bl[1],bl[2]]){segments(x,y-yR*2,x-xR,y-yR,lwd=lwd,col=col)}
				}
			}
			
			if(tr[1] <= nrow(m) & tr[2] <= ncol(m) & tr[1] > 0 & tr[2] > 0){
				if(m[center[1],center[2]] != m[tr[1],tr[2]]){
					d <- ifelse(center[2] %% 2 == TRUE,-0.5,0)
					x <- (center[1]+d)/(nrow(m)+c)
					y <- (center[2]-c)/ncol(m)
					segments(x,y+yR*2,x+xR,y+yR,lwd=lwd,col=col)
				}
			}
			
			if(tl[1] <= nrow(m) & tl[2] <= ncol(m) & tl[1] > 0 & tl[2] > 0){
				if(m[center[1],center[2]] != m[tl[1],tl[2]]){
					d <- ifelse(center[2] %% 2 == TRUE,-0.5,0)
					x <- (center[1]+d)/(nrow(m)+c)
					y <- (center[2]-c)/ncol(m)
					segments(x,y+yR*2,x-xR,y+yR,lwd=lwd,col=col)
				}
			}
		}
	}
}

prepareCrossTracks <- function(trackName, vType, res, cType=2){
	if(cType == 2){
		checkers <- list(i1=c(-1,0,0,1),i2=c(0,1,0,1),i3=c(-1,0,-1,0),i4=c(0,1,-1,0),o1=c(2,3,2,3),o2=c(-3,-2,-3,-2),o3=c(-3,-2,2,3),o4=c(2,3,-3,-2))
	}else{
		checkers <- list(i1=c(-1,0,0,1),i2=c(0,1,0,1),i3=c(-1,0,-1,0),i4=c(0,1,-1,0),o1=c(1,2,1,2),o2=c(-2,-1,-2,-1),o3=c(-2,-1,1,2),o4=c(1,2,-2,-1))
	}
	
	vNames <- c()
	for(vName in names(checkers)){
		print(vName)
		pattern <- checkers[[vName]]
		if(length(grep(vName,gvtrack.ls())) > 0){gvtrack.rm(vName)}
		gvtrack.create(vName, trackName, vType)
		
		gvtrack.iterator.2d(vName, 
		sshift1=pattern[1]*res, eshift1=pattern[2]*res,
		sshift2=pattern[3]*res, eshift2=pattern[4]*res)
		vNames <- append(vNames, vName)
		
	}
	return(vNames)
}

my_hex_grid <- function(dim=80, col=c("darkblue","blue","white", "red", "darkred"), border=FALSE){
	
	plot.new()
	
	centers.x <- seq(0, 1, length.out = dim*2+2)
	centers.x <- rep(c(centers.x[which(1:(dim*2+2) %% 2 ==0)[-(dim+1)]], centers.x[(which(1:(dim*2+2) %% 2 ==0)+1)[-(dim+1)]]), dim/2)
	centers.y <- rep(seq(0, 1, length.out = dim*2+1)[which(1:(dim*2) %% 2 ==0)], each=dim)
	size.x <- 1/(dim*2+1)
	size.y <- 1/(dim*3+2)
	
	for(i in 1:length(centers.x)){
		x <- c(centers.x[i], centers.x[i]+size.x, centers.x[i]+size.x, centers.x[i], centers.x[i]-size.x, centers.x[i]-size.x)
		y <- c(centers.y[i]+size.y*2, centers.y[i]+size.y, centers.y[i]-size.y, centers.y[i]-size.y*2, centers.y[i]-size.y, centers.y[i]+size.y)
		if(border==FALSE){
			polygon(x, y, border=NA, col=col[i])
		}else{
			polygon(x,y, col=col[i],border=col[i], ljoin=1,lwd=0.05)
		}
	}
# 	axis(1)
# 	axis(2)
}

addClusterLines <- function(data,rev=FALSE,image=TRUE,kmeans=TRUE,max=0, col='darkred',lwd=2.5,lty=1){
	if(kmeans == TRUE){
		start <- 0
		for (i in 1:max(data$cluster)){
			current <- start + length(data$cluster[data$cluster == i])
			segments(0,current/length(data$cluster),1,current/length(data$cluster), lwd=lwd, col=col, lty=lty, lend=1)
			start <- current
		}
	}else{
		start <- 0
		for (i in data){
			current <- start + i
			abline(h=current/max, lwd=lwd, col=col, lty=lty)
			start <- current
		}
	}
}

# mergeTouching2D <- function(data){
# 	
# 	forData <- c()
# 	for (chr in unique(as.character(data$chrom1))){
# 		inData <- data[data$chrom1 == chr,1:6]
# 		print(chr)
# 		for(set in unique(paste0(inData$start1,'_',inData$end1))){
# 			cData <- inData[paste0(inData$start1,'_',inData$end1) %in% set,]
# 			cData <- cData[order(cData$start2),]
# 			
# 			s1 <- cData[1,2]
# 			e1 <- cData[1,3]
# 			s2 <- cData[1,5]
# 			e2 <- cData[1,6]
# 			
# 			if(nrow(cData) > 1){
# 				for(i in 2:nrow(cData)){
# 					p1 <- cData[i-1,'start2']
# 					p2 <- cData[i-1,'end2']
# 					p3 <- cData[i,'start2']
# 					p4 <- cData[i,'end2']
# 					
# 					if(p2 != p3){
# 						forData <- rbind(forData,gintervals.2d(chr,s1,e1,chr,s2,p2))
# 						s2 <- p3
# 					}
# 					if(i == nrow(cData)){forData <- rbind(forData,gintervals.2d(chr,s1,e1,chr,s2,p4))}
# 				}
# 			}else{
# 				forData <- rbind(forData,gintervals.2d(chr,s1,e1,chr,s2,e2))
# 			}
# 		}
# 	}
# 	
# 	outData <- c()
# 	for (chr in unique(as.character(forData$chrom1))){
# 		inData <- forData[forData$chrom1 == chr,1:6]
# 		print(chr)
# 		for(set in unique(paste0(inData$start2,'_',inData$end2))){
# 			cData <- inData[paste0(inData$start2,'_',inData$end2) %in% set,]
# 			cData <- cData[order(cData$start1),]
# 			
# 			s1 <- cData[1,5]
# 			e1 <- cData[1,6]
# 			s2 <- cData[1,2]
# 			e2 <- cData[1,3]
# 			
# 			if(nrow(cData) > 1){
# 				for(i in 2:nrow(cData)){
# 					p1 <- cData[i-1,'start1']
# 					p2 <- cData[i-1,'end1']
# 					p3 <- cData[i,'start1']
# 					p4 <- cData[i,'end1']
# 					
# 					if(p2 != p3){
# 						outData <- rbind(outData,gintervals.2d(chr,s2,p2,chr,s1,e1))
# 						s2 <- p3
# 					}
# 					if(i == nrow(cData)){outData <- rbind(outData,gintervals.2d(chr,s2,p4,chr,s1,e1))}
# 				}
# 			}else{
# 				outData <- rbind(outData,gintervals.2d(chr,s2,e2,chr,s1,e1))
# 			}
# 		}
# 	}
# 	return(outData)
# }


peaksPlot <- function(peaks, intrv, col=NA, add=FALSE){
	yLim=c(-100,100)
	plotLim <- c(intrv[1,2],intrv[1,3])
	
	if(add==FALSE){
		plot(0, xlim=plotLim, ylim=yLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F, col='white')
		segments(plotLim[1],0,plotLim[2],0,lwd=1.75, col='darkgrey')
	}
	
	if(is.null(peaks)){return()}
	
	for(g in 1:nrow(peaks))
	{
		row <- peaks[g,]
		if(is.na(col)){
			rect(row$start, yLim[1], row$end, yLim[2], col=rgb(0.66,0.66,0.66, alpha=0.4), border = 0)
		}else{
			rect(row$start, yLim[1], row$end, yLim[2], col=col, border = 0)
		}
	}
}


### provide LIST
jitterPlot <- function(data,cols=NA,main='',ylim=NA,pTarget=NA,points=TRUE,cexPV=1,valPV=TRUE,jitWidth=0.5,pvalues=TRUE,test='wilcox'){
	
	if(!is.list(data)){
		print('PLEASE PROVIDE A NAMED LIST')
		return()
	}
	
	bCols <- rep('grey40', length(data))
	if(!is.na(pTarget)){bCols[which(names(data) == pTarget)] <- 'blue'}
	
	if(is.na(ylim[1])){
		plt <- boxplot(data, notch=T, frame=F, outline=F, main=main, col=bCols,las=2)
		yMax <- max(plt$stats)*0.85
	}else{
		boxplot(data, notch=T, frame=F, outline=F, main=main, col=bCols,ylim=ylim,las=2)
		yMax <- ylim[2]*0.85
	}
	abline(h=0, lwd=1.5,col='black',lty=2)
	for (i in 1:length(data)){
		Y <- data[[i]]
		X <- rep(i,length(Y))
		ifelse(!is.na(cols),col <- cols[i],col <- i+1)
		if(points ==TRUE){
			points(jitter(X,length(data)/i*jitWidth),Y,pch=19,cex=0.25, col=col)
		}
		
		if(pvalues == TRUE){
			pPlot <- FALSE
			if(is.na(pTarget)){
				
# 				if(i == 2){next}
				
				if(i < length(data)){
					if(test == 'wilcox'){
						pval <- wilcox.test(data[[i]],data[[i+1]])$p.value
					}else{
						pval <- ks.test(data[[i]],data[[i+1]])$p.value
					}
					pPlot <- TRUE
				}
			}else{
				if(test == 'wilcox'){
					pval <- wilcox.test(data[[i]],data[[pTarget]])$p.value
				}else{
					pval <- ks.test(data[[i]],data[[pTarget]])$p.value
				}
				if(names(data)[i] != pTarget){pPlot <- TRUE}
				i <- i-0.5
			}
			
			if(pPlot == TRUE){
				if(valPV == TRUE){text(x=i+0.5,y=yMax*1.1,signif(pval,2),cex=cexPV); ySig <- yMax
				}else{ySig <- yMax*1.1}
				
				if(pval <= 1e-6){text(x=i+0.5,y=ySig,'***',cex=cexPV*1.1,col=2)
				}else if(pval <= 1e-3){text(x=i+0.5,y=ySig,'**',cex=cexPV*1.1,col=4)
				}else if(pval <= 0.01){text(x=i+0.5,y=ySig,'*',cex=cexPV*1.1)
				}else{text(x=i+0.5,y=ySig,'NS',cex=cexPV*0.8)}
				
			}
		}
	}
}

g2anchors <- function(data,extB=1){
	data <- getAnchors(data)
	s1 <- gintervals(data$chrom,data$start,data$start+extB)
	s2 <- gintervals(data$chrom,data$end,data$end+extB)
	s <- gintervals.canonic(rbind(s1,s2))
	return(s)
}

g2b <- function(data,ext=1){
	data <- data[data$end > data$start,]
	keep <- apply(data,1,function(x){if(as.numeric(as.character(x['start']))-ext > 0 & as.numeric(as.character(x['end']))+ext < gintervals.all()[gintervals.all()$chrom==x['chrom'],'end']){return(1)}else{return(0)}})
	data <- data[which(keep == 1),]
	
	s1 <- gintervals(data$chrom,data$start-ext,data$start+ext)
	s2 <- gintervals(data$chrom,data$end-ext,data$end+ext)
	s <- gintervals.canonic(rbind(s1,s2))
	return(s)
}

plotLoopMaps <- function(loops, outFile, binSize=5e3, scoreCol=colorRampPalette(c('#428bca','#e2e3e4','#414040'))(200), hicDatasets=c('ED','LE'), scoreTracks=c('hic.ED_iced_5kb.ED_iced_5kb','hic.LE_iced_5kb.LE_iced_5kb'), extend=150e3, scale=TRUE){
	
	print(paste('Starting on',nrow(loops),'regions'))
	check <- 1:nrow(loops)
	
	geneCoordinates <- gintervals.load('intervals.ucscCanGenes')
	rownames(geneCoordinates) <- geneCoordinates$geneName
	geneCoordinates <- geneCoordinates[-grep(':|CG|CR',rownames(geneCoordinates)),]
	##############################################################################
	
	outputMaps <- list()
	
	for (i in check){
		progress(i, max.value=(max(check)-1), progress.bar = TRUE)
		chr <- as.character(loops[i,'chrom1'])
		res <- loops[i,'end1']-loops[i,'start1']
		dist <- (loops[i,'start2']-loops[i,'end1']) #+res/2
		mid <- loops[i,'end1']+dist/2
		start <- loops[i,'start1'] - extend
		end <- loops[i,'end2'] + extend
		
		if(start < 0){start <- 0}
		if(end > gintervals.all()[gintervals.all()$chrom == chr,'end']){end <- gintervals.all()[gintervals.all()$chrom == chr,'end']}
		
		mapLimits <- gintervals.2d(chr,start,end,chr,start,end)
		interval <- gintervals(chr,start,end)
		genes <- gintervals.neighbors(geneCoordinates, interval, maxdist=0)
		
		fMaps <- list()
		for (set in hicDatasets)
		{
			#### Scores
# 			tracks <- scoreTracks[grep(set,scoreTracks)]
			tracks <- scoreTracks[which(hicDatasets == set)]
			### vTrack
			vTrackScore <- paste0('v_',tracks,'_score')
			gvtrack.create(vTrackScore, tracks, 'avg')
			
			### iterator
			binnedIterator <- giterator.intervals(intervals=mapLimits, iterator=c(binSize,binSize))
			scores <- gextract(vTrackScore, binnedIterator, iterator=binnedIterator, colnames=c('score'))
			
			if(scale == TRUE){scores$score <- scaleData(scores$score,-99,99,1)}
			
# 			if(sum(is.na(scores$score))/nrow(scores) > 0.4){
# 				if(!is.null(outputMaps[[paste0('map',i)]][['pVal']])){
# 					outputMaps <- outputMaps[-length(outputMaps)]
# 				}
# 				next;
# 			}
			##############################################################
			outputMaps[[paste0('map',i)]][[set]] <- scores
			outputMaps[[paste0('map',i)]][['genes']] <- genes
			outputMaps[[paste0('map',i)]][['mapLimits']] <- mapLimits
			outputMaps[[paste0('map',i)]][['interval']] <- interval
			outputMaps[[paste0('map',i)]][['circle']] <- c(mid,dist,res)
			outputMaps[[paste0('map',i)]][['loop']] <- loops[i,]
			outputMaps[[paste0('map',i)]][['pVal']] <- loops[i,'pVal']
		}
	}
	
	print(paste('Found',length(names(outputMaps)),'regions'))
	
	allLoops <- c()
	for(l in names(outputMaps)){
		allLoops <- rbind(allLoops, outputMaps[[l]][['loop']])
	}
	
	print(allLoops)
	
	if(length(hicDatasets) == 2){
		pdf(outFile, width=32,height=10)
		layout(matrix(1:4,ncol=2), widths=c(1,1), heights=c(0.5,0.2,0.5,0.2), respect=T)
	}else{
		height <- 10
		width <- height*length(hicDatasets)
		
		pdf(outFile, width=width,height=height)
		layout(matrix(1:(length(hicDatasets)*2),ncol=length(hicDatasets), byrow=F), widths=rep(1,length(hicDatasets)), heights=c(0.5,0.2), respect=T)
		
# 		pdf(outFile, width=16,height=10)
# 		layout(matrix(1:2,ncol=1), widths=c(1,1), heights=c(0.5,0.2), respect=T)
	}
	par(mai=rep(0.3,4), xaxs="i", yaxs="i", family='sans', cex=0.8)

# 	peIterator <- c()
	for (i in names(outputMaps)){
		
		progress(which(names(outputMaps)==i), max.value=(length(names(outputMaps))-1), progress.bar = TRUE)
		
		mapLimits <- outputMaps[[i]][['mapLimits']]
		interval <- outputMaps[[i]][['interval']]
		genes <- outputMaps[[i]][['genes']]
		circle <- outputMaps[[i]][['circle']]
		loop <- outputMaps[[i]][['loop']]
		
# 		testArea <- expand2D(center2D(loop),25e3)
# 		if(is.null(testArea)){next}
		
		chr <- as.character(interval$chrom)
		start <- interval$start
		end <- interval$end
		
		for (set in hicDatasets){
			scores <- outputMaps[[i]][[set]]
			if(is.null(scores)){plot(0); plot(0); next}
			plotReadRhombus(scores,mapLimits,binSize,scoreCol,tri=TRUE,title=paste(set,chr,start/1e6,end/1e6,"\n", outputMaps[[i]][['pVal']]),axis=TRUE,plotClmn='score')
			abline(a=0,b=1, col = 'black', lty = 1, lwd = 4)
			draw.circle(circle[1],circle[2],circle[3], nv = 1000, border = 'black', col = NA, lty = 2, lwd = 3)
			plotRectangles2D(loop, mapLimits, col=4, lwd=2.3, image='TRI', side='top')
			plotRectangles2D(allLoops, mapLimits, col=1, lwd=1.3, image='TRI', side='top')
# 			plotRectangles2D(testArea, mapLimits, col=1, lwd=1.3, image='TRI', side='top')
			genesPlot(genes,c(start,end))
		}
# 		peIterator <- rbind(peIterator, loop)
	}
	dev.off()
}
colorRampPalette(c('#428bca','#e2e3e4','#414040'))(200)

plotMaps <- function(data, outFile, binSize=5e3, scoreCol=colorRampPalette(c('darkblue','lightgrey','darkred'))(200), hicDatasets=c('ED','LE'), scoreTracks=c('hic.ED_iced_2kb.ED_iced_2kb','hic.LE_iced_2kb.LE_iced_2kb'), doms=FALSE, markup=FALSE, markup2=FALSE, extend=150e3, scale=TRUE){
	
	print(paste('Starting on',nrow(data),'regions'))
	check <- 1:nrow(data)
	
	geneCoordinates <- gintervals.load('intervals.ucscCanGenes')
	rownames(geneCoordinates) <- geneCoordinates$geneName
# 	geneCoordinates <- geneCoordinates[-grep(':|CG|CR',rownames(geneCoordinates)),]
	##############################################################################
	
	outputMaps <- list()
	
	for (i in check){
		progress(i, max.value=(max(check)-1), progress.bar = TRUE)
		chr <- as.character(data[i,'chrom1'])
		res <- data[i,'end1']-data[i,'start1']
		dist <- (data[i,'start2']-data[i,'end1']) #+res/2
		mid <- data[i,'end1']+dist/2
		start <- data[i,'start1'] - extend
		end <- data[i,'end2'] + extend
		
		if(start < 0){start <- 0}
		if(end > gintervals.all()[gintervals.all()$chrom == chr,'end']){end <- gintervals.all()[gintervals.all()$chrom == chr,'end']}
		
		mapLimits <- gintervals.2d(chr,start,end,chr,start,end)
		interval <- gintervals(chr,start,end)
		genes <- gintervals.neighbors(geneCoordinates, interval, maxdist=0)
		
		fMaps <- list()
		for (set in hicDatasets)
		{
			#### Scores
# 			tracks <- scoreTracks[grep(set,scoreTracks)]
			tracks <- scoreTracks[which(hicDatasets == set)]
			### vTrack
			vTrackScore <- paste0('v_',tracks,'_score')
			gvtrack.create(vTrackScore, tracks, 'avg')
			
			### iterator
			binnedIterator <- giterator.intervals(intervals=mapLimits, iterator=c(binSize,binSize))
			scores <- gextract(vTrackScore, binnedIterator, iterator=binnedIterator, colnames=c('score'))
			
			if(scale == TRUE){scores$score <- scaleData(scores$score,-99,99,1)}
			
			if(doms[1] != FALSE){
				if(length(doms) == length(hicDatasets)){
					tads <- gintervals.load(doms[which(hicDatasets == set)])
				}else{
					tads <- gintervals.load(doms[1])
				}
				tads <- gintervals.neighbors(g2d(tads),mapLimits,maxdist1=0,maxdist2=0,mindist1=0,mindist2=0)[,1:6]
			}else{tads <- NA}
			
			if(is.data.frame(markup)){
				secInts <- gintervals.neighbors(markup,mapLimits,maxdist1=0,maxdist2=0,mindist1=0,mindist2=0)[,1:6]
			}else{secInts <- NA}
			
			if(is.data.frame(markup2)){
				secInts2 <- gintervals.neighbors(markup2,mapLimits,maxdist1=0,maxdist2=0,mindist1=0,mindist2=0)[,1:6]
			}else{secInts2 <- NA}
			
			##############################################################
			outputMaps[[paste0('map',i)]][[set]] <- scores
			outputMaps[[paste0('map',i)]][['genes']] <- genes
			outputMaps[[paste0('map',i)]][['mapLimits']] <- mapLimits
			outputMaps[[paste0('map',i)]][['interval']] <- interval
			outputMaps[[paste0('map',i)]][['tads']][[set]] <- tads
			outputMaps[[paste0('map',i)]][['markup']] <- secInts
			outputMaps[[paste0('map',i)]][['markup2']] <- secInts2
		}
	}
	
	print(paste('Plotting',length(names(outputMaps)),'regions'))
	
	allLoops <- c()
	for(l in names(outputMaps)){
		allLoops <- rbind(allLoops, outputMaps[[l]][['loop']])
	}
	
	if(length(hicDatasets) == 2){
		pdf(outFile, width=32,height=10)
		layout(matrix(1:4,ncol=2), widths=c(1,1), heights=c(0.5,0.2,0.5,0.2), respect=T)
	}else{
		height <- 10
		width <- height*length(hicDatasets)
		
		pdf(outFile, width=width,height=height)
		layout(matrix(1:(length(hicDatasets)*2),ncol=length(hicDatasets), byrow=F), widths=rep(1,length(hicDatasets)), heights=c(0.5,0.2), respect=T)
		
# 		pdf(outFile, width=16,height=10)
# 		layout(matrix(1:2,ncol=1), widths=c(1,1), heights=c(0.5,0.2), respect=T)
	}
	par(mai=rep(0.3,4), xaxs="i", yaxs="i", family='sans', cex=0.8)
	
	for (i in names(outputMaps)){
		
		progress(which(names(outputMaps)==i), max.value=(length(names(outputMaps))-1), progress.bar = TRUE)
		
		mapLimits <- outputMaps[[i]][['mapLimits']]
		interval <- outputMaps[[i]][['interval']]
		genes <- outputMaps[[i]][['genes']]
# 		tads <- outputMaps[[i]][['tads']]
		secInts <- outputMaps[[i]][['markup']]
		secInts2 <- outputMaps[[i]][['markup2']]
		
		chr <- as.character(interval$chrom)
		start <- interval$start
		end <- interval$end
		
		for (set in hicDatasets){
			scores <- outputMaps[[i]][[set]]
			tads <- outputMaps[[i]][['tads']][[set]]
			if(is.null(scores)){plot(0); plot(0); next}
			plotReadRhombus(scores,mapLimits,binSize,scoreCol,tri=TRUE,title=paste(set,chr,start/1e6,end/1e6,"\n", outputMaps[[i]][['pVal']]),axis=TRUE,plotClmn='score')
			abline(h=0, col = 'black', lty = 1, lwd = 1.5)
			
			if(is.data.frame(tads)){plotRectangles2D(tads, mapLimits, col=1, lwd=1.5, image='TRI', side='top')}
			if(is.data.frame(secInts)){plotRectangles2D(secInts, mapLimits, col='yellow', lty=1, lwd=1.5, image='TRI', side='top')}
			if(is.data.frame(secInts2)){plotRectangles2D(secInts2, mapLimits, col='black', lty=1, lwd=1.5, image='TRI', side='top')}
			
			genesPlot(genes,c(start,end))
		}
	}
	dev.off()
}

findSummit <- function(set,dataset,res=10){
	
	dataset <- gtrack.ls(dataset)
	
	if(length(dataset) > 1){
		dataset <- dataset[1]
		print(paste('Multiple datasets found, calculating for ',dataset))
	}
	
	if(is.na(dataset)){
		print('No dataset found')
		return('No dataset found')
	}
	
	data <- gextract(dataset,set,iterator=res,colnames='score')
	data[is.na(data)] <- 0
	
	retData <- ddply(data, .(intervalID), function(x){
		i <- which(x$score == max(x$score))[1]
		return(x[i,])
	})
	
	return(retData)
	
}

mergeTouching2D <- function(data){
	
	revData <- c()
	for (chr in unique(as.character(data$chrom1))){
		inData <- data[data$chrom1 == chr,1:6]
		print(chr)
		for(set in unique(paste0(inData$start2,'_',inData$end2))){
			cData <- inData[paste0(inData$start2,'_',inData$end2) %in% set,]
			cData <- cData[order(cData$start1),]
			
			s1 <- cData[1,5]
			e1 <- cData[1,6]
			s2 <- cData[1,2]
			e2 <- cData[1,3]
			
			if(nrow(cData) > 1){
				for(i in 2:nrow(cData)){
					p1 <- cData[i-1,'start1']
					p2 <- cData[i-1,'end1']
					p3 <- cData[i,'start1']
					p4 <- cData[i,'end1']
					
					if(p2 != p3){
						revData <- rbind(revData,gintervals.2d(chr,s2,p2,chr,s1,e1))
						s2 <- p3
					}
					if(i == nrow(cData)){revData <- rbind(revData,gintervals.2d(chr,s2,p4,chr,s1,e1))}
				}
			}else{
				revData <- rbind(revData,gintervals.2d(chr,s2,e2,chr,s1,e1))
			}
		}
	}
	
	forData <- c()
	for (chr in unique(as.character(revData$chrom1))){
		inData <- revData[revData$chrom1 == chr,1:6]
		print(chr)
		for(set in unique(paste0(inData$start1,'_',inData$end1))){
			cData <- inData[paste0(inData$start1,'_',inData$end1) %in% set,]
			cData <- cData[order(cData$start2),]
			
			s1 <- cData[1,2]
			e1 <- cData[1,3]
			s2 <- cData[1,5]
			e2 <- cData[1,6]
			
			if(nrow(cData) > 1){
				for(i in 2:nrow(cData)){
					p1 <- cData[i-1,'start2']
					p2 <- cData[i-1,'end2']
					p3 <- cData[i,'start2']
					p4 <- cData[i,'end2']
					
					if(p2 != p3){
						forData <- rbind(forData,gintervals.2d(chr,s1,e1,chr,s2,p2))
						s2 <- p3
					}
					if(i == nrow(cData)){forData <- rbind(forData,gintervals.2d(chr,s1,e1,chr,s2,p4))}
				}
			}else{
				forData <- rbind(forData,gintervals.2d(chr,s1,e1,chr,s2,e2))
			}
		}
	}
	
	return(forData)
# 	return(revData)
}



scaleData <- function(data,from=0,to=1,round=0.01){
	
	xM <- max(data,na.rm=T)
	xm <- min(data,na.rm=T)
	data <- sapply(data,function(x){(x-xm)/(xM-xm)})
	
	data <- round_any((data*(to-from)),round)
	data <- data + from
	
	data[data > to] <- to
	data[data < from] <- from
	
	return(data)
}

scaleMatrix <- function(data,from=0,to=1,round=0.01){
	lData <- as.vector(as.matrix(data))
	lData <- scaleData(lData,from,to,round)
	lData <- matrix(lData,ncol=ncol(data),byrow=F)
	return(lData)
}

prepareRadialTracks <- function(trackName, vType, res, scale){
	
	layers <- seq(res,scale,by=res)
	vNames <- c()
	for (shift in layers){
		
		vName <- paste0('l',which(layers == shift))
		print(vName)
		if(length(grep(vName,gvtrack.ls())) > 0){gvtrack.rm(vName)}
		gvtrack.create(vName, trackName, vType)
		
		gvtrack.iterator.2d(vName, 
		eshift1=shift, sshift1=-shift, 
		eshift2=shift, sshift2=-shift)
		
		vNames <- append(vNames, vName)
	}
	return(vNames)
}

plotVariableRhombus <- function(track, interval, colors=colorRampPalette(c('lightgrey','darkred'))(100), inverse=FALSE, mirror=FALSE, add=FALSE, tri=FALSE, axis=FALSE, band='false', plotClmn=7, title=''){
	
	if(plotClmn != 7){plotClmn <- which(colnames(track) == plotClmn)}
	
	if(tri == FALSE){
		### Squares
		track <- track[track$start2 > track$start1,]
		xAxis <- c(interval$start1,interval$end1)
		yAxis <- c(interval$start2,interval$end2)
		if(add == FALSE){
			plot(0,xlim=xAxis,ylim=yAxis, main=title, col='white',cex.main=2,xlab='', ylab='', yaxt='n', xaxt='n', frame=F)
		}
		
		apply(track, 1, function(x)
		{
			xx <- c(as.numeric(as.character(x['start1'])),as.numeric(as.character(x['start1'])),as.numeric(as.character(x['end1'])),as.numeric(as.character(x['end1'])))
			yy <- c(as.numeric(as.character(x['start2'])),as.numeric(as.character(x['end2'])),as.numeric(as.character(x['end2'])),as.numeric(as.character(x['start2'])))
			col <- colors[as.numeric(as.character(x[plotClmn]))+101]
			if(add == FALSE){
				polygon(xx,yy, lwd=1, col=col,border=col)
			}
			
			if(mirror == TRUE | inverse == TRUE){
				polygon(yy,xx, lwd=1, col=col,border=col)
			}
			
		})
		
		if(axis==TRUE){axis(1, cex=1.6); axis(2, cex=1.6);}
		
	}
# 	else{
# 		### Triangle
# 		track <- track[track$start2 >= track$start1,]
# 		xAxis <- c(interval$start1,interval$end2)
# 		winSize <- xAxis[2]-xAxis[1]
# 		yAxis <- c(0,winSize)
# 		
# 		if(band!='false'){yAxis[2]=yAxis[2]*band}
# 		
# 		
# 		track$dist <- abs(track$start2 - track$start1)
# 		track <- track[track$dist <= yAxis[2],]
# 		
# 		plot(0,xlim=xAxis,ylim=yAxis, col='white',cex.main=2,xlab='', ylab='', yaxt='n', xaxt='n', frame=F)
# # 		rect(xAxis[1],yAxis[1],xAxis[2],yAxis[2],col='lightgrey')
# # 		points(x=seq(xAxis[1],xAxis[2], length.out=length(colors)),y=rep(yAxis[1]-binSize*2,length(colors)),type='p', cex=3, col=colors[1:length(colors)], pch=19)
# 		apply(track, 1, function(x)
# 		{
# 			start <- as.numeric(as.character(x['start1']))+binSize/2
# 			end <- as.numeric(as.character(x['start2']))+binSize/2
# 			dist <- end-start
# 			
# 			mid <- min(c(start,end))+abs(dist)/2
# 			xx <- c(mid-binSize/2,mid,mid+binSize/2,mid)
# 			height <- abs(dist)
# 			yy <- c(height,height+binSize,height,height-binSize)
# 			
# 			if(yy[4] < 0){yy[4] <- 0}
# 			if(yy[2] > yAxis[2]+binSize){yy[2] <- yAxis[2]}
# 			col <- colors[as.numeric(as.character(x[plotClmn]))+101]
# 			
# 			polygon(xx,yy, lwd=1, col=col,border=col)
# 		})
# # 		if(title != ''){legend("topleft", legend=title, cex=1.6)}
# # 		if(title != ''){text(x=xAxis[1]+(winSize*0.025), y=yAxis[2]*0.8, labels=title, cex=1.6, pos=4)}
# 		if(title != ''){text(x=xAxis[1], y=yAxis[2]*0.8, labels=title, cex=1.6, pos=4)}
# 		if(axis==TRUE){axis(1, cex=1.2);}
# 	}
}



extractChIPData <- function(chipTarget,coords,extend=0,vType='avg'){
	
	chipTrack <- gtrack.ls(chipTarget)
	if(length(chipTrack) == 0){
		return('ChIP target not found in database...')
	}else if(length(chipTrack) > 1){
		return('Multiple ChIP targets found in database...')
	}else{
		gvtrack.create('vChIP',chipTrack,vType)
	}
	
	if(extend != 0){
		gvtrack.iterator("vChIP", sshift=-abs(extend), eshift=abs(extend))
	}
	
	chipData <- gextract("vChIP",coords,iterator=coords, colnames=chipTarget)
	return(chipData)
}

getAnchors <- function(set){
	
	s1 <- gintervals.canonic(gintervals(set$chrom1,set$start1,set$end1), unify_touching_intervals=FALSE)
	s2 <- gintervals.canonic(gintervals(set$chrom2,set$start2,set$end2), unify_touching_intervals=FALSE)
	anchorSites <- rbind(s1,s2)
	anchorSites <- gintervals.canonic(anchorSites, unify_touching_intervals=FALSE)
	return(anchorSites)
	
}


prepareMaps1D <- function(data, hicDatasets, scoreTracks, binSize=5e3, extend=0, center=FALSE, chipTracks=NULL, chipRes=50, chipVirt='max'){
	outputMaps <- list()
	
	if(!is.null(chipTracks)){
		vChipTracks <- paste0('v_',chipTracks)
		for(t in 1:length(chipTracks)){gvtrack.create(vChipTracks[t],chipTracks[t],chipVirt)}
	}
	
	for (i in 1:nrow(data)){
		print(i)
		chr <- as.character(data[i,'chrom'])
		
		if(center == TRUE){
			mid <- data[i,'start']+(data[i,'end']-data[i,'start'])/2
			start <- mid - extend
			end <- mid + extend
		}else{
			
			start <- data[i,'start'] - extend
			end <- data[i,'end'] + extend
		}
		
		if(start < 0){start <- 0}
		if(end > gintervals.all()[gintervals.all()$chrom == chr,'end']){end <- gintervals.all()[gintervals.all()$chrom == chr,'end']}
		
		mapLimits <- gintervals.2d(chr,start,end,chr,start,end)
		interval <- gintervals(chr,start,end)
		
		outputMaps[[paste0('map',i)]][['mapLimits']] <- mapLimits
		outputMaps[[paste0('map',i)]][['interval']] <- interval
		
		if(!is.null(chipTracks)){
			outputMaps[[paste0('map',i)]][['chipData']] <- gextract(vChipTracks,interval,iterator=chipRes)
		}else{
			outputMaps[[paste0('map',i)]][['chipData']] <- NULL
		}
		
		for (set in hicDatasets)
		{
			#### Scores
			tracks <- scoreTracks[grep(set,scoreTracks)]
			### vTrack
			vTrackScore <- paste0('v_',tracks,'_score')
			gvtrack.create(vTrackScore, tracks, 'avg')
			### iterator
			binnedIterator <- giterator.intervals(intervals=mapLimits, iterator=c(binSize,binSize))
			scores <- gextract(vTrackScore, binnedIterator, iterator=binnedIterator, colnames=c('score'))
			outputMaps[[paste0('map',i)]][[set]] <- scores
		}
	}
	return(outputMaps)
}

colorDomains <- function(targetDomains,somDomains){
	
	if(!'id' %in% colnames(targetDomains)){targetDomains$id <- paste0('id',1:nrow(targetDomains))}
	
	targetDomains$s1 <- targetDomains$end-targetDomains$start
	somDomains$s2 <- somDomains$end-somDomains$start
	
	colorTADs <- gintervals.neighbors(targetDomains, somDomains, maxdist=0, mindist=0, maxneighbors=5000)
	
# 	print(head(colorTADs))
	
	colorTADs$d1 <- colorTADs[,ncol(targetDomains)+2]-colorTADs[,2]
	colorTADs$d2 <- colorTADs[,ncol(targetDomains)+3]-colorTADs[,3]
	
	for (i in 1:nrow(colorTADs)){
		if(colorTADs[i,'d1'] < 0){ colorTADs[i,'s2'] <- colorTADs[i,'s2'] + colorTADs[i,'d1'] }
		if(colorTADs[i,'d2'] > 0){ colorTADs[i,'s2'] <- colorTADs[i,'s2'] - colorTADs[i,'d2'] }
	}
	
	colorTADs$id <- as.character(colorTADs$id)
	
	selectColor <- unlist(sapply(unique(colorTADs$id), function(x){
		data <- colorTADs[colorTADs$id == x,]
		overlaps <- sapply(unique(as.character(somDomains$cluster)),function(y){
			sum( data[data$cluster == y,'s2'] )/data[1,'s1']
		})
# 		return(names(overlaps)[which(overlaps == max(overlaps))[1]])
                return(c(names(overlaps)[which(overlaps == max(overlaps))[1]],max(overlaps)))
	}))
	
# 	print(head(selectColor))
# 	print(class(selectColor))
	
	colorTADs$color <- selectColor[1,as.character(colorTADs$id)]
	colorTADs$ratio <- as.numeric(as.character(selectColor[2,as.character(colorTADs$id)]))
	
# 	print(head(colorTADs))
	
	colorTADs <- colorTADs[!duplicated(colorTADs$id),c(1:3,which(colnames(colorTADs) == 'id'),which(colnames(colorTADs) == 'color'),which(colnames(colorTADs) == 'ratio'))]
	colorTADs$s1 <- colorTADs[,3]-colorTADs[,2]
	
# 	print(dim(colorTADs))
	
	return(colorTADs)
}

g2d <- function(intrv){
	return(gintervals.2d(intrv$chrom,intrv$start,intrv$end,intrv$chrom,intrv$start,intrv$end))
}

g1d <- function(intrv,ut=FALSE,method=c('both','left','right')[1]){
	s1 <- intrv[,1:3]
	colnames(s1) <- c('chrom','start','end')
	s2 <- intrv[,4:6]
	colnames(s2) <- c('chrom','start','end')
	if(method == 'both'){
		regions <- rbind(s1,s2)
	}else if(method == 'left'){
		regions <- s1
	}else if(method == 'right'){
		regions <- s2
	}
	
	regions <- gintervals.canonic(gintervals(regions$chrom,regions$start,regions$end), unify_touching_intervals=ut)
	return(regions)
}

exactNeighbors <- function(set1,set2,merge=FALSE){
	ids1 <- apply(set1,1,function(x)paste(as.character(x[1:3]),collapse='_'))
	ids1 <- gsub(' ','',ids1)
	ids2 <- apply(set2,1,function(x)paste(as.character(x[1:3]),collapse='_'))
	ids2 <- gsub(' ','',ids2)
	
	if(merge == TRUE){
		res <- gintervals.neighbors(set1,set2,maxneighbors=1)
		res <- res[which(ids1 %in% ids2),]
	}else{
		res <- set1[which(ids1 %in% ids2),]
	}
	return(res)
}

exactNeighbors2D <- function(set1,set2,merge=FALSE){
	ids1 <- apply(set1,1,function(x)paste(as.character(x[1:6]),collapse='_'))
	ids1 <- gsub(' ','',ids1)
	ids2 <- apply(set2,1,function(x)paste(as.character(x[1:6]),collapse='_'))
	ids2 <- gsub(' ','',ids2)
	
	if(merge == TRUE){
		res <- gintervals.neighbors(set1,set2,maxneighbors=1)
		res <- res[which(ids1 %in% ids2),]
	}else{
		res <- set1[which(ids1 %in% ids2),]
	}
	return(res)
}

mappaCheck <- function(set, thr=1, shift=5e3,mapTrack='redb.GATC_map'){
	
	
	gvtrack.create("vMap", mapTrack, "min")
	gvtrack.iterator("vMap", sshift=-shift, eshift=shift)
	data <- gextract('vMap',set,iterator=set,colnames='vMap')
	data[is.na(data)] <- 0
	data <- data[data$vMap > thr,]
# 	return(data[,1:ncol(set)])
	return(set[data$intervalID,])
}

clustersPlot <- function(set1, plotLim, cols=c('midnightblue','seagreen','lightgrey'), vertical=FALSE, main=''){
	
	if(vertical == TRUE)
	{
		xlim=c(-50,50)
		plot(0, xlab='', ylab='', yaxt='n', xaxt='n', frame=F, ylim=plotLim, xlim=xlim, col='white', main=main)
		segments(0,plotLim[1],0,plotLim[2],lwd=1.75, col='darkgrey')
		
		for(g in 1:nrow(set1))
		{
			row <- set1[g,]
			if(is.null(row$Cluster)){
				col <- cols[1]
			}else{
				col <- cols[row$Cluster] #+ 1
			}
			rect(-50,row$start, 50, row$end, col=col, border=col)
		}
	}else{
		yLim=c(-50,50)
		plot(0, xlim=plotLim, ylim=yLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F, col='white', main=main)
		segments(plotLim[1],0,plotLim[2],0,lwd=1.75, col='darkgrey')
		
		for(g in 1:nrow(set1))
		{
			row <- set1[g,]
			if(is.null(row$Cluster)){
				col <- cols[1]
			}else{
				col <- cols[row$Cluster] #+ 1
			}
			rect(row$start, -50, row$end, 50, col=col, border=col)
		}
	}
}

genCRF_Track <- function(tracks, trackName, vType='area', scale=25e3, res=2e3, minPairs=500, export=TRUE){
	
	iter_1d <- giterator.intervals(intervals=ALLGENOME[[1]], iterator=res)
	iter_2d = gintervals.2d(chroms1 = iter_1d$chrom, starts1=iter_1d$start, ends1=iter_1d$end, chroms2 = iter_1d$chrom, starts2=iter_1d$start, ends2=iter_1d$end)
	
	if(!is.null(tracks)){
		print('preparing tracks...')
		resData <- extractData(tracks, vType, scale, iter_2d)
	}else{
		resData <- gextract('obs_center','obs_upstream','obs_downstream', iter_2d, iterator=iter_2d, band=c(-scale,-1e3))
	}
	
	resData[is.na(resData)] <- 1
	
	#################
	corrF <- 1
	UDratio <- 2
	
	resData$obs_sides <- (resData$obs_upstream+resData$obs_downstream)/corrF
	resData$UDratio <- abs(log2((resData$obs_upstream+1)/(resData$obs_downstream+1)))
	
	resData$fc <- -log2(resData$obs_center/resData$obs_sides)
	
	resData[resData$UDratio >= UDratio | resData$obs_center <= minPairs | resData$obs_upstream <= minPairs | resData$obs_downstream <= minPairs,'fc'] <- NA
	
	
	if(min(resData$fc, na.rm=TRUE) < 0){
		resData$fc <- resData$fc + abs(min(resData$fc, na.rm=TRUE))
	}else{
		resData$fc <- resData$fc - abs(min(resData$fc, na.rm=TRUE))
	}
	
	#################
	if(export == TRUE){
		intrvs <- gintervals(resData[,1],resData[,2],resData[,6])
		gtrack.create_sparse(track=trackName, intrvs, resData$fc, description=paste0('Res: ',res,',Scale: ',scale))
# 		gtrack.create(track=trackName, intrvs, resData$fc, description=paste0('Res: ',res,',Scale: ',scale))
	}
	
	return(resData)
}


genProxInsTrack <- function(tracks, trackName, vType='area', scale=25e3, res=2e3, minPairs=500, export=TRUE){
	
	iter_1d <- giterator.intervals(intervals=ALLGENOME[[1]], iterator=res)
	iter_2d = gintervals.2d(chroms1 = iter_1d$chrom, starts1=iter_1d$start, ends1=iter_1d$end, chroms2 = iter_1d$chrom, starts2=iter_1d$start, ends2=iter_1d$end)
	
	if(!is.null(tracks)){
		print('preparing tracks...')
		resData <- extractData(tracks, vType, scale, iter_2d)
	}else{
		resData <- gextract('obs_center','obs_upstream','obs_downstream', iter_2d, iterator=iter_2d, band=c(-scale,-1e3))
	}
	
	resData[is.na(resData)] <- 1
	
	#################
	corrF <- 1
	UDratio <- 2
	
	resData$obs_sides <- (resData$obs_upstream+resData$obs_downstream)/corrF
	resData$UDratio <- abs(log2((resData$obs_upstream+1)/(resData$obs_downstream+1)))
	
	resData$fc <- -log2(resData$obs_center/resData$obs_sides)
	
	resData[resData$UDratio >= UDratio | resData$obs_center <= minPairs | resData$obs_upstream <= minPairs | resData$obs_downstream <= minPairs,'fc'] <- NA
	
	
	if(min(resData$fc, na.rm=TRUE) < 0){
		resData$fc <- resData$fc + abs(min(resData$fc, na.rm=TRUE))
	}else{
		resData$fc <- resData$fc - abs(min(resData$fc, na.rm=TRUE))
	}
	
	#################
	if(export == TRUE){
		intrvs <- gintervals(resData[,1],resData[,2],resData[,6])
		gtrack.create_sparse(track=trackName, intrvs, resData$fc, description=paste0('Res: ',res,',Scale: ',scale))
# 		gtrack.create(track=trackName, intrvs, resData$fc, description=paste0('Res: ',res,',Scale: ',scale))
	}
	
	return(resData)
}

computeProximalINS <- function(intrvs, tracks=NULL, vType='area', scale=25e3, res=2e3, minPairs=500){
	
# 	iter_1d = intrvs
	iter_1d = expand1D(center1D(intrvs),res/2)
# 	expand1D(center1D(intrvs),res)
	iter_2d = gintervals.2d(chroms1 = iter_1d$chrom, starts1=iter_1d$start, ends1=iter_1d$end, chroms2 = iter_1d$chrom, starts2=iter_1d$start, ends2=iter_1d$end)
	
	if(!is.null(tracks)){
		print('multi...')
		resData <- extractData(tracks, vType, scale, iter_2d)
	}else{
		resData <- gextract('obs_center','obs_upstream','obs_downstream', iter_2d, iterator=iter_2d, band=c(-scale,-1e3))
	}
	
	resData[is.na(resData)] <- 1
	
	#################
	corrF <- 2
	UDratio <- 1
	
	
	resData$obs_sides <- (resData$obs_upstream+resData$obs_downstream)/corrF
	resData$UDratio <- abs(log2((resData$obs_upstream+1)/(resData$obs_downstream+1)))
	
	resData$fc <- -log2(resData$obs_center/resData$obs_sides)
	
# 	resData <- resData[resData$fc > minFC,]
	
	resData <- resData[resData$UDratio < UDratio & resData$obs_center > minPairs & resData$obs_upstream > minPairs & resData$obs_downstream > minPairs & resData$fc > 0,]
	
	if(!is.null(tracks)){
		for (i in 1:length(tracks)){
			resData[,paste0('fc_',i)] <- -log2(resData[,paste0('oc_',i)]/resData[,paste0('os_',i)])
		}
	}
	
	print(dim(resData))
# 	if(min(resData$fc, na.rm=TRUE) < 0){
# 		resData$fc <- resData$fc + abs(min(resData$fc, na.rm=TRUE))
# 	}else{
# 		resData$fc <- resData$fc - abs(min(resData$fc, na.rm=TRUE))
# 	}
	
	return(resData)
}

prepareTracks <- function(trackName, vType, scale){
# 	print('Preparing V Tracks')
	if(length(grep('obs_center',gvtrack.ls())) > 0){gvtrack.rm('obs_center')}
	gvtrack.create("obs_center", trackName, vType)
	if(length(grep('obs_upstream',gvtrack.ls())) > 0){gvtrack.rm('obs_upstream')}
	gvtrack.create("obs_upstream", trackName, vType)
	if(length(grep('obs_downstream',gvtrack.ls())) > 0){gvtrack.rm('obs_downstream')}
	gvtrack.create("obs_downstream", trackName, vType)
	
	gvtrack.iterator.2d("obs_center", 
		eshift1=0, sshift1=-scale, 
		eshift2=scale, sshift2=0)
	
	gvtrack.iterator.2d("obs_upstream", 
		eshift1=0, sshift1=-scale, 
		eshift2=0, sshift2=-scale)
	
	gvtrack.iterator.2d("obs_downstream", 
		eshift1=scale, sshift1=0, 
		eshift2=scale, sshift2=0)
}

extractData <- function(tracks, vType, scale, iterator){
	
	prepareTracks(tracks[1],vType,scale)
	retData <- gextract('obs_center','obs_upstream','obs_downstream', iterator, iterator=iterator, band=c(-scale,-1e3))
	
	retData$oc_1 <- retData$obs_center
	retData$ou_1 <- retData$obs_upstream
	retData$od_1 <- retData$obs_downstream
	retData$os_1 <- (retData$obs_upstream+retData$obs_downstream)/2
	
	for (track in tracks[-1]){
		prepareTracks(track,vType)
		
		resData <- gextract('obs_center','obs_upstream','obs_downstream', iterator, iterator=iterator, band=c(-scale,-1e3))
		retData[,paste0('oc_',which(tracks == track))] <- resData$obs_center
		retData[,paste0('ou_',which(tracks == track))] <- resData$obs_upstream
		retData[,paste0('od_',which(tracks == track))] <- resData$obs_downstream
		retData[,paste0('os_',which(tracks == track))] <- (resData$obs_upstream+resData$obs_downstream)/2
		
		retData$obs_center = rowSums(cbind(retData$obs_center, resData$obs_center), na.rm=T)
		retData$obs_upstream = rowSums(cbind(retData$obs_upstream, resData$obs_upstream), na.rm=T)
		retData$obs_downstream = rowSums(cbind(retData$obs_downstream, resData$obs_downstream), na.rm=T)
	}
	return(retData)
}


plotFC <- function(intrvs, tracks=NULL, vType='area', scale=25e3, res=2e3, main='test', minFC=0){
	
	iter_1d = expand1D(center1D(intrvs),res/2)
# 	iter_1d = intrvs
	iter_2d = gintervals.2d(chroms1 = iter_1d$chrom, starts1=iter_1d$start, ends1=iter_1d$end, chroms2 = iter_1d$chrom, starts2=iter_1d$start, ends2=iter_1d$end)
	
	if(!is.null(tracks)){
		print('multi...')
		resData <- extractData(tracks, vType, scale, iter_2d)
	}else{
# 		resData <- gextract('obs_center','obs_sides', iter_2d, iterator=iter_2d, band=c(-scale,-1e3))
		resData <- gextract('obs_center','obs_upstream','obs_downstream', iter_2d, iterator=iter_2d, band=c(-scale,-1e3))
	}
	
	resData[is.na(resData)] <- 1
	
	#################
	minPairs <- 500
	resData <- resData[resData$obs_center > minPairs,]
	resData <- resData[resData$obs_upstream > minPairs,]
	resData <- resData[resData$obs_downstream > minPairs,]
	upDownRatio <- 1
	resData <- resData[abs(log2((resData$obs_upstream+1)/(resData$obs_downstream+1))) < upDownRatio,]
	corrF <- 2
	resData$obs_sides <- (resData$obs_upstream+resData$obs_downstream)/corrF
	
	resData$fc <- -log2((resData$obs_center+1)/(resData$obs_sides+1))
	allReads <- sum(sum(resData$obs_center,na.rm=T),sum(resData$obs_sides,na.rm=T))
	
# 	resData$fc <- -log2((resData$obs_center/allReads)/(resData$obs_sides/allReads))
	
# 	print(head(resData))
	print(allReads)
	
	resData <- resData[resData$fc > minFC,]
	
	all <- -log2(sum(resData$obs_center,na.rm=T)/sum(resData$obs_sides,na.rm=T))
	
	resData <- resData[order(resData$fc),]
	
	plot(resData$fc, pch=19, cex=0.5, main=main, ylim=c(0,4), frame=F)
	
	if(!is.null(tracks)){
		
		for (i in 1:length(tracks)){
			
# 			resData[,paste0('os_',i)] <- (resData[,paste0('os_',i)]-resData[,paste0('oc_',i)])/corrF
			fc <- -log2(resData[,paste0('oc_',i)]/resData[,paste0('os_',i)])
			resData[,paste0('fc_',i)] <- fc
			points(fc, col=i+1, pch=19, cex=0.4)
		}
	}
	legend('topleft',legend=c(signif(all,4),signif(allReads/1e6,4), nrow(resData)))
	
	return(resData)
	
}

meanNorm <- function(in_data)
{
	norm_dataset <- as.matrix(in_data)
	for(i in 1:ncol(norm_dataset))
	{
		current_set <- subset(norm_dataset, norm_dataset[,i] != 0)
		mean_value <- mean(current_set[,i], na.rm=T)
		for (j in 1:nrow(norm_dataset))
		{
			if (norm_dataset[j,i] != 0)
			{
				######## MEAN NORM
				norm_dataset[j,i] = norm_dataset[j,i]/mean_value
			}
		}
	}
	return(norm_dataset)
}

makeDomains <- function(intervals, shiftStart=0){
	
	colnames(intervals)[1:3] <- c('chrom','start','end')
	intervals <- gintervals(intervals$chrom, intervals$start, intervals$end)
	retDomains <- gintervals(intervals[1,'chrom'],intervals[1,'start']+shiftStart,intervals[2,'start']+shiftStart)
	for (i in 2:(nrow(intervals)-1)){
		if(intervals[i,'chrom']==intervals[i+1,'chrom']){
			retDomains <- rbind(retDomains, gintervals(intervals[i,'chrom'],intervals[i,'start']+shiftStart,intervals[i+1,'start']+shiftStart))
		}
	}
	retDomains$id <- paste0('tadID',1:nrow(retDomains))
	return(retDomains)
}

clusterHeatmap <- function(data,kClusters,method='kmeans',cl=NA,colors=c('darkblue','lightgrey','darkred'),title='',outfile='heatmapOut.pdf',size=c(3,7), clColumn=FALSE, ret=FALSE, showRows=FALSE){
	
	require(pheatmap)
	set.seed(2)
	
	data[is.na(data)] <- 0
	if(method != 'none'){
		clusters <- kmeans(data, kClusters)
		kData <- cbind(data, clusters$cluster[match(rownames(data), names(clusters$cluster))])
		colnames(kData)[ncol(kData)] <- 'Cluster'
		kData <- as.data.frame(kData)
		rowOrder <- order(kData$Cluster, rowSums(kData[,colnames(data)]), decreasing=T)
		kData <- kData[rowOrder,]
		
		annotation_row <- data.frame(Cluster = kData$Cluster)
		rownames(annotation_row) <- rownames(kData)
		aCol <- colorRampPalette(c('grey85','darkorange'))(kClusters)
		ann_col <- list(Cluster=aCol[unique(kData$Cluster)])
		
		gaps <- c()
		for (k in 1:kClusters){gaps <- append(gaps, max(which(kData$Cluster == k)))}
		
		data <- data[rowOrder,]
		
		breaks <- seq(min(data),max(data),length.out=200)
		cC <- colorRampPalette(colors)(length(breaks)-1)
		
		
		dev.new()
		pheatmap(data, scale="none", breaks=breaks, color=cC, main=title, show_rownames=showRows, cluster_cols=clColumn, cluster_rows=F , gaps_row=gaps, 
		annotation_row=annotation_row,annotation_colors=ann_col, annotation_legend = FALSE, width=size[1],height=size[2], filename=outfile, border_color=NA)
# 		dev.off()
		
		
	}else{
		breaks <- seq(min(data),max(data),length.out=200)
		cC <- colorRampPalette(colors)(length(breaks)-1)
		
		annotation_row <- data.frame(Cluster = cl$order)
		rownames(annotation_row) <- rownames(data)
		aCol <- colorRampPalette(c('grey85','darkorange'))(kClusters)
		ann_col <- list(Cluster=aCol[unique(cl$clusters$cluster[cl$order])])
		
		gaps <- c()
		for (k in 1:kClusters){gaps <- append(gaps, max(which(cl$clusters$cluster[cl$order] == k)))}
		
# 		pheatmap(data, scale="none", breaks=breaks, color=cC, main=title, show_rownames=F, cluster_cols=clColumn, cluster_rows=F , gaps_row=gaps, 
# 		annotation_row=annotation_row,annotation_colors=ann_col, annotation_legend = FALSE, width=size[1],height=size[2], filename=outfile, border_color=NA)
		
		pheatmap(data, scale="none", breaks=breaks, color=cC, main=title, show_rownames=F, cluster_cols=clColumn, cluster_rows=F , gaps_row=gaps, 
		annotation_legend = FALSE, width=size[1],height=size[2], filename=outfile, border_color=NA)
		
# 		pheatmap(data, scale="none", breaks=breaks, color=cC, main=title, show_rownames=F, cluster_cols=clColumn, cluster_rows=F , 
# 		annotation_legend = FALSE, width=size[1],height=size[2], border_color=cC, filename=outfile)
	}
	
	
	
	
	if(ret==TRUE){
		return(list(order=rowOrder,clusters=clusters))
	}
}

clusterMotifs <- function(data,kClusters,cl=NA,colors=c('darkblue','lightgrey','darkred'),colorCuts=200,title='',outfile='heatmapOut.pdf',size=c(3,7), clColumn=FALSE, ret=FALSE, showRows=FALSE){
	
	require(pheatmap)
	set.seed(2)
	
	data[is.na(data)] <- 0
	clusters <- kmeans(data, kClusters)
	
	#### ORDER BASED ON MOTIF POSITION
	mData <- data
	mData[mData != 0] <- 1
	
	mData <- cbind(mData, clusters$cluster[match(rownames(mData), names(clusters$cluster))])
	colnames(mData) <- paste0('c',1:ncol(mData))
	colnames(mData)[ncol(mData)] <- 'k'
	mData <- as.data.frame(mData)
	mData <- mData[order(mData$k),]
	
	clOrder <- c()
	for (k in 1:max(mData$k)){
		cData <- mData[mData$k == k,1:(ncol(mData)-1)]
		clOrder <- append(clOrder, which(colSums(cData) == max(colSums(cData)))[1])
	}
	mData$l <- sapply(mData$k,function(x)which(order(clOrder) == x))
	
	mData$k <- mData$l
	mData <- mData[,-ncol(mData)]
	
	kData <- cbind(data, mData[match(rownames(data), rownames(mData)),'k'])
	colnames(kData)[ncol(kData)] <- 'Cluster'
	kData <- as.data.frame(kData)
	rowOrder <- order(kData$Cluster, rowSums(kData[,colnames(data)]), decreasing=T)
	kData <- kData[rowOrder,]
	
	annotation_row <- data.frame(Cluster = kData$Cluster)
	rownames(annotation_row) <- rownames(kData)
	aCol <- colorRampPalette(c('grey85','darkorange'))(kClusters)
	ann_col <- list(Cluster=aCol[unique(kData$Cluster)])
	
	gaps <- c()
	for (k in 1:kClusters){gaps <- append(gaps, max(which(kData$Cluster == k)))}
	
	data <- data[rowOrder,]
	
	breaks <- seq(min(data),max(data),length.out=colorCuts)
	cC <- colorRampPalette(colors)(length(breaks)-1)
	
	
	dev.new()
	pheatmap(data, scale="none", breaks=breaks, color=cC, main=title, show_rownames=showRows, cluster_cols=clColumn, cluster_rows=F , gaps_row=gaps, 
	annotation_row=annotation_row,annotation_colors=ann_col, annotation_legend = FALSE, width=size[1],height=size[2], filename=outfile, border_color=NA)
# 	dev.off()
	
	if(ret==TRUE){
		return(list(order=rowOrder,data=kData))
	}
}

plotChIP <- function(chipData, plotOrder=NA, plotCols=NA, vertical=FALSE, chipRes=10, plotAnn=F, plotStart=F, inputType='avg', cex=1, plotMax=NA, main=TRUE, plotAxis=TRUE){
	
	cNames <- colnames(chipData)[-grep('chrom|start|end|intervalID', colnames(chipData))]
	cNames <- gsub('chipseq.|v_chipseq.','',cNames)
	if(!is.na(plotOrder)){
		cNames <- cNames[plotOrder]
	}
	
	for (chip in cNames)
	{
		if (length(grep(chip,colnames(chipData))) == 0){next}
		
		if (length(grep(chip,colnames(chipData))) > 1){
			allChIP <- colnames(chipData)[grep(chip,colnames(chipData))]
			chip <- allChIP[which(nchar(allChIP) == min(nchar(allChIP)))]
			currentData <- chipData[,which(colnames(chipData) == chip)]
			chip <- gsub('chipseq.|v_chipseq.','',chip)
		}else{
			currentData <- chipData[,grep(chip, colnames(chipData))]
		}
		
		if(is.na(plotCols[1])){col='darkslategrey'
		}else{
			col <- plotCols[which(chip == cNames)]
		}
		
		currentData[is.na(currentData)] <- 0
		if(inputType == 'log'){currentData <- log2(currentData+1)}
		if(inputType == 'gpm'){currentData <- -log2(1-currentData)}
		
		if(length(grep('.ins_',chip)) > 0){currentData <- -currentData}
		
		if(main==TRUE){
# 			title <- gsub('chipseq.|v_chipseq.|v_bb_chip.|v_hic.|_KCL|_R1|_R2','',chip)
			title <- gsub('chipseq.|v_chipseq.|v_bb_chip.|v_hic.|_KCL','',chip)
		}else{
			title <- ''
		}
		
		if(vertical == TRUE){
			barplot(currentData, cex.main=1.6, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, main=gsub('v_chipseq.|v_bb_chip|v_hic.','',chip), horiz=T, space=0)
			axis(1, cex.axis=cex, las=1)
			
			if(is.data.frame(plotAnn))
			{
				for(g in 1:nrow(plotAnn))
				{
					row <- plotAnn[g,]
					rect(0,(row$start-plotStart)/chipRes, 1000, (row$end-plotStart)/chipRes, col=rgb(0.66,0.66,0.66, alpha=0.45), border = 0)
				}
			}
		}else{
			if(!is.na(plotMax[1])){
				yLim <- c(0,plotMax[which(cNames == chip)])
				barplot(currentData, cex.main=cex*1.2, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, main=title, space=0, ylim=yLim)
				
				if(plotAxis == TRUE){
					axis(2, at=yLim, labels=signif(yLim,2), cex.axis=cex, las=1)
				}
			}else{
				barplot(currentData, cex.main=cex*1.2, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, col=col, border=col, main=title, space=0)
				if(plotAxis == TRUE){
					axis(2, at=c(0,max(currentData,na.rm=TRUE)), labels=c(0,signif(max(currentData,na.rm=TRUE),2)), cex.axis=cex, las=3)
				}
			}
# 		axis(2, cex.axis=cex, las=2)
		
		if(is.data.frame(plotAnn)){
# 			print('plotting annotation')
			
			for(g in 1:nrow(plotAnn))
			{
				row <- plotAnn[g,]
# 				rect((row$start-plotStart)/chipRes, 0, (row$end-plotStart)/chipRes, 1000, col=rgb(0.66,0.66,0.66, alpha=0.5), border = 0)
# 				rect((row$start-plotStart)/chipRes, 0, (row$end-plotStart)/chipRes, 1000, col=rgb(0,0.33,0.66, alpha=0.2), border = 0)
				rect((row$start-plotStart)/chipRes, 0, (row$end-plotStart)/chipRes, 1000, col=rgb(0.83,0.83,0.83, alpha=0.2), border = 0)
# 				abline(v=row$start,col=4)
# 				abline(h=row$start,col=4)
			}
		}
		}
	}
}


sampleMatchedPeaks <- function(peaks1, peaks2, eMargin=0.1, keepChr=FALSE){
	nonMatched <- 0
	set.seed(13)
	if(keepChr == TRUE){
		rPeaks <- c()
		for (chr in unique(peaks1$chrom)){
			selSet <- peaks1[peaks1$chrom == chr,]
			fullSet <- peaks2[peaks2$chrom == chr,]
			fullSet$peakID <- paste(fullSet$chrom,fullSet$start,fullSet$end, sep='_')
			
			for (i in 1:(nrow(selSet)))
			{
				eLevel <- as.numeric(as.character(selSet[i,'level']))
				
				### Select set of suitable peaks
				availablePeaks <- fullSet[fullSet$level >= eLevel-eMargin*eLevel & fullSet$level <= eLevel+eMargin*eLevel, ]
				
				if(nrow(availablePeaks) == 0)
				{
					availablePeaks <- fullSet
					nonMatched <- nonMatched + 1
					rPeaks <- rbind(rPeaks, fullSet[sample(nrow(fullSet), 1),])
				}else{
					rPeaks <- rbind(rPeaks, availablePeaks[sample(nrow(availablePeaks), 1),])
				}
				
				### Remove selected genes from pool
				fullSet <- fullSet[-which(fullSet$peakID %in% rPeaks$peakID),]
			}
		}
		return(list(peaks=rPeaks,nonMatched=nonMatched))
		
	}else{
		selSet <- peaks1
		fullSet <- peaks2
		fullSet$peakID <- paste(fullSet$chrom,fullSet$start,fullSet$end, sep='_')
		
		rPeaks <- c()
		
		for (i in 1:(nrow(selSet)))
		{
			eLevel <- as.numeric(as.character(selSet[i,'level']))
			
			### Select set of suitable peaks
			availablePeaks <- fullSet[fullSet$level >= eLevel-eMargin*eLevel & fullSet$level <= eLevel+eMargin*eLevel, ]
			if(nrow(availablePeaks) == 0)
			{
				rPeaks <- rbind(rPeaks, fullSet[sample(nrow(fullSet), 1),])
				nonMatched <- nonMatched + 1
				
			}else{
				rPeaks <- rbind(rPeaks, availablePeaks[sample(nrow(availablePeaks), 1),])
			}
			### Remove selected genes from pool
			fullSet <- fullSet[-which(fullSet$peakID %in% rPeaks$peakID),]
		}
		return(list(peaks=rPeaks,nonMatched=nonMatched))
	}
}


sampleMatchedExpression <- function(genes, expression, eMargin=0.1){
	
	allGenes <- gintervals.load('glp.intervals.ucscCanTSS')
	rownames(allGenes) <- allGenes$geneName
	
	genes[is.na(genes)] <- 0
	
	chr <- as.character(genes[1,1])
	
	fullSet <- allGenes[allGenes$chrom == chr,]
	fullSet <- fullSet[-which(fullSet$geneName %in% genes$geneName),]
	fullSet$expression <- expression[fullSet$geneName]
	
	
	fullSet[is.na(fullSet)] <- 0
	
	rGenes <- c()
	
	for (i in 1:nrow(genes))
	{
		eLevel <- as.numeric(as.character(genes[i,'expression']))
		### Select set of suitable genes
		availableGenes <- fullSet[fullSet$expression >= eLevel-eMargin*eLevel & fullSet$expression <= eLevel+eMargin*eLevel, ]
		if(nrow(availableGenes) == 0)
		{
			availableGenes <- fullSet[fullSet$expression >= eLevel-eMargin*eLevel*3 & fullSet$expression <= eLevel+eMargin*eLevel*3, ]
			rGenes <- rbind(rGenes, availableGenes[sample(nrow(availableGenes), 1),])
		}else{
			rGenes <- rbind(rGenes, availableGenes[sample(nrow(availableGenes), 1),])
		}
		
		### Remove selected genes from pool
		fullSet <- fullSet[-which(fullSet$geneName %in% rGenes$geneName),]
	}
	return(rGenes)
}


clipQuantile <- function(set, thr, up=TRUE, down=TRUE)
{
	set <- as.numeric(as.character(set))
# 	set[is.na(set)] <- 0
	upLim <- quantile(set,thr,na.rm=TRUE)
	downLim <- quantile(set,1-thr,na.rm=TRUE)
	
# 	print(paste(downLim,upLim))
	
	set[set <= downLim] <- downLim
	set[set >= upLim] <- upLim
	
	return(set)
}

clipMatrix <- function(mtx, thr, up=TRUE, down=TRUE){
	
	mtx[is.na(mtx)] <- 0
	upLim <- quantile(mtx,thr,na.rm=TRUE)
	downLim <- quantile(mtx,1-thr,na.rm=TRUE)
	
	if(down){mtx[mtx <= downLim] <- downLim}
	if(up){mtx[mtx >= upLim] <- upLim}
	
	return(mtx)
}

get_borders <- function(iterator=2e3, ins_track=ins_track, ins_dom_thresh=5e-2, extend_border_region=2e3, mapFilter=FALSE, posValues=FALSE){
	
	if(posValues == TRUE){
		ins <- gscreen(sprintf("(!is.na(%s)) & (%s >= %f) ", ins_track, ins_track, gquantiles(ins_track, ins_dom_thresh)), iterator=iterator)
	}else{
		ins <- gscreen(sprintf("(!is.na(%s)) & (%s <= %f) ", ins_track, ins_track, gquantiles(ins_track, ins_dom_thresh)), iterator=iterator)
	}
	
	
	chroms <- gintervals.all()
	rownames(chroms) <- as.character(chroms$chrom)
	ins$max <- chroms[ as.character(ins$chrom), 'end']
	ins$start <- pmax(ins$start - extend_border_region, 0)
	ins$end  <- pmin(ins$end   + extend_border_region, ins$max)
	insc <- gintervals.canonic(ins, unify=T)
# 	print(head(insc))
	
	if(posValues == TRUE){
		gvtrack.create("ins_max", ins_track, "max")
		bords <- gextract("ins_max", intervals=insc, iterator=insc)
		bords <- bords[(!is.na(bords$ins_max))&(bords$ins_max > gquantiles(ins_track, ins_dom_thresh)),]
		gvtrack.rm("ins_max")
	}else{
		gvtrack.create("ins_min", ins_track, "min")
		bords <- gextract("ins_min", intervals=insc, iterator=insc)
		bords <- bords[(!is.na(bords$ins_min))&(bords$ins_min < gquantiles(ins_track, ins_dom_thresh)),]
		gvtrack.rm("ins_min")
	}
	
	return(bords[,1:3])
	
# 	if(mapFilter == TRUE)
# 	{
# 		print('Filtering...')
# 		mapTrack <- 'redb.feGATC_map'
# 		gvtrack.create("vMap", mapTrack, "avg")
# 		mData <- gextract("vMap", insc, iterator=insc)
# 		mData <- mData[complete.cases(mData),]
# 		mData <- mData[mData[,"vMap"] >= 2.5,]
# 		insc <- mData[,1:3]
# 		gvtrack.rm("vMap")
# 		print(nrow(insc))
# 		bords <- gextract("ins_min", intervals=insc, iterator=insc)
# 		print(nrow(bords))
# 		bords <- bords[complete.cases(bords),]
# 		print(nrow(bords))
# 		bords <- bords[bords$ins_min < quantile(bords$ins_min, ins_dom_thresh),]
# 		print(nrow(bords))
# 		print(quantile(bords$ins_min, ins_dom_thresh))
# 		print(gquantiles(ins_track, ins_dom_thresh))
# 		
# 	}else{
# 		bords <- gextract("ins_min", intervals=insc, iterator=insc)
# 		bords <- bords[(!is.na(bords$ins_min))&(bords$ins_min < gquantiles(ins_track, ins_dom_thresh)),]
# 	}
	
# 	gvtrack.rm("ins_min")
	
}


callPeaks <- function(data, probThr=0.01,scoreThr=5, smooth=FALSE,setName='res',tolerance=0.025,plot=FALSE,eqPlot=TRUE){
	
	if(smooth==TRUE){d <- smooth.spline(data, tol=tolerance, keep.data=T)
	}else{d <- spline(data, n=length(data))}
	
	tp <- turnpoints(d$y)
	tp.pos <- tp$tppos
	tp.proba <- tp$proba
	keep <- which(tp$peaks == TRUE)
	tp.pos <- tp.pos[tp.proba < probThr]
	tp.proba <- tp.proba[tp.proba < probThr]
	tp.pos <- intersect(tp.pos, keep)
	tp.proba <- tp.proba[which(tp.pos %in% keep)]
	
	posBins <- data[d$x[tp.pos]]
	posBins <- cbind(posBins,tp.proba)
	colnames(posBins) <- c(setName,'pVal')
	posBins <- as.data.frame(posBins[posBins[,setName] > scoreThr,])
	if(ncol(posBins) == 1){posBins <- as.data.frame(t(posBins))}
	
	posBins$refID <- sapply(posBins$res, function(x){which(data == x)[1]})
	
	if(plot == TRUE){
		if(min(data) >= 0){eqPlot == FALSE}
		if(eqPlot == TRUE){
			y <- max(abs(data))
			plot(data, type='l', lwd=2, col=4, ylim=c(-y,y))
			abline(h=0,lty=2,col='lightgrey')
		}else{
			plot(data, type='l', lwd=2, col=4)
		}
		lines(d$y, lwd=1.5, col=1, lty=2)
		points(x=posBins$refID, y=posBins$res, pch=19, col=2)
	}
	
	
	return(list(peaks=posBins,data=d,nP=length(tp.pos)))
}


bonevScale <- function(){
### boyan
	palette.breaks = function(n, colors, breaks){
		colspec = colorRampPalette(c(colors[1],colors[1]))(breaks[1])
		for(i in 2:(length(colors)) ){
			colspec = c(colspec, colorRampPalette(c(colors[i-1], colors[i]))(abs(breaks[i]-breaks[i-1])))
		}
		colspec = c( colspec,
			colorRampPalette(c(colors[length(colors)],colors[length(colors)]))(n-breaks[length(colors)])
			)
		colspec
	}
	col.pbreaks <- c(20,35,50,65,75,85,95)
	col.pos <- palette.breaks(100 , c("lightgrey","lavenderblush2","#f8bfda","lightcoral","red","orange","yellow"), col.pbreaks)
	col.nbreaks <- c(20,35,50,65,75,85,95)
	col.neg <- rev(palette.breaks(100 , c("lightgrey", "powderblue", "cornflowerblue", "blue","blueviolet", "mediumpurple4", "black"), col.nbreaks ))
	scoreCol <- c(col.neg, col.pos)
	return(scoreCol)
}

build2D <- function(set1,set2,names1='',names2='',minDist=0,maxDist=10e10,band=TRUE,diag=TRUE,retNames=TRUE, trans=FALSE){
# 	print(c(nrow(set1),nrow(set2)))
	####################################################
# 	chrs <- unique(as.character(rbind(set1,set2)$chrom))
	chrs <- unique(c(as.character(set1$chrom),as.character(set2$chrom)))
	ints <- c()
	for(chr in chrs){
		d1 <- set1[as.character(set1$chrom) == chr,]
		d2 <- set2[as.character(set2$chrom) == chr,]
		
		if(nrow(d1) == 0 | nrow(d2) == 0){next}
		
		d_inds = as.data.frame(expand.grid(1:nrow(d1), 1:nrow(d2)))
		d_ints_for = data.frame(chrom1=d1[d_inds[,1],'chrom'], start1=d1[d_inds[,1],'start'], end1=d1[d_inds[,1],'end'], chrom2=d2[d_inds[,2],'chrom'], start2=d2[d_inds[,2],'start'], end2=d2[d_inds[,2],'end'])
		d_ints_rev = data.frame(chrom1=d2[d_inds[,2],'chrom'], start1=d2[d_inds[,2],'start'], end1=d2[d_inds[,2],'end'], chrom2=d1[d_inds[,1],'chrom'], start2=d1[d_inds[,1],'start'], end2=d1[d_inds[,1],'end'])
		d_ints <- rbind(d_ints_for,d_ints_rev)
		d_ints <- d_ints[!duplicated(paste(d_ints$start1,d_ints$end1,d_ints$start2,d_ints$end2,sep='_')),]
		##################
		d_ints <- gintervals.canonic(d_ints)
		
		if(is.null(ints)){ints <- d_ints
		}else{ints <- rbind(ints,d_ints)}
	}
# 	print(nrow(ints))
	if(is.null(ints)){return(NA)}
	####################################################
	if(band==TRUE){ints <- ints[ints$start2 >= ints$start1,]}
	
	ints$dist <- abs(ints$start2 - ints$start1)
	ints <- ints[ints$dist >= minDist & ints$dist <= maxDist,]
	
	if(diag==FALSE){ints <- ints[ints$dist != 0,]}
	
	if(trans==TRUE){
		for(chr in chrs){
		print(chr)
			d1 <- set1[as.character(set1$chrom) == chr,]
			d2 <- set2[as.character(set2$chrom) != chr,]
			
			if(nrow(d1) == 0 | nrow(d2) == 0){next}
			
			d_inds = as.data.frame(expand.grid(1:nrow(d1), 1:nrow(d2)))
			d_ints_for = data.frame(chrom1=d1[d_inds[,1],'chrom'], start1=d1[d_inds[,1],'start'], end1=d1[d_inds[,1],'end'], chrom2=d2[d_inds[,2],'chrom'], start2=d2[d_inds[,2],'start'], end2=d2[d_inds[,2],'end'])
			d_ints_rev = data.frame(chrom1=d2[d_inds[,2],'chrom'], start1=d2[d_inds[,2],'start'], end1=d2[d_inds[,2],'end'], chrom2=d1[d_inds[,1],'chrom'], start2=d1[d_inds[,1],'start'], end2=d1[d_inds[,1],'end'])
			d_ints <- rbind(d_ints_for,d_ints_rev)
			d_ints <- d_ints[!duplicated(paste(d_ints$chrom1,d_ints$start1,d_ints$end1,d_ints$chrom2,d_ints$start2,d_ints$end2,sep='_')),]
			d_ints$dist <- NA
			##################
# 			d_ints <- gintervals.canonic(d_ints)
			
			print(head(d_ints))
			
			if(is.null(ints)){ints <- d_ints
			}else{
				ints <- rbind(ints,d_ints)
				ints <- ints[!duplicated(paste(ints$chrom1,ints$start1,ints$end1,ints$chrom2,ints$start2,ints$end2,sep='_')),]
			}
		}
	}
	
	if(retNames == TRUE){
		
		if(names1[1]==''){names1 <- paste(paste0('p',1:nrow(set1)),set1$chrom,set1$start,set1$end,sep='_')}
		if(names2[1]==''){names2 <- paste(paste0('p',1:nrow(set2)),set2$chrom,set2$start,set2$end,sep='_')}
		
		outNames1 <- c()
		outNames2 <- c()
		for (i in 1:nrow(ints)){
			x <- ints[i,]
			s1 <- paste(x['start1'],x['end1'],sep='_')
			s2 <- paste(x['start2'],x['end2'],sep='_')
			oN1 <- unlist(strsplit(names1[grep(s1,names1)],'_'))[1]
			oN2 <- unlist(strsplit(names2[grep(s2,names2)],'_'))[1]
			outNames1 <- append(outNames1,oN1)
			outNames2 <- append(outNames2,oN2)
		}
		
		ints$name1 <- outNames1
		ints$name2 <- outNames2
# 		ints$name1 <- apply(ints,1,function(x){ unlist(strsplit(names1[grep(paste(x['start1'],x['end1'],sep='_'),names1)],'_'))[1] })
# 		ints$name2 <- apply(ints,1,function(x){ unlist(strsplit(names2[grep(paste(x['start2'],x['end2'],sep='_'),names2)],'_'))[1] })
	}
# 	print(nrow(ints))
	return(ints)
}

buildTrans <- function(set1,set2){
	####################################################
	chrs <- unique(c(as.character(set1$chrom),as.character(set2$chrom)))
	ints <- c()
	for(chr in chrs){
	print(chr)
		d1 <- set1[as.character(set1$chrom) == chr,]
		d2 <- set2[as.character(set2$chrom) != chr,]
		
		if(nrow(d1) == 0 | nrow(d2) == 0){next}
		
		d_inds = as.data.frame(expand.grid(1:nrow(d1), 1:nrow(d2)))
		d_ints_for = data.frame(chrom1=d1[d_inds[,1],'chrom'], start1=d1[d_inds[,1],'start'], end1=d1[d_inds[,1],'end'], chrom2=d2[d_inds[,2],'chrom'], start2=d2[d_inds[,2],'start'], end2=d2[d_inds[,2],'end'])
		d_ints_rev = data.frame(chrom1=d2[d_inds[,2],'chrom'], start1=d2[d_inds[,2],'start'], end1=d2[d_inds[,2],'end'], chrom2=d1[d_inds[,1],'chrom'], start2=d1[d_inds[,1],'start'], end2=d1[d_inds[,1],'end'])
		d_ints <- rbind(d_ints_for,d_ints_rev)
		d_ints <- d_ints[!duplicated(paste(d_ints$chrom1,d_ints$start1,d_ints$end1,d_ints$chrom2,d_ints$start2,d_ints$end2,sep='_')),]
		##################
		print(head(d_ints))
		
		if(is.null(ints)){ints <- d_ints
		}else{
			ints <- rbind(ints,d_ints)
# 			ints <- ints[!duplicated(paste(ints$chrom1,ints$start1,ints$end1,ints$chrom2,ints$start2,ints$end2,sep='_')),]
		}
	}
	return(ints)
}

powerLawFit <- function(g){
	dP <-  degree.distribution(g, mode = "all", cumulative = FALSE)
	dP <- dP[-1]
	x <- 1:max(degree(g))
	obsDegrees = which(dP != 0)
	dP <- dP[obsDegrees]
	x <- x[obsDegrees]
	y <- dP
	fit <- lm(log(y) ~ log(x))
	params <- coef(fit)
	exps <- append(exps,params[2])
	rsq <- signif(summary(fit)$adj.r.squared, 5)
	
	return(list(model=fit,rsq=rsq))
}


plotPoints <- function(track, interval, colors=colorRampPalette(c('lightgrey','darkred'))(100), inverse=FALSE, mirror=FALSE, add=FALSE, tri=FALSE, axis=FALSE, band='false', title='',pointCEX=0.25){
	
	
# 	if(nrow(track) < 50000){pointCEX <- 0.5;}
	
	if(tri == FALSE){
		### Squares
		track <- track[track$start2 > track$start1,]
		track <- track[order(track[,7]),]
		xAxis <- c(interval$start1,interval$end1)
		yAxis <- c(interval$start2,interval$end2)
		if(add == FALSE){
			plot(0,xlim=xAxis,ylim=yAxis, main=title, col='white',cex.main=2,xlab='', ylab='', yaxt='n', xaxt='n', frame=F)
			if(axis == TRUE){axis(1, cex=1.4); axis(2, cex=1.4); }
		}
		
		apply(track, 1, function(x)
		{
			col <- colors[as.numeric(as.character(x[7]))+101]
			xx <- as.numeric(as.character(x['start1']))
			yy <- as.numeric(as.character(x['start2']))
			
			if(add == FALSE){
				points(x=xx, y=yy, type='p', cex=pointCEX, col=col, pch=19)
# 				polygon(xx,yy, lwd=1, col=col,border=col)
			}
			
			if(mirror == TRUE | inverse == TRUE){
				points(x=yy, y=xx, type='p', cex=pointCEX, col=col, pch=19)
			}
			
		})
	}else{
		### Triangle
		track <- track[track$start2 >= track$start1,]
		xAxis <- c(interval$start1,interval$end2)
		winSize <- xAxis[2]-xAxis[1]
		yAxis <- c(0,winSize)
		plot(0,xlim=xAxis,ylim=yAxis, col='white',cex.main=2,xlab='', ylab='', yaxt='n', xaxt='n', frame=F)
# 		par(mai=rep(0.25,4))
		if(axis == TRUE){axis(1, cex=1.4)}
		apply(track, 1, function(x)
		{
			col <- colors[as.numeric(as.character(x[7]))+101]
			
			start <- as.numeric(as.character(x['start1']))
			end <- as.numeric(as.character(x['start2']))
			dist <- end-start
			mid <- min(c(start,end))+abs(dist)/2
			
			xx <- mid
			yy <- abs(dist)
			
			bandShift <- 3.6
			if(band == 'true'){yy <- yy*bandShift}
			
			points(x=xx, y=yy, type='p', cex=pointCEX, col=col, pch=19)
		})
		if(title != ''){legend("topleft", legend=title, cex=1.6)}
		
	}
}



colScalePlot <- function(x, labels=c(-1,0,1), pal=colorRampPalette(c("blue","#87FFFF", "white", "#FF413D", "black")), pos=1, vert=FALSE, main=''){
	if(vert==TRUE){
		plot(y=1:x,x=rep(0,x),col=pal(x), pch=15, xlab='', ylab='', yaxt='n', xaxt='n',frame=F,ylim=c(-x*0.2,x*1.2),xlim=c(-5,5), main=main);
		sapply(1:length(labels), function(y){text(y=seq(0,1,length.out=length(labels))[y]*x,x=0,labels=labels[y],pos=pos)});
	}else{
		plot(x=1:x,y=rep(0,x),col=pal(x), pch=15, xlab='', ylab='', yaxt='n', xaxt='n',frame=F,xlim=c(-x*0.2,x*1.2),ylim=c(-5,5), main=main);
		sapply(1:length(labels), function(y){text(x=seq(0,1,length.out=length(labels))[y]*x,y=0,labels=labels[y],pos=pos)});
	}
}


plotRectangles2D <- function(data, int2d=limits2d, fc=0.5, col='yellow', lwd=0.5, lty=1, image=TRUE, side='both', shiftEnd=0, shiftStart=0,shade=FALSE){
	####
# 	data <- gintervals.intersect(data,int2d)
# 	data <- gintervals.neighbors(data,int2d,maxdist1=0,maxdist2=0,mindist1=0,mindist2=0)[,1:ncol(data)]
	
	
	
# 	print(nrow(data))
	
	if(is.null(data)){return()}
	
	#### Normalize @ 1 by 1 for image...
	xLims <- c(int2d[2:3])
	xL <- as.numeric(as.character(int2d[3] - int2d[2]))
	yLims <- c(int2d[5:6])
	yL <- as.numeric(as.character(int2d[6] - int2d[5]))
	
	start <- as.numeric(as.character(xLims[1]))
	end <- as.numeric(as.character(xLims[2]))
	
# 	if('fc' %in% colnames(data)){
# 		data <- data[data$fc > fc,]
# 		if(nrow(data) < 1){return()}
# 	}
# 	print(nrow(data))
	for(g in 1:nrow(data))
	{
		row <- data[g,]
# 		if(row$start1 < start | row$end1 > end){next}
		if(image==FALSE)
		{
			x1 <- row$start1
			y1 <- row$end1
			x2 <- row$start2
			y2 <- row$end2
		}else if (image==TRUE){
			x1 <- (row$start1 - as.numeric(as.character(xLims[1])))/xL
			y1 <- (row$end1 - as.numeric(as.character(xLims[1])))/xL
			x2 <- (row$start2 - as.numeric(as.character(yLims[1])))/yL
			y2 <- (row$end2 - as.numeric(as.character(yLims[1])))/yL
		}else if (image=='TRI'){
			
			if(row$start1 != row$start2){
				
				x1 <- row$start1 + (row$start2-row$start1)/2
				x2 <- row$start1 + (row$end2-row$start1)/2
				x3 <- row$end1 + (row$end2-row$end1)/2
				x4 <- row$end1 + (row$start2-row$end1)/2
				
				y1 <- row$start2-row$start1
				y2 <- row$end2-row$start1
				y3 <- row$end2-row$end1
				y4 <- row$start2-row$end1
				
				
				segments(x1,y1,x2,y2,col=col,lwd=lwd,lty=lty)
				segments(x2,y2,x3,y3,col=col,lwd=lwd,lty=lty)
				
				segments(x3,y3,x4,y4,col=col,lwd=lwd,lty=lty)
				segments(x4,y4,x1,y1,col=col,lwd=lwd,lty=lty)
				if(shade == TRUE){
					print('filled...')
					polygon(c(x1,x2,x3,x4),c(y1,y2,y3,y4),col=1)
				}
			}else{
				x1 <- row$start1
				x2 <- row$end2
				mid <- x1 + (x2-x1)/2
				y1 <- 0
				y2 <- (x2-x1)
				
				
				if(nrow(data) > 1){
					if(g == 1){
						x1 <- row$end2
						x2 <- row$end2
						mid <- start + (x2-start)/2
						y1 <- 0
						y2 <- (x2-start)
					}else if(g == nrow(data)){
						x1 <- row$start2
						x2 <- row$start2
						mid <- x1 + (end-x1)/2
						y1 <- 0
						y2 <- (end-x1)
					}
				}
				
				if(side=='top'){
					segments(x1,y1,mid,y2,col=col,lwd=lwd,lty=lty)
					segments(mid,y2,x2,y1,col=col,lwd=lwd,lty=lty)
				}else if(side=='bottom'){
					segments(x2,y1,mid,-y2,col=col,lwd=lwd,lty=lty)
					segments(mid,-y2,x1,y1,col=col,lwd=lwd,lty=lty)
				}else{
					segments(x1,y1,mid,y2,col=col,lwd=lwd,lty=lty)
					segments(mid,y2,x2,y1,col=col,lwd=lwd,lty=lty)
					segments(x2,y1,mid,-y2,col=col,lwd=lwd,lty=lty)
					segments(mid,-y2,x1,y1,col=col,lwd=lwd,lty=lty)
				}
			}
		}
		
		if (image!='TRI'){
			
			if(side=='top'){
				segments(x1,x2,x1,y2,col=col,lwd=lwd,lty=lty)
				segments(x1,y2,y1,y2,col=col,lwd=lwd,lty=lty)
			}else if(side=='bottom'){
				segments(y1,y2,y1,x2,col=col,lwd=lwd,lty=lty)
				segments(x1,x2,y1,x2,col=col,lwd=lwd,lty=lty)
			}else{
				segments(x1,x2,x1,y2,col=col,lwd=lwd,lty=lty)
				segments(x1,y2,y1,y2,col=col,lwd=lwd,lty=lty)
				segments(y1,y2,y1,x2,col=col,lwd=lwd,lty=lty)
				segments(x1,x2,y1,x2,col=col,lwd=lwd,lty=lty)
			}
		}
	}
}

expand1D <- function(set,expand,center=1){
	require(plyr)
	set$dist <- set$end - set$start
	
	set$mid <- set$start + round_any((set$end-set$start)/2,center)
	
	set <- set[order(set$mid),]
# 	print('here')
	
	mids <- set$mid
	min <- 10e10
	
	if(nrow(set) > 1){
		for (i in 1:(length(mids)-1)){
			dist <- mids[i+1] - mids[i]
			if(dist < min ){min <- dist}
		}
		
		if(min < expand*2){
			print(paste('************************************************',min))
		}
	}
	
	set$start <- set$mid - expand
	set$end <- set$mid + expand
	
	set <- set[set$start > 0,]
	
	for (chr in as.character(gintervals.all()$chrom)){
		maxSize <- as.numeric(as.character(gintervals.all()[gintervals.all()$chrom == chr,'end']))
		if(length( which(set$chrom == chr & set$end > maxSize) ) > 0){
			set <- set[-which(set$chrom == chr & set$end > maxSize),]
		}
	}
	
	set <- gintervals.canonic(set)
	
	return(set)
}

expand2D <- function(set,expand,selectOn=NULL,higher=TRUE){
	
	set[is.na(set)] <- 0
	
	set$start1 <- set$start1 - expand
	set$end1 <- set$end1 + expand
	set$start2 <- set$start2 - expand
	set$end2 <- set$end2 + expand
	
# 	set$end1 <- set$end1 - 1
# 	set$end2 <- set$end2 - 1
	
	set <- set[set$start1 > 0,]
	
	if(nrow(set) == 0){
		return(NULL)
	}
	
	for (chr in as.character(gintervals.all()$chrom)){
		maxSize <- as.numeric(as.character(gintervals.all()[gintervals.all()$chrom == chr,'end']))
		if(length( which(set$chrom1 == chr & set$end2 > maxSize) ) > 0){
			set <- set[-which(set$chrom1 == chr & set$end2 > maxSize),]
		}
	}
	
	set$mid1 <- set$start1 + round_any((set$end1-set$start1)/2,1)
	set$mid2 <- set$start2 + round_any((set$end2-set$start2)/2,1)
	
	if(nrow(set) == 1){
		return(set)
	}
	
	print('starting')
	print(nrow(set))
	
	done <- 1
	while(done > 0){
		
# 		set <- set[order(set$chrom1,set$mid2,set$mid1),]
# 		set <- set[order(set$chrom1,set$mid1,set$mid2),]
		set <- set[order(set$chrom1,set$mid1+set$mid2),]
		print(nrow(set))
		remove <- c()
		d1 <- 0
		d2 <- 0
		
		for (i in 1:(nrow(set)-1)){
			if(set[i,'chrom1'] != set[i+1,'chrom1']){next}
			
			d1 <- set[i+1,'mid1'] - set[i,'mid1']
			d2 <- set[i+1,'mid2'] - set[i,'mid2']
			
			if(d1 <= expand*2 & d2 <= expand*2){
				
				if(!is.null(selectOn)){
					if(set[i,which(colnames(set)==selectOn)] >= set[i+1,which(colnames(set)==selectOn)]){
						
						if(higher == TRUE){remove <- append(remove,i+1)
						}else{remove <- append(remove,i)}
					}else{
						if(higher == TRUE){remove <- append(remove,i)
						}else{remove <- append(remove,i+1)}
					}
				}else{
					remove <- append(remove,i)
				}
			}
# 			if(length(remove) > 0){break}
		}
# 		print(remove)
		
		remove <- unique(remove)
		
		set$d1 <- d1
		set$d2 <- d2
		done <- length(remove)
# 		print(done)
		
		if(done != 0){
			set <- set[-remove,]
		}
	}
	
	done <- 1
	while(done > 0){
		
		set <- set[order(set$chrom1,set$mid2,set$mid1),]
# 		set <- set[order(set$chrom1,set$mid1,set$mid2),]
# 		set <- set[order(set$chrom1,set$mid1+set$mid2),]
		print(nrow(set))
		remove <- c()
		d1 <- 0
		d2 <- 0
		
		for (i in 1:(nrow(set)-1)){
			if(set[i,'chrom1'] != set[i+1,'chrom1']){next}
			
			d1 <- set[i+1,'mid1'] - set[i,'mid1']
			d2 <- set[i+1,'mid2'] - set[i,'mid2']
			
			if(d1 <= expand*2 & d2 <= expand*2){
				
				if(!is.null(selectOn)){
					if(set[i,which(colnames(set)==selectOn)] >= set[i+1,which(colnames(set)==selectOn)]){
						
						if(higher == TRUE){remove <- append(remove,i+1)
						}else{remove <- append(remove,i)}
					}else{
						if(higher == TRUE){remove <- append(remove,i)
						}else{remove <- append(remove,i+1)}
					}
				}else{
					remove <- append(remove,i)
				}
			}
# 			if(length(remove) > 0){break}
		}
# 		print(remove)
		
		remove <- unique(remove)
		
		set$d1 <- d1
		set$d2 <- d2
		done <- length(remove)
# 		print(done)
		
		if(done != 0){
			set <- set[-remove,]
		}
	}
	
# 	done <- 1
# 	while(done > 0){
# 		
# # 		set <- set[order(set$chrom1,set$mid2,set$mid1),]
# 		set <- set[order(set$chrom1,set$mid1,set$mid2),]
# # 		set <- set[order(set$chrom1,set$mid1+set$mid2),]
# 		print(nrow(set))
# 		remove <- c()
# 		d1 <- 0
# 		d2 <- 0
# 		
# 		for (i in 1:(nrow(set)-1)){
# 			if(set[i,'chrom1'] != set[i+1,'chrom1']){next}
# 			
# 			d1 <- set[i+1,'mid1'] - set[i,'mid1']
# 			d2 <- set[i+1,'mid2'] - set[i,'mid2']
# 			
# 			if(d1 <= expand*2 & d2 <= expand*2){
# 				
# 				if(!is.null(selectOn)){
# 					if(set[i,which(colnames(set)==selectOn)] >= set[i+1,which(colnames(set)==selectOn)]){
# 						
# 						if(higher == TRUE){remove <- append(remove,i+1)
# 						}else{remove <- append(remove,i)}
# 					}else{
# 						if(higher == TRUE){remove <- append(remove,i)
# 						}else{remove <- append(remove,i+1)}
# 					}
# 				}else{
# 					remove <- append(remove,i)
# 				}
# 			}
# # 			if(length(remove) > 0){break}
# 		}
# # 		print(remove)
# 		
# 		remove <- unique(remove)
# 		
# 		set$d1 <- d1
# 		set$d2 <- d2
# 		done <- length(remove)
# # 		print(done)
# 		
# 		if(done != 0){
# 			set <- set[-remove,]
# 		}
# 	}
# 	
	return(set)
}

# To compute the length of the peaks
length1D <- function(set,accuracy=1){
	require(plyr)
	set$length <- round_any(set$end - set$start,accuracy)
	set$start <- set$length
	set$end <- set$length
	
	return(set)
}

# To center the peaks
center1D <- function(set,center=1){
	require(plyr)
	set$mid <- round_any(set$start + (set$end-set$start)/2,center)
	set$start <- set$mid
	set$end <- set$mid + 1
	
	return(set)
}

center2D <- function(set,center=1){
	require(plyr)
	set$mid1 <- round_any(set$start1 + (set$end1-set$start1)/2,center)
	set$mid2 <- round_any(set$start2 + (set$end2-set$start2)/2,center)
	
	set$start1 <- set$mid1
	set$end1 <- set$mid1 + 1
	
	set$start2 <- set$mid2
	set$end2 <- set$mid2 + 1
	
	return(set)
}

plotReadRhombus <- function(track, interval, binSize=2.5e+3, colors=colorRampPalette(c('lightgrey','darkred'))(100), inverse=FALSE, mirror=FALSE, add=FALSE, tri=FALSE, axis=FALSE, band='false', plotClmn=7, title='', cex=1.6){
	
	if(plotClmn != 7){plotClmn <- which(colnames(track) == plotClmn)}
	
	print('plotting')
	
	if(tri == FALSE){
		### Squares
# 		track <- track[track$start2 > track$start1,]
		track <- track[track$start2 > track$start1 | track$start2 == track$start1,]
		xAxis <- c(interval$start1,interval$end1)
		yAxis <- c(interval$start2,interval$end2)
		if(add == FALSE){
			plot(0,xlim=xAxis,ylim=yAxis, main=title, col='white',cex.main=2,xlab='', ylab='', yaxt='n', xaxt='n', frame=F)
		}
		
		apply(track, 1, function(x)
		{
			xx <- c(as.numeric(as.character(x['start1'])),as.numeric(as.character(x['start1'])),as.numeric(as.character(x['end1'])),as.numeric(as.character(x['end1'])))
			yy <- c(as.numeric(as.character(x['start2'])),as.numeric(as.character(x['end2'])),as.numeric(as.character(x['end2'])),as.numeric(as.character(x['start2'])))
			col <- colors[as.numeric(as.character(x[plotClmn]))+101]
			if(add == FALSE){
				polygon(xx,yy, lwd=1, col=col,border=col)
			}
			
			if(mirror == TRUE | inverse == TRUE){
				polygon(yy,xx, lwd=1, col=col,border=col, ljoin=1,lwd=0.05)
			}
			
		})
		
		if(axis==TRUE){axis(1, cex=1.6); axis(2, cex=1.6);}
		
	}else{
		### Triangle
		track <- track[track$start2 >= track$start1,]
		xAxis <- c(interval$start1,interval$end2)
		winSize <- xAxis[2]-xAxis[1]
		yAxis <- c(0,winSize)
		
		if(band!='false'){yAxis[2]=yAxis[2]*band}
		
		
		track$dist <- abs(track$start2 - track$start1)
		track <- track[track$dist <= yAxis[2],]
		
		plot(0,xlim=xAxis,ylim=yAxis, col='white',cex.main=2,xlab='', ylab='', yaxt='n', xaxt='n', frame=F)
# 		rect(xAxis[1],yAxis[1],xAxis[2],yAxis[2],col='lightgrey')
# 		points(x=seq(xAxis[1],xAxis[2], length.out=length(colors)),y=rep(yAxis[1]-binSize*2,length(colors)),type='p', cex=3, col=colors[1:length(colors)], pch=19)
		apply(track, 1, function(x)
		{
			start <- as.numeric(as.character(x['start1']))+binSize/2
			end <- as.numeric(as.character(x['start2']))+binSize/2
			dist <- end-start
			
			mid <- min(c(start,end))+abs(dist)/2
			xx <- c(mid-binSize/2,mid,mid+binSize/2,mid)
			height <- abs(dist)
			yy <- c(height,height+binSize,height,height-binSize)
			
			if(yy[4] < 0){yy[4] <- 0}
			if(yy[2] > yAxis[2]+binSize){yy[2] <- yAxis[2]}
			col <- colors[as.numeric(as.character(x[plotClmn]))+101]
			
			polygon(xx,yy, lwd=1, col=col,border=col,ljoin=1,lwd=0.05)
# 			polygon(xx,yy, lwd=1, col=col,border=NA)
		})
# 		if(title != ''){legend("topleft", legend=title, cex=1.6)}
# 		if(title != ''){text(x=xAxis[1]+(winSize*0.025), y=yAxis[2]*0.8, labels=title, cex=1.6, pos=4)}
		if(title != ''){text(x=xAxis[1], y=yAxis[2]*0.8, labels=title, cex=cex, pos=4)}
		if(axis==TRUE){axis(1, cex=1.2);}
	}
}


genesPlot <- function(set1, plotLim, vertical=FALSE, header='test', cex=1.4, rHeight=40){
# 	print(dim(set1))
	if(vertical == TRUE)
	{
# 		par(mar = c(plotMar[2], 0, plotMar[1], 0))
		xlim=c(-150,150)
		plot(0, xlab='', ylab='', yaxt='n', xaxt='n', frame=F, ylim=plotLim, xlim=xlim)
		if(is.null(set1)){return()}
# 		print(plotLim)
		segments(0,plotLim[1],0,plotLim[2],lwd=1.75, col='darkgrey')

		for(g in 1:nrow(set1))
		{
			row <- set1[g,]
			if(row$strand == -1){rect(5,row$start, 50, row$end, col=rgb(1,1,0, alpha=0.9), border = 1); textY=row$end; currentX <- 125; textPos <- 3;
			}else{rect(-5, row$start, -50, row$end, col=rgb(0,1,0, alpha=0.5), border = 1); textY =row$start; currentX <- -125; textPos <- 1;}
			
			if(nrow(set1) < 20){
# 				text(x=currentX,y=textY, labels=row$geneName, pos=textPos, cex=cex)
				text(x=currentX,y=textY, labels=row$geneName, pos=textPos, cex=1.2, srt = 90)
			}
		}
		
		text(x=0,y=plotLim[2], labels=unique(set1$chrom), pos=3, cex=cex, col='grey15')
	}else{
# 		par(mar = c(0, plotMar[1], 0, plotMar[2]))
		yLim=c(-150,150)
		plot(0, xlim=plotLim, ylim=yLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F)
		if(is.null(set1)){return()}
		segments(plotLim[1],0,plotLim[2],0,lwd=1.75, col='darkgrey')
		
		minH <- 5
		arExt <- (plotLim[2]-plotLim[1])/125
		arH <- rHeight*1.3
		
		for(g in 1:nrow(set1))
		{
			row <- set1[g,]
			if(row$geneName == header){rect(row$start, yLim[1], row$end, yLim[2], col=rgb(0.66,0.66,0.66, alpha=0.4), border = 0)}
			
			if(row$strand == 1){
				rect(row$start, minH, row$end, rHeight, col=rgb(0,1,0, alpha=0.9), border = 1); currentX=row$start; textY <- arH; textPos <- 3;
				segments(row$start,minH,row$start,arH)
				segments(row$start,arH,row$start+arExt,arH)
				
				segments(row$start+arExt*0.75,arH+arH*0.15,row$start+arExt,arH)
				segments(row$start+arExt*0.75,arH-arH*0.15,row$start+arExt,arH)
				
			}else{
				rect(row$start, -minH, row$end, -rHeight, col=rgb(1,1,0, alpha=0.9), border = 1); currentX=row$end; textY <- -arH; textPos <- 1;
				segments(row$end,-minH,row$end,-arH)
				segments(row$end,-arH,row$end-arExt,-arH)
				
				segments(row$end-arExt*0.75,-arH+arH*0.15,row$end-arExt,-arH)
				segments(row$end-arExt*0.75,-arH-arH*0.15,row$end-arExt,-arH)
			}
			
			if(nrow(set1) < 50){
				text(x=currentX,y=textY, labels=row$geneName, pos=textPos, cex=cex)
			}
		}
# 		text(x=plotLim[1],y=0, labels=unique(set1$chrom), pos=2, cex=1.4, col='grey15')
	}
}

genesExpPlot <- function(genes, expData, plotLim, cells, cex=1.6, scaleAll=TRUE, log=FALSE){
	
	yLim=c(0,150)
	plot(0, xlim=plotLim, ylim=yLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F)
	if(is.null(genes)){return()}
	segments(plotLim[1],0,plotLim[2],0,lwd=1.75, col='darkgrey')
	
	set1 <- expData[rownames(expData) %in% genes$geneName,cells]
	if(log==TRUE){set1 <- log10(set1+1)}
	me <- max(set1)
	barWidth <- (plotLim[2]-plotLim[1])/300
	
	for(gene in rownames(set1)){
		g <- which(genes$geneName == gene)
		barWidth <- (genes[g,'end']-genes[g,'start'])/length(cells)
		
		pos <- genes[g,'start'] + (genes[g,'end']-genes[g,'start'])/2 - barWidth*(length(cells)/2)
		
		if(scaleAll == TRUE){
			exp <- (set1[gene,cells]/me)*max(yLim)
		}else{
			exp <- scaleData(set1[gene,cells],0,max(yLim),1)
		}
		
		for(c in cells){
			x <- pos + barWidth*(which(cells == c)-1)
			y <- as.numeric(as.character(exp[which(cells == c)]))
# 			rect(x, 0, x+barWidth, y, col=rgb(0.66,0.66,0.99, alpha=0.4), border = 1)
			rect(x, 0, x+barWidth, y, col=which(cells == c), border = 0)
			if(y < 1){segments(x,0,x+barWidth,0,lwd=1.75, col=1)}
		}
		
	}
# 	axis(2)
	axis(2, at=c(0,max(yLim,na.rm=TRUE)), labels=c(0,signif(max(set1,na.rm=TRUE),2)), cex.axis=cex, las=3)
# 	print9)
}

plotAxis <- function(intrv, orientation='DOWN', vertical=FALSE, cex=1.2){
	
	plotLim <- c(intrv$start,intrv$end)
	chr <- as.character(intrv$chrom)
	
	print(plotLim)
	
	if(vertical == TRUE)
	{
# 		par(mar = c(plotMar[2], 0, plotMar[1], 0))
		xlim=c(0,0)
		plot(0, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F, ylim=plotLim, xlim=xlim)
		axis(2, cex.axis=cex)
# 		text(x=0,y=plotLim[2], labels=chr, pos=3, cex=1.4, col='grey15')
	}else{
		if(orientation == 'DOWN'){orientation=1; axMar <- c(2.3,0.3)
		}else{orientation=3; axMar <- c(0.3,2.3)}

# 		par(mar = c(axMar[1], plotMar[1], axMar[2], plotMar[2]))
		yLim=c(0,100)
		plot(0, xlim=plotLim, ylim=yLim, xlab='', ylab='', yaxt='n', xaxt='n', lwd=1.5, frame=F)
		axis(orientation, cex.axis=cex)
# 		axis(2, labels=chr, at=0, las=2, tick=FALSE, cex.axis=1.6)
		text(chr, x=plotLim[1] + (plotLim[2] - plotLim[1])/2, y=100, pos=1, cex=cex)
	}
}





### Insulation
gtrack.2d.gen_insu_track = function(track_nm, scale, res, min_diag_d=1000, new_track, description=""){
	
	k_reg = 10
	if (length(track_nm) > 1){
		prof = gtrack.2d.gen_joint_insu_prof(track_nm, scale, res, min_diag_d)
	}else{
		print("single source track")
		prof = gtrack.2d.gen_insu_prof(track_nm, scale, res, min_diag_d)
	}
	
	message("names ", paste(names(prof), collapse=","))
	names(prof)[1] = "chrom"
	names(prof)[2] = "start"
	names(prof)[3] = "end"
	
	gtrack.create_sparse(track=new_track, prof[,c(1,2,3)], log2(prof$obs_ins/(prof$obs_big+k_reg)), description=description)
}

gtrack.2d.gen_joint_insu_prof = function(track_nms, scale, res, min_diag_d = 1000){
	
	ins = gtrack.2d.gen_insu_prof(track_nms[1], scale, res, min_diag_d)
	
	for (track in track_nms[-1]){
		t_ins =  gtrack.2d.gen_insu_prof(track, scale, res, min_diag_d)
		ins$obs_big = rowSums(cbind(ins$obs_big, t_ins$obs_big), na.rm=T)
		ins$obs_ins = rowSums(cbind(ins$obs_ins, t_ins$obs_ins), na.rm=T)
	}
	
	ins$obs_big[ins$obs_big == 0] = NA
	ins$obs_ins[ins$obs_ins == 0] = NA
	return(ins)
}

gtrack.2d.gen_insu_prof = function(track_nm, scale, res, min_diag_d=1000){
	
	iter_1d = giterator.intervals(intervals=ALLGENOME[[1]], iterator=res)
	
	iter_2d = gintervals.2d(chroms1 = iter_1d$chrom, starts1=iter_1d$start, ends1=iter_1d$end,
				chroms2 = iter_1d$chrom, starts2=iter_1d$start, ends2=iter_1d$end)
	
	if(length(gvtrack.ls("obs_big")) == 1){gvtrack.rm("obs_big")}
	
	if(length(gvtrack.ls("obs_ins")) == 1){gvtrack.rm("obs_ins")}
	
	gvtrack.create("obs_big", track_nm, "weighted.sum")
	gvtrack.create("obs_ins", track_nm, "weighted.sum")
	gvtrack.iterator.2d("obs_big", 
		eshift1=scale, sshift1=-scale, 
		eshift2=scale, sshift2=-scale)
	gvtrack.iterator.2d("obs_ins", 
		eshift1=0, sshift1=-scale, 
		eshift2=scale, sshift2=0)
	
	message("will iter on ", dim(iter_2d)[1])
	ins = gextract("obs_big", "obs_ins", gintervals.2d.all(), iterator=iter_2d, band=c(-scale*2,0))
	ins_diag = gextract("obs_big", "obs_ins", gintervals.2d.all(), iterator=iter_2d, band=c(-min_diag_d,0))
	
	###fixing NA bug
	ins[is.na(ins)] = 0
	ins_diag[is.na(ins_diag)] = 0
	##
	ins$obs_big = ins$obs_big - ins_diag$obs_big
	ins$obs_ins = ins$obs_ins - ins_diag$obs_ins
	message("will retrun ins with ", dim(ins)[1], " rows")
	return(ins)
}

### Borders / Domains
gtrack.2d.get_insu_doms = function(insu_track, thresh, iterator=500){
	doms = gscreen(sprintf("is.na(%s) | %s > %f", insu_track, insu_track, thresh), iterator=iterator)
	return(doms)
}

gtrack.2d.get_insu_borders = function(insu_track, thresh, iterator=500){
	bords = gscreen(sprintf("!is.na(%s) & %s < %f", insu_track, insu_track, thresh), iterator=iterator)
	return(bords)
}




