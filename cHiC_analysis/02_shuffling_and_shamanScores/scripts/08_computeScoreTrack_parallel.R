require(misha)
require(shaman)

require(doParallel)

args = commandArgs(trailingOnly=TRUE)

sample = args[1]

#  chmod -R 777 cHiC_scoreTracks/

##########################################
### misha working DB
mDBloc <-  '/zdata/data/mishaDB/trackdb/'
db <- 'mm10'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

source('/zdata/data/auxFunctions/auxFunctions.R')
#options(scipen=20,gmax.data.size=0.5e8,shaman.sge_support=1)
options(scipen=20,gmax.data.size=1e9,shaman.sge_support=1)

step   <- 250000    # Size of the chunk to divide the calculation
#step   <- 1000000    # Size of the chunk to divide the calculation
expand <- step/4 # Size of skin around the main chunk to look for expected counts
#k      <- 100    # Number of neighbours to look for
k      <- 250    # Number of neighbours to look for

for(chr in c("chr18"))
#for(chr in gintervals.all()$chrom)
{
    print(chr)

    #start=52598845
    #end=57599408
    start = 52000000
    end   = 58000000    

    #chrIntrv <- gintervals.all()[gintervals.all()$chrom == chr,]
    chrIntrv <- gintervals(chr,start,end)
    print(chrIntrv)

    chrIter <- giterator.intervals(intervals=g2d(chrIntrv),iterator=c(step,step))
    print(nrow(chrIter))

    ### SCORE MAPS
    obsTrack <- c(paste0("chic.chic_",sample,"_mm10_IJ"))
    #obsTrack <- c("chic.chic_NPC_DAB_Rep1_mm10_IJ","chic.chic_NPC_DAB_Rep2_mm10_IJ")    
    expTrack <- paste0(obsTrack,"_shuffle_500Small_1000High")
    setName  <- paste0("chicScores_",sample,"_mm10_IJ")
    #setName  <- paste0("chicScores_NPC_DAB_merge_mm10_IJ")    
    work_dir <- paste0("./_tmp_",sample,"_",chr,"/")
    #work_dir <- paste0("./_tmp_NPC_DAB_merge_",chr,"/")    

    if(dir.exists(work_dir)){next}
    if(!dir.exists(work_dir)){dir.create(work_dir, mode="7777", recursive=TRUE);}

    print(obsTrack)
    print(setName)
    print(expTrack)

    hicName <- "chic"
	
    trackName <- paste0(setName,"_",chr,"_score_k",k,"_",step/1e3,"kb")
    #trackName <- paste0(setName,"_",chr,"_score_k",k,"_1000kb")    
    outTrack <- paste0(hicName,".",setName,".",trackName)

    print(outTrack)

    if(gtrack.exists(outTrack))
    {
        next
    }

    registerDoParallel(cores=1)

    foreach(i=1:nrow(chrIter)) %dopar% {
        int1 <- chrIter[i,]	
	int2 <- expand2D(int1,expand)[,1:6]
	print(paste0("Interval to score"))
	print(int1)	
	print(paste0("Interval + skin"))	
	print(int2)

	outFile <- paste0(work_dir,chr,"_",i,".scores")
	if(file.exists(outFile))
        {
	    print("done...")
	    return(i)
	}

	genDist <- abs(int1$start2-int1$start1)
	if(genDist <= 500e3)
	{
	    print(paste0("Closer than ",500e3,"bp. Leave it for later!"))	
	    return(i)
        }

	print(paste0("scoring portion ",i,"..."))
	scores <- shaman_score_hic_mat(obs_track_nms=obsTrack,exp_track_nms=expTrack,focus_interval=int1,regional_interval=int2,k=k)
	data <- scores$points
	write.table(data,file=outFile,sep="\t",quote=F,row.names=F,)

	return(i)
    }

    registerDoParallel(cores=1)	

    foreach(i=1:nrow(chrIter)) %dopar% {
        int1 <- chrIter[i,]	
	int2 <- expand2D(int1,expand)[,1:6]
	print(paste0("Interval to score"))
	print(int1)	
	print(paste0("Interval + skin"))	
	print(int2)

	outFile <- paste0(work_dir,chr,"_",i,".scores")
	if(file.exists(outFile))
        {
	    print("done...")
	    return(i)
	}

	genDist <- abs(int1$start2-int1$start1)
	if(genDist > 500e3)
	{
	    print(paste0("Further than ",500e3,"bp. Already done!"))
	    return(i)
        }

	print(paste0("scoring portion ",i,"..."))
	scores <- shaman_score_hic_mat(obs_track_nms=obsTrack,exp_track_nms=expTrack,focus_interval=int1,regional_interval=int2,k=k)
	data <- scores$points
	write.table(data,file=outFile,sep="\t",quote=F,row.names=F,)

	return(i)
    }

    files <- list.files(work_dir, full.names=T,pattern="scores")
	
    trackFolder <- paste0(dbDir,"tracks/",gsub("\\.","/",hicName),"/",setName,"/")
    if(!dir.exists(trackFolder)){dir.create(trackFolder, mode="7777", recursive=TRUE);}

    ### Single Import
    gtrack.2d.import(outTrack, paste("Score track for", hicName," - ",chr), files)
    ### Add Int AND Reverse INT
    # gtrack.2d.import_contacts(outTrack, paste("import diffMaps scores for ", trackName), files)

    ### Remove temp files....
    for (f in files){file.remove(f)}
}
