#!/bin/env Rscript

options(warn=-1)

suppressWarnings(library(misha))
mDBloc <-  '../mishaDB/trackdb/'
db <- 'mm10'
dbDir <- paste0(mDBloc,db,'/')
gdb.init(dbDir)
gdb.reload()

cellType <- "XXXcellTypeXXX"
print(cellType)

hicTracks <- c((paste0('chic.models_',cellType,'.',gsub('.track','',list.files(paste0(dbDir,paste0('tracks/chic/models_',cellType,'/')))))))
obsTracks <- hicTracks[grep("_XXXdirXXX_",hicTracks)]
print(obsTracks)

### BEGIN Insulation AuxFunctions ###
gtrack.2d.gen_insu_track = function(track_nm, scale, res, min_diag_d=1000, new_track, description="")
{
    k_reg = 10
    if (length(track_nm) > 1)
    {
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

gtrack.2d.gen_joint_insu_prof = function(track_nms, scale, res, min_diag_d = 1000)
{
    ins = gtrack.2d.gen_insu_prof(track_nms[1], scale, res, min_diag_d)
	
    for (track in track_nms[-1])
    {
	t_ins =  gtrack.2d.gen_insu_prof(track, scale, res, min_diag_d)
	ins$obs_big = rowSums(cbind(ins$obs_big, t_ins$obs_big), na.rm=T)
	ins$obs_ins = rowSums(cbind(ins$obs_ins, t_ins$obs_ins), na.rm=T)
    }
	
    ins$obs_big[ins$obs_big == 0] = NA
    ins$obs_ins[ins$obs_ins == 0] = NA
    return(ins)
}

gtrack.2d.gen_insu_prof = function(track_nm, scale, res, min_diag_d=1000)
{
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
	

    ins[is.na(ins)] = 0
    ins_diag[is.na(ins_diag)] = 0

    ins$obs_big = ins$obs_big - ins_diag$obs_big
    ins$obs_ins = ins$obs_ins - ins_diag$obs_ins
    message("will retrun ins with ", dim(ins)[1], " rows")
    return(ins)
}
### END Insulation AuxFunctions ###

for(k in seq(400, 400, by=200))
{
    insScale <- k * 1e+3
    insResolution <- 5e+3
	
    for(obsTrack in obsTracks)
    {
        insTrack <- obsTrack		
        insTrackName <- paste0('insulation.models_',cellType,'.',gsub(paste0('chic.models_',cellType,'.'),'',obsTrack),'_',insScale/1000,'kb')
        if(gtrack.exists(insTrackName)){print(paste0(insTrackName,' exists')); next;}

        gtrack.2d.gen_insu_track(insTrack, insScale, insResolution, min_diag_d=1000, insTrackName, description="")
    }    
}
