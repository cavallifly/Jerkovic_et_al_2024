library(doParallel)
library(misha)
library(shaman)
options(scipen=20,gmax.data.size=1e9,shaman.mc_support=1)
args = commandArgs(trailingOnly=TRUE)

mDBloc <- '/zdata/data/mishaDB/trackdb/'
#mDBloc <- '/media/data/mishaDB/trackdb/'
db <- 'mm10'
dbDir <- paste0(mDBloc,db,'/')

gsetroot(dbDir)
gdb.reload()

# Create shuffled track
sample = args[1]
currentTracks <- paste0('chic.',sample)
currentWD   <- '.'
grdSmall    <- 500e3
grdHigh     <- 1000e3
outName <- paste0(currentTracks,'_shuffle_',grdSmall/1e3,'Small_',grdHigh/1e3,'High')

if(dir.exists(outName)){next;}


print(max(gintervals.all()$end))
print(currentTracks)
print(outName)
print(grdSmall)
print(grdHigh)

for (currentTrack in currentTracks)
{
    shaman_shuffle_hic_track(track_db=dbDir, obs_track_nm=currentTrack, exp_track_nm=outName, work_dir=currentWD, max_jobs=4, grid_small=grdSmall, grid_high=grdHigh)
}
