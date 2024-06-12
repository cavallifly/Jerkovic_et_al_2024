options(warn=-1)

library(misha)

mDBloc <- '../mishaDB/trackdb/'
db <- 'mm10'
dbDir <- paste0(mDBloc,db,'/')
print(mDBloc)
gsetroot(dbDir)
gdb.reload()

dataDir <- 'virtual_cHiC_matrices/'
cellType <- "XXXcellTypeXXX"
flPattern <- c("*_P_*")[1]
print(dataDir)
print(flPattern)

dataTracks <- list.files(dataDir, pattern=glob2rx(flPattern))

for (track in dataTracks)
{
    trackName=gsub('\\.txt$','', track)
    trackName=gsub("\\.","_", trackName)

    if(gtrack.exists(paste0('chic.models_',cellType,'.',trackName))){print(paste0(trackName," exists")); next;}

    gtrack.2d.import(track=paste0('chic.models_',cellType,'.',trackName), description="Models contacts for chr18:53700000-56700000", file=c(paste0('virtual_cHiC_matrices/',track)))
}
quit()
