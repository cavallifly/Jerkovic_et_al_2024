library(strawr)
library(hicrep)
options("scipen"=999, max.print=1000000, pillar.colnames = "full")

#HiCRep has several user-adjustable parameters:

# resol is the resolution of the Hi-C data, i.e. bin size. For example,
# resol= 40000 for data with 40kb resolution. lbr and ubr are the lower
# and upper bounds of the genomic distance between interaction
# loci. Using these bounds, users can limit the reproducibility
# assessment to the interactions within a given range of genomic
# distance. In this example, the interactions whose interacting loci
# are no more than 5M apart are considered.

# h is the smoothing parameter that controls the smoothing level. A
# larger value of h leads to a higher level of smoothing. When h=0, no
# smoothing is applied. The optimal parameter choice can be obtained by
# running the function htrain() on a pair of reasonably deeply
# sequenced data. Users may use these values for their analyses if
# training is not feasible. Below is an example code of htrain() using
# training samples mat1 and mat2. It will select the best h in the
# range of 0-10.

# h_value <- htrain(mat1, mat2, resol = 40000, lbr = 0, ubr = 5000000,
# range = 0:10) Training the smoothing parameter h will increase the
# computational time. However, it is not necessary to train h every
# time. In general, for a given resolution, the value of h trained from
# a pair of deeply sequenced biological replicates can be used for other
# datasets with the same resolution. Here we provide the h values
# trained from two hESC replicates in Dixon et al 2015 (GEO accession:
# GSE52457), which were sequenced at 500 million reads, for a range of
# resolutions. Users may directly use these values.

#Resolution ----- h
#10kb ----- 20
#25kb ----- 10
#40kb ----- 5
#100kb ----- 3
#500kb ----- 1 or 2
#1Mb ----- 0 or 1

# For a given HiC dataset, a higher resolution matrix usually requires
# more smoothing, i.e. a higher h value, to enhance its domain
# structures, due to the increasing level of sparsity in the data. To
# compare reproducibility between samples with the same resolution, the
# same smoothing parameter should be used.

inDir='./'

resolutions = c(5000,10000) #,20000)
#lbrs        = c(0,1000000,2000000,3000000,4000000)
#ubrs        = c(0,1000000,2000000,3000000,4000000)
lbrs        = c(0)
ubrs        = c(4000000)

### Code to train h value ###
for(resolution in resolutions)
{
    outFile=paste0("01_hicRep_analysis_captureRegion_get_smoothing_parameter_res",resolution,"bp_hicFiles.out")
    for(lbr in lbrs)
    {
        for(ubr in ubrs)
        {
            if(ubr<lbr+1000000){next}
            chr = '18'

            hic1=list.files("../01_hic_files/",pattern=paste0(".*NPC_WT_.*ep1.*hic"),full.names=T,recursive=T)
            hic2=list.files("../01_hic_files/",pattern=paste0(".*NPC_WT_.*ep2.*hic"),full.names=T,recursive=T)	    

            mat1 <- hic2mat(hic1, chromosome1 = chr, chromosome2 = chr, resol = resolution, method = "NONE")
            mat2 <- hic2mat(hic2, chromosome1 = chr, chromosome2 = chr, resol = resolution, method = "NONE")

            h_value <- htrain(mat1, mat2, resol = resolution, lbr = lbr, ubr = ubr, range = 0:30)
	    line = paste(lbr,ubr,resolution,h_value)
	    write.table(line, outFile, append = TRUE, sep = "\t", dec = ".", row.names = F, col.names = F)
        }
    }
}
quit()
### Code to train h value ###
