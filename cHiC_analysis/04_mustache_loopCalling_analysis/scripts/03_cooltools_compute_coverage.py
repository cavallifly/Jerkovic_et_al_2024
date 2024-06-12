import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import cooltools
import cooler


mainDir="./"



#for sample in ["chic_mESC_DAB_merge_mm10_IJ","chic_mESC_DA_merge_mm10_IJ","chic_mESC_DB_merge_mm10_IJ","chic_mESC_dCTCF_merge_mm10_IJ","chic_mESC_Dprom_merge_mm10_IJ","chic_mESC_polyA_merge_mm10_IJ","chic_mESC_UPPax6_merge_mm10_IJ","chic_mESC_WT_merge_mm10_IJ","chic_NPC_DAB_merge_mm10_IJ","chic_NPC_DA_merge_mm10_IJ","chic_NPC_DB_merge_mm10_IJ","chic_NPC_dCTCF_merge_mm10_IJ","chic_NPC_Dprom_merge_mm10_IJ","chic_NPC_polyA_merge_mm10_IJ","chic_NPC_WT_merge_mm10_IJ"]:
for sample in ["chic_mESC_DAB_merge_onchr18_mm10_IJ_20kb","chic_mESC_DAB_merge_onchr18CaptureDS_mm10_IJ_20kb"]: #,"chic_mESC_DA_merge_onchr18_mm10_IJ_20kb","chic_mESC_DB_merge_onchr18_mm10_IJ_20kb","chic_mESC_dCTCF_merge_onchr18_mm10_IJ_20kb","chic_mESC_Dprom_merge_onchr18_mm10_IJ_20kb","chic_mESC_polyA_merge_onchr18_mm10_IJ_20kb","chic_mESC_UPPax6_merge_onchr18_mm10_IJ_20kb","chic_mESC_WT_merge_onchr18_mm10_IJ_20kb","chic_NPC_DAB_merge_onchr18_mm10_IJ_20kb","chic_NPC_DA_merge_onchr18_mm10_IJ_20kb","chic_NPC_DB_merge_onchr18_mm10_IJ_20kb","chic_NPC_dCTCF_merge_onchr18_mm10_IJ_20kb","chic_NPC_Dprom_merge_onchr18_mm10_IJ_20kb","chic_NPC_polyA_merge_onchr18_mm10_IJ_20kb","chic_NPC_WT_merge_onchr18_mm10_IJ_20kb"]:
    #for resolution in [5000,10000,20000]:
    for resolution in [40000,100000]:        
        for coolFile,outFile in [('./%s.mcool' % (sample),'MADmax_badColumns_%s_at_%dbp_resolution.txt' % (sample, resolution))]:
            
            if os.path.exists(mainDir + outFile):
                print(outFile)
                continue

            print("%d %s %s" % (resolution, coolFile, outFile))
            fpOut=open(mainDir + outFile,"w")

            # to print which resolutions are stored in the mcool, use list_coolers
	    #print(cooler.fileops.list_coolers(mainDir + coolFile))

            print("Loading %s at %dbp resolution" % (coolFile,resolution))
            clr = cooler.Cooler(mainDir + coolFile + '::resolutions/%d' % resolution)
            chromstarts = []
            #for i in clr.chromnames:
            for i in ['chr18']:
                print(f'{i} : {clr.extent(i)}')
                bin1 = clr.extent(i)[0] + int(52580000/resolution)
                bin2 = clr.extent(i)[0] + int(57600000/resolution) 
                for bin in clr.bins()[bin1:bin2].itertuples():
                    #print(bin['chrom','start','end','weight'])
                    #print(bin[1],bin[2],bin[3],bin[4],bin[5])
                    if(np.isnan(bin[4])):
                        print(bin[1:5],file=fpOut)
                    #print(bin[1:2])

	    #clr.bins().fetch(115540)
            #exit(1)
            #region = ('18',52760000,57280000)
            #clr = clr.chromnames.fetch('18')
            #print(clr)

            print("Computing cis_coverage and tot_coverage")
            cis_coverage, tot_coverage = cooltools.coverage(clr)
            for i in range(len(tot_coverage)):
                print(clr.bins()[i],cis_coverage[i],tot_coverage[i], file=fpOut)
            continue
