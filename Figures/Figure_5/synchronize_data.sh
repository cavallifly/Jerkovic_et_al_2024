mkdir -p panelA panelB panelE panelF

### Score maps panel A
rsync -avz nikita:/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/03_plots_cHiC_and_ChIPseq/scoreMaps*NPC_dCTCFD* ./panelA/

### Insulation profiles panel B
rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/03_insulation_analysis_cooler/insulation_profiles_cooler/*NPC_WT*bedGraph ./panelB/

rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/03_insulation_analysis_cooler/insulation_profiles_cooler/*NPC_dCTCF*bedGraph ./panelB/

rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/03_insulation_analysis_cooler/insulation_profiles_cooler/*NPC_DA_*bedGraph ./panelB/

### Differential score maps panel F
rsync -avz nikita:/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/02_differential_maps/ShamanScores/diffMaps_chr18_54680000_56000000_r5000bp_k250_kexp250_NPC_dCTCF_vs_NPC_DA_k250_kexp250.tsv ./panelF/

rsync -avz nikita:/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/02_differential_maps/ShamanScores/diffMaps_zmin-100_zmax100/*NPC_dCTCF_vs_NPC_DA* ./panelF/

### ATAC-seq quantification panelE
rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/ChIPseq_quantification/Ivana_ATACseq/meanAndStdDev_normalizedCounts_using_DEseq2_for_ANOVA_wt_dCTCF_with_stats.tsv ./PanelE/

rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/ChIPseq_quantification/Ivana_ATACseq/results_from_DEseq2_analysis.tsv ./PanelE/

rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/ChIPseq_quantification/scripts_MDS/barplot_with_pvalues_from_DEseq2.R ./PanelE/


# Remove not needed files
rm -fvr ./*/*r2kb* ./*/*w200kb* ./*/*~ *~ ./*/score*_dCTCF_* ./*/score*_DA_* ./*/Insul*merge*bedGraph ./*/*chr18_53859000*_56456400bp* ./*/*toModify*
