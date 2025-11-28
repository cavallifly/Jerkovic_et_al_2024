mkdir -p panelA panelB panelC

### Score maps panelA
rsync -avz nikita:/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/03_plots_cHiC_and_ChIPseq/scoreMaps*NPC_WT* ./panelA/

rsync -avz nikita:/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/03_plots_cHiC_and_ChIPseq/scoreMaps*NPC_DAD* ./panelA/

### Insulation profiles panelB
rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/03_insulation_analysis_cooler/insulation_profiles_cooler/*NPC_WT*bedGraph ./panelB/

rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/03_insulation_analysis_cooler/insulation_profiles_cooler/*NPC_DA_*bedGraph ./panelB/

### Insulation profiles panelC
rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/qPCR_statistical_analysis/ZFP608_ddCt_values_dC2-C3-AC2C3_for_anova/*.pdf ./panelE/

rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/qPCR_statistical_analysis/ZFP608_ddCt_values_dC2-C3-AC2C3_for_anova/*.bed ./panelE/

rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/qPCR_statistical_analysis/ZFP608_ddCt_values_dC2-C3-AC2C3_for_anova/*.out ./panelE/

# Remove not needed files
rm -fvr ./*/*r2kb* ./*/*w200kb* ./*/*~ *~ ./*/score*_WT_* ./*/score*_polyA_* ./*/*_NPC_WTDStomESCWT_* ./*/Insul*merge*bedGraph ./*/*chr18_53859000*_56456400bp* ./*/*toModify* ./*/*ex?.pdf ./*/*ex?-?.pdf ./*/qq*.pdf
