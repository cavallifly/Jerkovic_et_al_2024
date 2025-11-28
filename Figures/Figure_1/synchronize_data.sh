mkdir -p panelA panelB

### Score maps
rsync -avz nikita:/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/03_plots_cHiC_and_ChIPseq/scoreMaps*mESC_WT* ./panelA/

rsync -avz nikita:/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/03_plots_cHiC_and_ChIPseq/scoreMaps*NPC_WT* ./panelA/

### Insulation profiles
rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/03_insulation_analysis_cooler/insulation_profiles_cooler/*mESC_WT*bedGraph ./panelB/

rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/03_insulation_analysis_cooler/insulation_profiles_cooler/*NPC_WT*bedGraph ./panelB/


rm ./*/*r2kb* ./*/*w200kb* ./*/*chr18_54680000*_56000000bp* ./*/*~ *~ ./*/score*_WT_* ./*/*_NPC_WTDStomESCWT_* ./*/Insul*merge*bedGraph
