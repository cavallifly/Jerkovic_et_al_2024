### Score maps
rsync -avz nikita:/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/03_plots_cHiC_and_ChIPseq/scoreMaps*NPC_WT* .

rsync -avz nikita:/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/03_plots_cHiC_and_ChIPseq/scoreMaps*NPC_polyA* .

### Insulation profiles
rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/03_insulation_analysis_cooler/insulation_profiles_cooler/*NPC_WT*bedGraph .

rsync -avz /run/user/1001/gvfs/smb-share\:server\=archive.igh.internal\,share\=cavalli/commun/Marco/2021_04_13_Project_with_Ivana/03_insulation_analysis_cooler/insulation_profiles_cooler/*NPC_polyA*bedGraph .

# Remove not needed files
rm -fr *r2kb* *w200kb* *~ score*_WT_* score*_polyA_* *_NPC_WTDStomESCWT_* Insul*merge*bedGraph *chr18_53859000*_56456400bp* *toModify*
