bedtools merge -d 300 -i <(cat NPC_wt_Pax6_peaks_q0.01_summits.bed ESC_P_Pax6_peaks_q0.01_summits.bed | sort -k 1,1d -k 2,2n -k 3,3n) > merged_and_filtered_Pax6_peaks_q0.01_summits_bedtools.bed
