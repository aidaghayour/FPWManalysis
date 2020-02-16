############ GEt Bed files of common peaks in HepG2 for CEBPB with CEBPD/ATF4

cebpb_atf4_hepg2 <- commonPeaks(target_peak_id="MM1_HSA_HepG2_CEBPB",
                                       motif_only_for_target_peak=TRUE,
                                       compared_peak_id="MM1_HSA_HepG2_ATF4",
                                       motif_only_for_compared_peak=TRUE,
                                       methylation_profile_in_narrow_region=FALSE)
cebpb_cebpd_hepg2 <- commonPeaks(target_peak_id="MM1_HSA_HepG2_CEBPB",
                                motif_only_for_target_peak=TRUE,
                                compared_peak_id="MM1_HSA_HepG2_CEBPD",
                                motif_only_for_compared_peak=TRUE,
                                methylation_profile_in_narrow_region=FALSE)
                                
cebpb_atf4_peaks <-cebpb_atf4_hepg2[,1][[1]]@common_peak
cebpb_atf4_peaks[,2] <- cebpb_atf4_peaks[,2]-100
cebpb_atf4_peaks[,3] <- cebpb_atf4_peaks[,3]+100
write.table(cebpb_atf4_peaks[,1:3],"cebpb_atf4_peaks_HepG2.bed",sep="\t",quote=F,row.names=F,col.names=F)

cebpb_cebpd_peaks <-cebpb_cebpd_hepg2[,1][[1]]@common_peak
cebpb_cebpd_peaks[,2] <- cebpb_cebpd_peaks[,2]-100
cebpb_cebpd_peaks[,3] <- cebpb_cebpd_peaks[,3]+100
write.table(cebpb_cebpd_peaks[,1:3],"cebpb_cebpd_peaks_HepG2.bed",sep="\t",quote=F,row.names=F,col.names=F)

################### Peaks in K562 for MAFF with MAFG/NFE2


maff_nfe2_k562 <- commonPeaks(target_peak_id="MM1_HSA_K562_MAFF",
                                motif_only_for_target_peak=TRUE,
                                compared_peak_id="MM1_HSA_K562_NFE2",
                                motif_only_for_compared_peak=TRUE,
                                methylation_profile_in_narrow_region=FALSE)
maff_mafg_k562 <- commonPeaks(target_peak_id="MM1_HSA_K562_MAFF",
                                 motif_only_for_target_peak=TRUE,
                                 compared_peak_id="MM1_HSA_K562_MAFG",
                                 motif_only_for_compared_peak=TRUE,
                                 methylation_profile_in_narrow_region=FALSE)

maff_nfe2_peaks <-maff_nfe2_k562[,1][[1]]@common_peak
maff_nfe2_peaks[,2] <- maff_nfe2_peaks[,2]-100
maff_nfe2_peaks[,3] <- maff_nfe2_peaks[,3]+100
write.table(maff_nfe2_peaks[,1:3],"maff_nfe2_peaks_K562.bed",sep="\t",quote=F,row.names=F,col.names=F)

maff_mafg_peaks <-maff_mafg_k562[,1][[1]]@common_peak
maff_mafg_peaks[,2] <- maff_mafg_peaks[,2]-100
maff_mafg_peaks[,3] <- maff_mafg_peaks[,3]+100
write.table(maff_mafg_peaks[,1:3],"maff_mafg_peaks_K562.bed",sep="\t",quote=F,row.names=F,col.names=F)
