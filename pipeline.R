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



