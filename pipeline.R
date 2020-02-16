############ GEt Bed files of common peaks in K562 for CEBPB with CEBPD/ATF4
cebpb_atf4_K562 <- commonPeaks(target_peak_id="MM1_HSA_K562_CEBPB",
                                motif_only_for_target_peak=TRUE,
                                compared_peak_id="MM1_HSA_K562_ATF4",
                                motif_only_for_compared_peak=TRUE,
                                methylation_profile_in_narrow_region=FALSE)
cebpb_cebpd_K562 <- commonPeaks(target_peak_id="MM1_HSA_K562_CEBPB",
                                 motif_only_for_target_peak=TRUE,
                                 compared_peak_id="MM1_HSA_K562_CEBPD",
                                 motif_only_for_compared_peak=TRUE,
                                 methylation_profile_in_narrow_region=FALSE)

cebpb_atf4_peaks <-cebpb_atf4_K562[,1][[1]]@common_peak
cebpb_atf4_peaks[,2] <- cebpb_atf4_peaks[,2]-100
cebpb_atf4_peaks[,3] <- cebpb_atf4_peaks[,3]+100
write.table(cebpb_atf4_peaks[,1:3],"cebpb_atf4_peaks_K562.bed",sep="\t",quote=F,row.names=F,col.names=F)

cebpb_cebpd_peaks <-cebpb_cebpd_K562[,1][[1]]@common_peak
cebpb_cebpd_peaks[,2] <- cebpb_cebpd_peaks[,2]-100
cebpb_cebpd_peaks[,3] <- cebpb_cebpd_peaks[,3]+100
write.table(cebpb_cebpd_peaks[,1:3],"cebpb_cebpd_peaks_K562.bed",sep="\t",quote=F,row.names=F,col.names=F)
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

################### Peaks of MAFF in K562 compared to HepG2/HeLa-S3

maff_maff_k562 <- commonPeaks(target_peak_id="MM1_HSA_K562_MAFF",
                              motif_only_for_target_peak=TRUE,
                              compared_peak_id="MM1_HSA_HepG2_MAFF",
                              motif_only_for_compared_peak=TRUE,
                              methylation_profile_in_narrow_region=FALSE)

maff_maff_peaks <-maff_maff_k562[,1][[1]]@common_peak
maff_maff_peaks[,2] <- maff_maff_peaks[,2]-100
maff_maff_peaks[,3] <- maff_maff_peaks[,3]+100
write.table(maff_maff_peaks[,1:3],"maff_maff_peaks_K562_HelG2.bed",sep="\t",quote=F,row.names=F,col.names=F)
##
maff_maff_k562 <- commonPeaks(target_peak_id="MM1_HSA_K562_MAFF",
                              motif_only_for_target_peak=TRUE,
                              compared_peak_id="MM1_HSA_HeLa-S3_MAFF",
                              motif_only_for_compared_peak=TRUE,
                              methylation_profile_in_narrow_region=FALSE)

maff_maff_peaks <-maff_maff_k562[,1][[1]]@common_peak
maff_maff_peaks[,2] <- maff_maff_peaks[,2]-100
maff_maff_peaks[,3] <- maff_maff_peaks[,3]+100
write.table(maff_maff_peaks[,1:3],"maff_maff_peaks_K562_HeLa-S3.bed",sep="\t",quote=F,row.names=F,col.names=F)

################### Peaks of MAFF in HepG2 compared to HeLa-S3/K562
maff_maff_HepG2_HeLaS3 <- commonPeaks(target_peak_id="MM1_HSA_HepG2_MAFF",
                                      motif_only_for_target_peak=TRUE,
                                      compared_peak_id="MM1_HSA_HeLa-S3_MAFF",
                                      motif_only_for_compared_peak=TRUE,
                                      methylation_profile_in_narrow_region=FALSE)

maff_HepG2_HeLaS3_peaks <- maff_maff_HepG2_HeLaS3[,1][[1]]@common_peak
maff_HepG2_HeLaS3_peaks[,2] <- maff_HepG2_HeLaS3_peaks[,2]-100
maff_HepG2_HeLaS3_peaks[,3] <- maff_HepG2_HeLaS3_peaks[,3]+100
write.table(maff_HepG2_HeLaS3_peaks[,1:3],"maff_maff_peaks_HepG2_HeLaS3.bed",sep="\t",quote=F,row.names=F,col.names=F)

maff_maff_HepG2_K562 <- commonPeaks(target_peak_id="MM1_HSA_HepG2_MAFF",
                                    motif_only_for_target_peak=TRUE,
                                    compared_peak_id="MM1_HSA_K562_MAFF",
                                    motif_only_for_compared_peak=TRUE,
                                    methylation_profile_in_narrow_region=FALSE)

maff_HepG2_K562_peaks <- maff_maff_HepG2_K562[,1][[1]]@common_peak
maff_HepG2_K562_peaks[,2] <- maff_HepG2_K562_peaks[,2]-100
maff_HepG2_K562_peaks[,3] <- maff_HepG2_K562_peaks[,3]+100
write.table(maff_HepG2_K562_peaks[,1:3],"maff_maff_peaks_HepG2_K562.bed",sep="\t",quote=F,row.names=F,col.names=F)

################### Peaks of MAFF in HeLa-S3 compared to HepG2/K562

maff_maff_HeLaS3_K562 <- commonPeaks(target_peak_id="MM1_HSA_HeLa-S3_MAFF",
                                      motif_only_for_target_peak=TRUE,
                                      compared_peak_id="MM1_HSA_K562_MAFF",
                                      motif_only_for_compared_peak=TRUE,
                                      methylation_profile_in_narrow_region=FALSE)

maff_HeLaS3_K562_peaks <- maff_maff_HeLaS3_K562[,1][[1]]@common_peak
maff_HeLaS3_K562_peaks[,2] <- maff_HeLaS3_K562_peaks[,2]-100
maff_HeLaS3_K562_peaks[,3] <- maff_HeLaS3_K562_peaks[,3]+100
write.table(maff_HeLaS3_K562_peaks[,1:3],"maff_maff_peaks_HeLaS3_K562.bed",sep="\t",quote=F,row.names=F,col.names=F)
##
maff_maff_HeLaS3_HepG2 <- commonPeaks(target_peak_id="MM1_HSA_HeLa-S3_MAFF",
                                      motif_only_for_target_peak=TRUE,
                                      compared_peak_id="MM1_HSA_HepG2_MAFF",
                                      motif_only_for_compared_peak=TRUE,
                                      methylation_profile_in_narrow_region=FALSE)

maff_HeLaS3_HepG2_peaks <- maff_maff_HeLaS3_HepG2[,1][[1]]@common_peak
maff_HeLaS3_HepG2_peaks[,2] <- maff_HeLaS3_HepG2_peaks[,2]-100
maff_HeLaS3_HepG2_peaks[,3] <- maff_HeLaS3_HepG2_peaks[,3]+100
write.table(maff_HeLaS3_HepG2_peaks[,1:3],"maff_maff_peaks_HeLaS3_HepG2.bed",sep="\t",quote=F,row.names=F,col.names=F)
######################################################################### Exporting Target motif's matrix
CEBPBCEBPB_K562 <- TFregulomeR::intersectPeakMatrix(peak_id_x = "MM1_HSA_K562_CEBPB", motif_only_for_id_y = T,motif_only_for_id_x = T, peak_id_y = "MM1_HSA_K562_CEBPB") #intersection of CEBPB-CEBPB

TFregulomeR::exportMMPFM(fun_output = CEBPBCEBPB_K562, 
                         fun = "intersectPeakMatrix", 
                         save_motif_PFM = TRUE, 
                         save_betaScore_matrix = FALSE)
##
CEBPBCEBPB_HepG2 <- TFregulomeR::intersectPeakMatrix(peak_id_x = "MM1_HSA_HepG2_CEBPB", motif_only_for_id_y = T,motif_only_for_id_x = T, peak_id_y = "MM1_HSA_HepG2_CEBPB") #intersection of CEBPB-CEBPB
TFregulomeR::exportMMPFM(fun_output = CEBPBCEBPB_HepG2, 
                         fun = "intersectPeakMatrix", 
                         save_motif_PFM = TRUE, 
                         save_betaScore_matrix = FALSE)
##
MAFFMAFF_K562 <- TFregulomeR::intersectPeakMatrix(peak_id_x = "MM1_HSA_K562_MAFF", motif_only_for_id_y = T,motif_only_for_id_x = T, peak_id_y = "MM1_HSA_K562_MAFF") #intersection of MAFF-MAFF
TFregulomeR::exportMMPFM(fun_output = MAFFMAFF_K562, 
                         fun = "intersectPeakMatrix", 
                         save_motif_PFM = TRUE, 
                         save_betaScore_matrix = FALSE)
##
MAFF_HepG2 <- TFregulomeR::intersectPeakMatrix(peak_id_x = "MM1_HSA_HepG2_MAFF", motif_only_for_id_y = T,motif_only_for_id_x = T, peak_id_y = "MM1_HSA_HepG2_MAFF") 
TFregulomeR::exportMMPFM(fun_output = MAFF_HepG2, 
                         fun = "intersectPeakMatrix", 
                         save_motif_PFM = TRUE, 
                         save_betaScore_matrix = FALSE)
##
MAFF_HeLaS3 <- TFregulomeR::intersectPeakMatrix(peak_id_x = "MM1_HSA_HeLa-S3_MAFF", motif_only_for_id_y = T,motif_only_for_id_x = T, peak_id_y = "MM1_HSA_HeLa-S3_MAFF") #intersection of MAFF-MAFF
TFregulomeR::exportMMPFM(fun_output = MAFF_HeLaS3, 
                         fun = "intersectPeakMatrix", 
                         save_motif_PFM = TRUE, 
                         save_betaScore_matrix = FALSE)
######################################## FPWMs exported
#case 1
source("./R/createFPWM.R")
fpwm_CEBPB_K562_F5 <- createFPWM(mainTF ="CEBPB",
                    partners = c("ATF4","ATF7","ATF3","JUND","FOS","CEBPD"),
                    cell = "K562", 
                    forkPosition = 5)
source("./R/plotFPWM.R")
plotFPWM(fpwm_CEBPB_K562_F5,pdfName="fpwm_CEBPB_K562_F5.pdf")
source("./R/write.FPWM.R")
write.FPWM(FPWM = fpwm_CEBPB_K562_F5,fileName="fpwm_CEBPB_K562_F5.transfac")
write.FPWM(FPWM = fpwm_CEBPB_K562_F5,fileName="fpwm_CEBPB_K562_F5.FPWMtransfac",format="FPWMtransfac")
## case 2
source("./R/createFPWM.R")
fpwm_CEBPB_HepG2_F4 <- createFPWM(mainTF ="CEBPB",
                    partners = c("ATF4","ATF7","ATF3","JUND","FOS","CEBPD"),
                    cell = "HepG2", 
                    forkPosition = 4)
source("./R/plotFPWM.R")
plotFPWM(fpwm_CEBPB_HepG2_F4,pdfName="fpwm_CEBPB_HepG2_F4_plot.pdf")
source("./R/write.FPWM.R")
write.FPWM(FPWM = fpwm_CEBPB_K562_F4,fileName="fpwm_CEBPB_K562_F4.transfac")
write.FPWM(FPWM = fpwm_CEBPB_K562_F4,fileName="fpwm_CEBPB_K562_F4.FPWMtransfac",format="FPWMtransfac")
## case 3

source("./R/createFPWM.R")
fpwm_MAFF_K562_F3_flipped <- createFPWM(mainTF ="MAFF",
                    partners = c("MAFG","C11orf30","NFE2","MAFK","ZNF316","NFE2l2"),
                    cell = "K562", 
                    forkPosition = 3,flipMatrix = TRUE)
source("./R/plotFPWM.R")
plotFPWM(fpwm_MAFF_K562_F3_flipped,pdfName="fpwm_MAFF_K562_F3_Flipped_plot.pdf")
source("./R/write.FPWM.R")
write.FPWM(FPWM = fpwm_MAFF_K562_F3_flipped,fileName="fpwm_MAFF_K562_F3_flipped.transfac")
write.FPWM(FPWM = fpwm_MAFF_K562_F3_flipped,fileName="fpwm_MAFF_K562_F3_flipped.FPWMtransfac",format="FPWMtransfac")
## case 4
