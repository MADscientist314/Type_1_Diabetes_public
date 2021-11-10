
sam<-data.frame(sample_data(ASV_physeq_core))
sam<-distinct(select(data.frame(sam),c(patient, disease,SampleType,Antibiotics_2)))
euc_dist_2 <- dist(t(vst_trans_count_tab[ , colnames(vst_trans_count_tab) %in% rownames(sam)]))

# and now making a sample info table with just the basalts
sample_info_tab_2 <- sample_info_tab[row.names(sample_info_tab) %in% rownames(sam), ]
sample_info_tab_2
# running betadisper on just these based on level of alteration as shown in the images above:
beta<-betadisper(euc_dist_2, sample_info_tab_2$Antibiotics_2) # 0.7
permutest(beta)
capture.output(permutest(beta),file = "euc_vst_energ_bialko_proc_2_within_group_permutest_betadisper.tsv")


b<-adonis(formula = euc_dist_2~disease*Antibiotics_2,
          data =  sample_info_tab_2,method = "bray",
          strata = sample_info_tab_2$SampleType, permutations = 9999) # 0.003
b
capture.output(b, file="euc_dist_2_adonis_disease_Antibiotics_2.txt")






sam<-data.frame(sample_data(Delivery_physeq_core))
sam<-distinct(select(data.frame(sam),c(patient, disease,SampleType,Delivery)))
euc_dist_2 <- dist(t(vst_trans_count_tab[ , colnames(vst_trans_count_tab) %in% rownames(sam)]))

# and now making a sample info table with just the basalts
sample_info_tab_2 <- sample_info_tab[row.names(sample_info_tab) %in% rownames(sam), ]
sample_info_tab_2
# running betadisper on just these based on level of alteration as shown in the images above:
beta<-betadisper(euc_dist_2, sample_info_tab_2$Delivery) # 0.7
permutest(beta)



bdeliv<-adonis(formula = euc_dist_2~disease*Delivery,
          data =  sample_info_tab_2,method = "bray",
          strata = sample_info_tab_2$SampleType, permutations = 9999) # 0.003
bdeliv
capture.output(b, file="euc_dist_2_adonis_disease_Delivery.txt")