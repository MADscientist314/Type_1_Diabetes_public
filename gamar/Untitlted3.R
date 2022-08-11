#load the libraries
library(phyloseq)
library(tidyverse)
library(doParallel)
library(DESeq2)
library(vegan)
library(ggpubr)
library(ggsci)
library(scales)
#set the working directory
setwd( "/media/jochum00/Jochum_Raid/T1d/gamar")
#set the seed
#register the cluster (Will not work like this on Windows)
registerDoParallel(cores = 48)
##########################################################################
ASV_physeq

meta<-data.frame(read.table("meta_pt+disease+Delivery+HBA1c_1trym.tsv",header = T,sep = "\t"))
rownames(meta)<-meta$SampleID



infant<-subset_samples(physeq = ASV_physeq,Host=="Mother")
infant<-subset_taxa(physeq = infant,taxa_sums(infant)>1)
infant




##########################################################################
# Run DESEQ2 VST normalaiztion
deseq_counts<- phyloseq_to_deseq2(physeq = infant,design = ~disease)
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")

# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(infant)
vst_physeq <- phyloseq(vst_count_phy, tax_table(infant),sample_info_tab_phy,phy_tree(infant))




# MAke some plots
# # generating and visualizing the PCoA with phyloseq
# vst_pcoa <- ordinate(vst_physeq, 
#                      method="PCoA",
#                      distance=vst_physeq_wunifrac)
# 
# vst_wunifrac_pcoa_eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
# 
# # and making our new ordination of just basalts with our adonis statistic
# plot_ordination(vst_physeq, vst_pcoa, 
#                 color="disease",
#                 shape="disease") +
#   scale_color_manual(values=disease_pal)+
#   geom_smooth(method = "glm",se = F) + 
#   #stat_ellipse()+
#   # coord_fixed(sqrt(vst_wunifrac_pcoa_eigen_vals[2]/vst_wunifrac_pcoa_eigen_vals[1])) + ggtitle("PCoA") + 
#   facet_wrap(SampleType~.,nrow = 3,scales = "free")+
#   theme_pubr()+
#   #xscale("log10")+
#   #yscale("log10")+
#   theme(legend.position="right")
# 





##########################################################################
#run the weighted Unifrac
ASV_physeq_wunifrac<-UniFrac(physeq = infant,
                             normalized = T,
                             parallel = T,
                             weighted = T)

vst_physeq_wunifrac<-UniFrac(physeq = vst_physeq,
                             normalized = T,
                             parallel = T,
                             weighted = T)

#run the unweighted Unifrac
ASV_physeq_unifrac<-UniFrac(physeq = infant,
                            parallel = T,
                            weighted = F)
vst_physeq_unifrac<-UniFrac(physeq = vst_physeq,
                            parallel = T,
                            weighted = F)







#Make some plots
# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, 
                     method="PCoA",
                     distance=vst_physeq_wunifrac)

vst_wunifrac_pcoa_eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

# and making our new ordination of just basalts with our adonis statistic
plot_ordination(vst_physeq, vst_pcoa, 
                color="disease",
                shape="disease") +
  scale_color_manual(values=disease_pal)+
  geom_smooth(method = "glm",se = F) + 
  #stat_ellipse()+
  # coord_fixed(sqrt(vst_wunifrac_pcoa_eigen_vals[2]/vst_wunifrac_pcoa_eigen_vals[1])) + ggtitle("PCoA") + 
  facet_wrap(SampleType~.,nrow = 3,scales = "free")+
  theme_pubr()+
  #xscale("log10")+
  #yscale("log10")+
  theme(legend.position="right")
















##########################################################################
# run the adonis premanova analyses
sam<-data.frame(sample_data(infant))
##########################################################################
##########################################################################
# I stopped here on [1] "2022-07-20 16:14:34 CDT"

Sys.time()
##########################################################################
##########################################################################
infant
#sam%>%filter(patient==19)
sam<-sam[sample_names(infant),]
# sam
# sam2<-sam%>%group_by(patient,disease)%>%
#   select(patient,disease,Delivery,HbA1C_1trym_2)%>%
#   mutate(across(where(is.character),as.factor))%>%distinct_all()
# sam2  
# data.frame(sam2)
# 
# sam<-sam%>%
#   mutate(SampleID=rownames(sam))%>%
#   select(SampleID,patient,disease,Host,SampleType)%>%
#   mutate(across(where(is.character),as.factor))%>%
#   distinct_all()
# sam
# sam3<-full_join(sam,sam2)
# sam3
# #   mutate(across(where(is.character),factor))%>%
# #   mutate(Delivery=Delivery,
# #          HbA1C_1trym_2=HbA1C_1trym_2)%>%distinct_all()
# data.frame(sam)
# # sam
# # sam<-sam%>%group_by(patient,disease)%>%
# #   summarise(Delivery=Delivery,
# #             HbA1C_1trym_2=HbA1C_1trym_2)
# # 
# # sam
# # sam<-sam%>%group_by(patient,disease)%>%summarise_all(.funs = distinct)
sam<-sam%>%group_by(patient,disease)%>%
  select(patient,disease,SampleType,Delivery,HbA1C_1trym_2)%>%
  mutate(across(where(is.character),factor))%>%
  ungroup()
sam

with(sam, adonis2(ASV_physeq_wunifrac ~ disease, data = sam, permutations = 999, strata = patient))

adonis2(ASV_physeq_wunifrac~.,data = sam,
        permutations = 9999,
        parallel = 48) # 0.003

# adonis2(ASV_physeq_core_wunifrac~disease*SampleType*Host,data = sam,
#         permutations = 9999,
#         parallel = 48) # 0.003
adonis2(vst_physeq_unifrac~disease+HbA1C_1trym_2+Delivery+SampleType+patient,data = sam,
        permutations = 9999,
        parallel = 48) # 0.003


# adonis2(vst_physeq_wunifrac~disease*SampleType*Delivery,data = sam,
#         permutations = 9999,
#         parallel = 48) # 0.003
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = vst_physeq_wunifrac ~ disease * SampleType * Delivery, data = sam, permutations = 9999, parallel = 48)
# Df SumOfSqs      R2       F Pr(>F)    
# disease                       1 0.002702 0.01373  2.4494 0.0156 *  
#   SampleType                    1 0.017218 0.08751 15.6082 0.0001 ***
#   Delivery                      1 0.002601 0.01322  2.3580 0.0176 *  
#   disease:SampleType            1 0.002070 0.01052  1.8762 0.0598 .  
# disease:Delivery              1 0.000591 0.00300  0.5354 0.8674    
# SampleType:Delivery           1 0.003892 0.01978  3.5285 0.0015 ** 
#   disease:SampleType:Delivery   1 0.001107 0.00563  1.0038 0.4078    
# Residual                    151 0.166574 0.84660                   
# Total                       158 0.196755 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# > # adonis2(ASV_physeq_core_wunifrac~disease*SampleType*Host,data = sam,
#   > #         permutations = 9999,
#   > #         parallel = 48) # 0.003
#   > adonis2(vst_physeq_wunifrac~.,data = sam,
#             +         permutations = 9999,
#             +         parallel = 48) # 0.003
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = vst_physeq_wunifrac ~ ., data = sam, permutations = 9999, parallel = 48)
# Df SumOfSqs      R2       F Pr(>F)    
# patient        52 0.071149 0.36161  1.3432 0.0135 *  
#   disease         1 0.003284 0.01669  3.2234 0.0031 ** 
#   SampleType      1 0.016375 0.08323 16.0752 0.0001 ***
#   Delivery        1 0.000779 0.00396  0.7649 0.6083    
# HbA1C_1trym_2   1 0.001262 0.00642  1.2391 0.2626    
# Residual      102 0.103905 0.52809                   
# Total         158 0.196755 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

vst_physeq_Cervix<-subset_samples(SampleType=="Cervix",physeq = vst_physeq)
vst_physeq_wunifrac_Cervix<-UniFrac(physeq = vst_physeq_Cervix,
                             normalized = T,
                             parallel = T,
                             weighted = T)

sam_Cervix<-sam%>%filter(SampleType=="Cervix")
adonis2(vst_physeq_wunifrac_Cervix~disease,
        data = sam_Cervix,
        permutations = 9999,
        parallel = 48) # 0.003



# run the beta diversity analyses
sam<-data.frame(sample_data(ASV_physeq))

adonis2(vst_physeq_wunifrac_Cervix~disease,data = sam,
        permutations = 9999,
        parallel = 48) # 0.003

my_pal<-pal_aaas()(10)
disease_pal<-c(my_pal[10],my_pal[2])
show_col(disease_pal)
# generating and visualizing the PCoA with phyloseq
vst_pcoa_cervix <- ordinate(vst_physeq_Cervix, 
                     method="PCoA",
                     distance=vst_physeq_wunifrac_Cervix)

vst_wunifrac_pcoa_eigen_vals_Cervix <- vst_pcoa_cervix$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

# and making our new ordination of just basalts with our adonis statistic
plot_ordination(vst_physeq_Cervix, vst_pcoa_cervix, 
                color="disease",
                shape="disease") +
  scale_color_manual(values=disease_pal)+
  geom_smooth(method = "glm",se = F) + 
  #stat_ellipse()+
  coord_fixed(sqrt(vst_wunifrac_pcoa_eigen_vals[2]/vst_wunifrac_pcoa_eigen_vals[1])) + ggtitle("PCoA") + 
#  facet_wrap(SampleType~.,nrow = 3,scales = "free")+
  theme_pubr()+
  #xscale("log10")+
  #yscale("log10")+
  theme(legend.position="right")


##########################################################################
# next steps
#TODO: do this for breakdown with each maternal and infant / SampleType
##########################################################################

#########################################################################
# populate the delivery and A1c values accross the var

sam2<-sam%>%
  select(patient,disease,Host,SampleType,Delivery,HbA1C_1trym_2)%>%
  mutate(SampleID=rownames(sam))
sam2<-as_tibble(sam2)# A tibble: 527 × 6

deliv<-sam2%>%
  select(patient,disease,Delivery)%>%
  filter(!is.na(Delivery))%>%
  distinct_all()
deliv# A tibble: 89 × 3

#generate the A1c for the T1d population
A1C_T1d<-sam2%>%
  select(patient,disease,HbA1C_1trym_2)%>%
  filter(!is.na(HbA1C_1trym_2))%>%
  distinct_all()
A1C_T1d# A tibble: 47 × 3

#generate the A1c for the T1d population
A1C_T1d_missing<-sam2%>%
  select(patient,disease,HbA1C_1trym_2)%>%
  filter(disease=="T1D")%>%
  filter(is.na(HbA1C_1trym_2))%>%
  distinct_all()
A1C_T1d_missing# A tibble: 47 × 3

sam$Host
# impute the A1c for the missing T1D samples
mom<-as_tibble(sam)%>%filter(Host=="Mother")
mom
mom <- mom[,colSums(is.na(mom))<nrow(mom)]
mom
deliv
library(caret)
library(doParallel)
cl<-makeCluster(detectCores())
registerDoParallel(cl)
mom<-as_tibble(sam)%>%
  select(patient,
         disease,
         Host,
         tydzien_porodu_2,
         newborn_weight,
         energ_tluszcz_proc_2,
         energ_bialko_proc_2,
         HbA1C_1trym_2,
         HbA1C_2trym_2,
         HBA1C_POROD_2,
         BMI_after_2,
         BMI_before_2,
         insulin_age,
         insulin_how_long)%>%
  filter(Host=="Mother")%>%
  distinct_all()
mom
mom<-mom%>%mutate(across(where(is.character),factor))%>%
  mutate(patient=as.factor(patient))
mom_T1d<-mom%>%filter(disease=="T1D")


knn_preproc<-preProcess(mom_T1d,method = c("knnImpute",
                                           #"bagImpute",
                                           #"medianImpute",
                                           #"zv",
                                           "nzv"#,
                                           #"conditionalX"
))
center_preproc<-preProcess(mom_T1d,method = c("center","scale"),verbose = T,outcome = "disease")

preproc<-predict(center_preproc,mom_T1d_not_missing)
mom_T1d_train<-preproc%>%slice_sample(prop = 0.9,weight_by = patient)
mom_T1d_train
mom_T1d_test<-preproc%>%filter(!patient%in%mom_T1d_train$patient)

mom_T1d_test

knn_preproc<-preProcess(mom_T1d_train,method = c("knnImpute"),verbose = T,outcome = "disease")

mom_T1d_test2<-predict(object = knn_preproc,newdata = mom_T1d_test)

# bag_preproc<-preProcess(mom_T1d,method = c("bagImpute"))
# median_preproc<-preProcess(mom_T1d,method = c("medianImpute"))
# 

# mom_T1d_knn<-predict(object = preproc,newdata = mom_T1d_test)

postResample(obs = mom_T1d_test$HbA1C_1trym_2,pred = mom_T1d_test2$HbA1C_1trym_2)

mom_T1d_bag<-predict(object = preproc,newdata = mom_T1d)
mom_T1d_median<-predict(object = preproc,newdata = mom_T1d)


mom_T1d<-mom_T1d%>%
  mutate(knn=mom_T1d_knn$HbA1C_1trym_2,
         bag=mom_T1d_bag$HbA1C_1trym_2,
         median=mom_T1d_median$HbA1C_1trym_2)

model <- lm(HbA1C_1trym_2~knn+bag+median, 
            data = mom_T1d%>%select(where(is.numeric)))
msummary(model)
#Prediction
lm_fun <- makeFun(model)
mom_T1d_missing<-mom_T1d%>%filter(is.na(HbA1C_1trym_2))
mom_T1d_not_missing<-mom_T1d%>%filter(!is.na(HbA1C_1trym_2))%>%
  mutate(imputed="actual")

mom_T1d_missing<-mom_T1d_missing%>%
  #  group_by(patient,disease)%>%
  rowwise()%>%
  mutate(HbA1C_1trym_2=lm_fun(knn = knn,bag = bag,median = median))%>%
  mutate(imputed="imputed")
summary(knn_preproc)

mom_T1d_2<-full_join(mom_T1d_missing,mom_T1d_not_missing)%>%
  mutate(imputed=factor(imputed))
show_col(pal_aaas()(10))
aaas<-pal_aaas()(10)
my_pal<-c(aaas[6],aaas[10])

p<-ggscatter(data = mom_T1d_2,
             x = "HbA1C_1trym_2",
             y = 0,repel = T,#add.params = list(max.overlaps=10),
             group = "patient",color = "imputed",palette = my_pal)+
  ggrepel::geom_label_repel(label=mom_T1d_2$patient,mapping = aes(color=imputed),
                            max.overlaps = 100)

ggpar(p = p,main = "Distibution of 1st trimester HBA1C values",
      caption = "missing values were imputed using a 
      k nearest neighbors (knn) algorithm (RMSE=0.032,Rsquared=1.0,MAE=0.02)",
      legend = "right")


ggscatter(data = mom_T1d,x = "HbA1C_1trym_2",y = "knn")
ggscatter(data = mom_T1d,x = "HbA1C_1trym_2",y = "bag")
ggscatter(data = mom_T1d,x = "HbA1C_1trym_2",y = "median")




library(vegan)
library(cluster)
mom_T1d_knn
gower_dist <- daisy(mom_T1d_knn,
                    metric = "gower")
gow<-vegdist(x = mom_T1d_knn,method = "gower")
mom_T1d

mom$HbA1C_1trym_2


#generate the A1c for the control population
A1C_control<-sam2%>%
  select(patient,disease,HbA1C_1trym_2)%>%
  filter(is.na(HbA1C_1trym_2))%>%
  distinct_all()%>%
  mutate(HbA1C_1trym_2=round(rnorm(45, 5, 1.2),2))

A1c<-full_join(A1C_T1d,A1C_control)
A1c
# impute A1c values for the control patients
#generate a normal distribution of A1C values around 5 with SD 1.2,
#for the control dataset
A1C
norm <- round(rnorm(47, 5, 1.2),2)
norm
norm <- rnorm(5)
norm
df<-full_join(deliv,A1C)
library(mosaic)
favstats(A1C$HbA1C_1trym_2)
#########################################################################

runme<-function(HOST,ST){
  
  
}