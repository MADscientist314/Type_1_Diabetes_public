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
setwd( "k:/github/Type_1_Diabetes_public/gamar/")
#set the seed
#register the cluster (Will not work like this on Windows)
registerDoParallel(cores = 48)
##########################################################################
#Import the data
ASV_physeq<-readRDS("ASV_physeq_with_tree_07112022.RDS")
ASV_phyeq
ASV_physeq_core<-readRDS("ASV_physeq_core_with_tree_07112022.RDS")
ASV_physeq_core
#a<-readRDS("../ASV_physeq_core.RDS")
ASV_physeq<-prune_samples(samples = sample_names(a),ASV_physeq)
ASV_physeq_core<-prune_samples(samples = sample_names(a),ASV_physeq_core)

##########################################################################
# Run DESEQ2 VST normalaiztion
deseq_counts<- phyloseq_to_deseq2(physeq = ASV_physeq,design = ~disease)
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")

# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(ASV_physeq)
vst_physeq <- phyloseq(vst_count_phy, tax_table(ASV_physeq),sample_info_tab_phy,phy_tree(ASV_physeq))

##########################################################################
#run the weighted Unifrac
ASV_physeq_wunifrac<-UniFrac(physeq = ASV_physeq,
                             normalized = T,
                             parallel = T,
                             weighted = T)
ASV_physeq_core_wunifrac<-UniFrac(physeq = ASV_physeq_core,
                             normalized = T,
                             parallel = T,
                             weighted = T)
vst_physeq_wunifrac<-UniFrac(physeq = vst_physeq,
                             normalized = T,
                             parallel = T,
                             weighted = T)

#run the unweighted Unifrac
ASV_physeq_unifrac<-UniFrac(physeq = ASV_physeq,
                            parallel = T,
                            weighted = F)
ASV_physeq_core_unifrac<-UniFrac(physeq = ASV_physeq_core,
                            parallel = T,
                            weighted = F)
vst_physeq_unifrac<-UniFrac(physeq = vst_physeq,
                                 parallel = T,
                                 weighted = F)

##########################################################################
# run the adonis premanova analyses
sam<-data.frame(sample_data(ASV_physeq_core))
sam
sam%>%filter(patient==19)
sam<-sam[sample_names(ASV_physeq),]
sam2<-sam%>%group_by(patient,disease)%>%
  select(patient,disease,Delivery,HbA1C_1trym_2)%>%
  filter(!is.na(HbA1C_1trym_2))%>%
  filter(!is.na(Delivery))%>%
  ungroup()%>%
  mutate(across(where(is.character),as.factor))
  
data.frame(sam2)
sam
sam<-sam%>%
  mutate(SampleID=rownames(sam))%>%
  select(SampleID,patient,disease,Host,SampleType)%>%
  mutate(across(where(is.character),as.factor))%>%
  distinct_all()
sam
sam3<-full_join(sam,sam2)
sam3
#   mutate(across(where(is.character),factor))%>%
#   mutate(Delivery=Delivery,
#          HbA1C_1trym_2=HbA1C_1trym_2)%>%distinct_all()
data.frame(sam)
sam
sam<-sam%>%group_by(patient,disease)%>%
  summarise(Delivery=Delivery,
            HbA1C_1trym_2=HbA1C_1trym_2)

sam
sam<-sam%>%group_by(patient,disease)%>%summarise_all(.funs = distinct)
sam<-sam%>%group_by(patient,disease)%>%
  select(patient,disease,Host,SampleType,Delivery,HbA1C_1trym_2)%>%
  mutate(across(where(is.character),factor))%>%
  mutate(Delivery=Delivery,
         HbA1C_1trym_2=HbA1C_1trym_2)%>%ungroup()
sam
adonis2(ASV_physeq_wunifrac~.,data = sam,
        permutations = 9999,
        parallel = 48) # 0.003

adonis2(ASV_physeq_core_wunifrac~disease*SampleType*Host,data = sam,
        permutations = 9999,
        parallel = 48) # 0.003

adonis2(vst_physeq_wunifrac~disease*SampleType*Host,data = sam,
        permutations = 9999,
        parallel = 48) # 0.003


adonis2(ASV_physeq_unifrac~disease*SampleType*Host,data = sam,
        permutations = 9999,
        parallel = 48) # 0.003

adonis2(ASV_physeq_core_unifrac~disease*SampleType*Host,data = sam,
        permutations = 9999,
        parallel = 48) # 0.003

adonis2(vst_physeq_unifrac~disease*SampleType*Host,data = sam,
        permutations = 9999,
        parallel = 48) # 0.003




# run the beta diversity analyses
sam<-data.frame(sample_data(ASV_physeq))

adonis2(vst_physeq_wunifrac~disease*SampleType*Host,data = sam,
        permutations = 9999,
        parallel = 48) # 0.003
my_pal<-pal_aaas()(10)
disease_pal<-c(my_pal[2],my_pal[10])
show_col(disease_pal)
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
  geom_point(size=1) + 
  #coord_fixed(sqrt(vst_wunifrac_pcoa_eigen_vals[2]/vst_wunifrac_pcoa_eigen_vals[1])) + ggtitle("PCoA") + 
  facet_wrap(Host+SampleType~.,nrow = 3,scales = "free")+
  theme_pubclean()+
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



################################################################################
# missing 1st trimester A1c value imputations
# objective: impute missing values for 1st trimester hemologbin A1C
# and make a normal distribution of values for the controls
################################################################################

library(caret)
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

################################################################################




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