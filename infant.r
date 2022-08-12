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
setwd( "k:/github/Type_1_Diabetes_public/gamar")
#set the seed
#register the cluster (Will not work like this on Windows)
registerDoParallel(cores = detectCores())
##########################################################################
#Import the data
ASV_physeq<-readRDS("ASV_physeq_with_tree_07112022.RDS")
ASV_physeq_core<-readRDS("ASV_physeq_core_with_tree_07112022.RDS")
#Import the updated metadata
meta<-data.frame(read.table("meta_pt+disease+Delivery+HBA1c_1trym.tsv",header = T,sep = "\t"))
#annotate the rownames
rownames(meta)<-meta$SampleID

#prune the samples
ASV_physeq<-prune_samples(samples = as.character(meta$SampleID),x = ASV_physeq)
ASV_physeq_core<-prune_samples(samples = as.character(meta$SampleID),x = ASV_physeq_core)


#overwrite the sampel data
sample_data(ASV_physeq)<-sample_data(meta)
sample_data(ASV_physeq_core)<-sample_data(meta)

meta
################################################################################
######################### Neonatal samples #####################################
################################################################################
infant_adonis<-function(physeq)
  {
  #subset the samples to only infant
  infant<-subset_samples(physeq = physeq,Host=="Infant")
  #remove any taxa with zero counts
  infant<-subset_taxa(physeq = infant,taxa_sums(infant)>1)
  
}






##########################################################################
# next steps
#TODO: do this for breakdown with each maternal and infant / SampleType
##########################################################################

#########################################################################
# populate the delivery and A1c values accross the var
sam<-meta
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
mom<-mom[,colSums(is.na(mom))<nrow(mom)]
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
mom_T1d_2
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

#### Impute control A1c values
A1C_control<-sam2%>%
  select(patient,disease,HbA1C_1trym_2)%>%
  filter(disease=="Control")%>%
  #  filter(is.na(HbA1C_1trym_2))%>%
  distinct_all()%>%
  mutate(HbA1C_1trym_2=round(rnorm(42, 5, 0.6),2),
         imputed="imputed")%>%
  mutate(patient=factor(patient))

range(A1C_control$HbA1C_1trym_2)

A1c<-full_join(mom_T1d_2,A1C_control)

################################################################################



m<-data.frame(read.table("../metadata_07202022.txt",header = T,sep = "\t",row.names = 1))%>%
  select(patient,disease,Host,SampleType,Delivery_Mode_Vaginal_or_Caesarian,HbA1C_1trym_2)%>%
#  mutate(SampleID=rownames(m))%>%
  distinct_all()
m$SampleID<-rownames(m)
n<-as_tibble(m)%>%select(patient,disease,Delivery_Mode_Vaginal_or_Caesarian)%>%distinct_all()%>%mutate(patient=factor(patient))
n
dim(n)
A1c
o<-full_join(A1c,n)
o
p<-o%>%select(patient,disease,Delivery_Mode_Vaginal_or_Caesarian,HbA1C_1trym_2)%>%distinct_all()%>%
  mutate(Delivery=gsub("vaginal","Vaginal",Delivery_Mode_Vaginal_or_Caesarian),
         Delivery=gsub("caesarian","Csection",Delivery),
         Delivery=factor(Delivery))%>%
  select(-Delivery_Mode_Vaginal_or_Caesarian)%>%
  mutate(across(where(is.character),factor))%>%
  distinct_all()%>%
  as_tibble()
p
dim(sam)
sam3<-sam%>%
  select(patient,disease,Host,SampleType)%>%
  mutate(SampleID=rownames(sam),
         patient=factor(patient))%>%
  distinct_all()%>%
  mutate(across(where(is.character),factor))%>%
  as_tibble()
sam3

p%>%filter(disease=="T1D")%>%filter(patient%in%c("20","49","50"))
sam4<-inner_join(sam3,p)%>%distinct_all()%>%filter(!duplicated(SampleID))
sam4%>%arrange(SampleID)


write.table(sam4,"meta_pt+disease+Delivery+HBA1c_1trym.tsv",sep = "\t",row.names = F)

###############################################################################
sample_info_tab<-data.frame(sam4)
rownames(sample_info_tab)<-sample_info_tab$SampleID
sample_data(ASV_physeq_core)<-sample_data(sample_info_tab)
sample_data(ASV_physeq)<-sample_data(sample_info_tab)
################################################################################

