library(phyloseq)
library(vegan)
library(ggpubr)
library(ggsci)
library(mosaic)
library(microbiome)
library(tidyverse)
library(DESeq2)
library(pairwiseAdonis)
# The objective of this script is to get adonis permanova scores 
# for the facet plots in Supplemental figure 9.

# Objective
# I need to code a way where the neonatal and maternal dyad samples all have filled in  delivery, HbA1c>7.2, and disease variables.


##################### Day 00 ###########################

deseq_counts <- phyloseq_to_deseq2(ASV_physeq_core,design = ~disease)
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# now followed by the transformation function:
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
vst_physeq <-ASV_physeq_core
otu_table(vst_physeq)<-vst_count_phy

my_pal<-c("#EE0000FF","#1B1919FF")
# generating and visualizing the PCoA with phyloseq

T1d_physeq<-subset_samples(physeq = ASV_physeq_core,disease=="T1D")
T1d_physeq<-subset_samples(physeq = T1d_physeq,patient!="20")
T1d_physeq<-subset_samples(physeq = T1d_physeq,patient!="49")
T1d_physeq<-subset_samples(physeq = T1d_physeq,patient!="50")
T1d_physeq
mother<-meta%>%
  filter(disease=="T1D")%>%
  filter(Host=="Mother")%>%filter(!patient%in%c(20,49,50))%>%
  mutate(HBA1c_cutoff=HbA1C_1trym_2>7.2)%>%
  select(patient,disease,HBA1c_cutoff)
mother2<-data.frame(inner_join(mother,meta)%>%distinct_all())
mother2$Delivery
nsamples(T1d_physeq)
dim(mother2)
rownames(mother2)<-mother2$SampleID
sample_data(T1d_physeq)<-sample_data(mother2)

mother2$Host

deseq_counts <- phyloseq_to_deseq2(T1d_physeq,design = ~HBA1c_cutoff)
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# now followed by the transformation function:
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
vst_physeq <-T1d_physeq
otu_table(vst_physeq)<-vst_count_phy

T1d_physeq_comp<-microbiome::transform(T1d_physeq,transform = "compositional")
vst_bray<-vegdist(x = t(otu_table(T1d_physeq_comp)),method = "bray")

infant<-mother2%>%filter(Host=="Infant")
infant_T1d_physeq<-subset_samples(T1d_physeq,Host=="Infant")
infant_T1d_physeq_comp<-microbiome::transform(infant_T1d_physeq,transform = "compositional")
infant_T1d_physeq_bray<-vegdist(x = t(otu_table(infant_T1d_physeq)),method = "bray")
infant_res<-pairwise.adonis(x = infant_T1d_physeq_bray,factors = infant$HBA1c_cutoff,perm =  9999)

mother<-mother2%>%filter(Host=="Mother")
mother_T1d_physeq<-subset_samples(T1d_physeq,Host=="Mother")
mother_T1d_physeq_comp<-microbiome::transform(mother_T1d_physeq,transform = "compositional")
mother_T1d_physeq_bray<-vegdist(x = t(otu_table(mother_T1d_physeq)),method = "bray")
mother_res<-pairwise.adonis(x = mother_T1d_physeq_bray,factors = mother$HBA1c_cutoff,perm =  9999)
mother_res

######################## Delivery #######################
infant<-

summary(res)
eigen_vals <- vst_bray_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_bray_pcoa, color="disease") + 
  geom_point(size=1) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values = my_pal)+
  theme_bw() + theme(legend.position="right")
#remove all columns with na
not_any_na <- function(x) all(!is.na(x))
sam<-meta%>%filter(SampleType=="Introitus")
sam<-as.character(sam$SampleID)
Day<-prune_samples(x = vst_physeq,sam)
Day<-core(x = Day,detection = 0/100,prevalence = 1/100)
Day<-prune_taxa(x = Day,taxa = taxa_sums(Day)>0)
Day_unifrac<-vegdist(x = t(otu_table(Day)),method = "bray")
Day_meta<-meta%>%
  filter(SampleID%in%sample_names(Day))%>%
  select(disease)
dim(Day_meta)
dim(as.matrix(Day_unifrac))

res<-adonis2(Day_unifrac~disease,data = Day_meta,perm = 9999,parallel = T)
res
return(res)
}

day_adonis<-function(X){
  sam<-meta%>%filter(SampleType==X)
  sam<-as.character(sam$SampleID)
  Day<-prune_samples(x = vst_physeq,sam)
  Day<-core(x = Day,detection = 0/100,prevalence = 1/100)
  Day<-prune_taxa(x = Day,taxa = taxa_sums(Day)>0)
  Day_unifrac<-vegdist(x = t(otu_table(Day)),method = "bray")
  Day_meta<-meta%>%
    filter(SampleID%in%sample_names(Day))%>%
    select(disease)
  res<-adonis2(Day_unifrac~disease,data = Day_meta,perm = 9999,parallel = T)
  return(res)
}
list<-meta%>%select(SampleType)%>%distinct_all%>%arrange(SampleType)
list
lapply(list$SampleType, function(X) day_adonis(X = X))
