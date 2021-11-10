# Load libraries
library(phyloseq)
library(DESeq2)
library(microbiome)
library(ggplot2)
library(dplyr)
library(vegan)
getwd()

## REMOVE CHLOROPLAST AND MITOCHONDRIAL OTU'S using phyloseq if not done using QIIME [Current data: DONE IN QIIME WITH filter_taxa_from_otu_table.py]
ASV_physeq_core <- subset_taxa(ASV_physeq_core, Class != "c__Chloroplast")
ASV_physeq_core <- subset_taxa(ASV_physeq_core, Family != "f__mitochondria")
sam<-read.csv(file = "vst_adonis_08102020_BACKUP/sample_coreAug7.txt", header = T, row.names = 1, sep = "\t")
sam<-sample_data(sam)
sample_data(ASV_physeq_core)<-sam
sam<-as.tibble(sam)
#############################################################################
#ASV_physeq_core_genus<-tax_glom(physeq = ASV_physeq_core,taxrank = "Genus")
deseq_counts <- phyloseq_to_deseq2(ASV_physeq_core, design = ~ disease)
deseq_counts
#deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~ disease) 
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
plot(euc_dend, ylab="VST Euc. dist.")

# generating and visualizing the PCoA with phyloseq
vst_physeq<-phyloseq(otu_table(vst_trans_count_tab, taxa_are_rows = T), sample_data(ASV_physeq_core))
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="disease", shape="disease") + 
  geom_point(size=2) + labs(col="disease") + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + 
  ggtitle(paste("SampleType","PCoA"))+
  facet_wrap(~SampleType)+
  theme_minimal()  


# Probiotics intervention example data 
# Pick relative abundances (compositional) and sample metadata 
# pseq.rel <- microbiome::transform(ASV_physeq_core, "compositional")
# otu <- abundances(pseq.rel)
# meta <- meta(pseq.rel)
# library(vegan)
# permanova <- adonis(t(otu) ~disease *SampleType* tydzien_porodu_2,strata =meta$SampleType,
#                     data = meta, permutations=999)
# 
# 
# permanova
# P-value
#print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])

#######BETADISPER and ADONIS  
#Test for homogeneity within groups
# sam<-select_(.data = sam, c(patient, disease, Delivery))
q<-distinct(select(data.frame(sam),c(patient, disease)))
q1 <- q %>% filter(!is.na(q$disease))
q1
##
beta<-betadisper(euc_dist,group =meta(vst_physeq)$disease)
permutest(beta)
TukeyHSD(beta)
capture.output(permutest(beta),file = "euc_vst_disease_within_group_anova_betadisper.tsv")
#+sample_info_tab$energ_bialko_proc_2
#Perumtational test for significance

capture.output(adonis(euc_dist~disease*SampleType*tydzien_porodu_2,strata =meta(vst_physeq)$SampleType,data = meta(vst_physeq), permutations = 9999),file = "euc_vst_disease_permanova.tsv")

#############I DONT KNOW ABOUT THIS, BUT I THINK I NEED TO NORMALIZE FOR DUPLICATE SAMPLING##########
#############THIS NEXT PART IS MY ATTEMPT TO DO SO FOR THE EUC DIST VST COUNT TAB ADONIS#############

sam<-data.frame(sample_data(ASV_physeq_core))
sam<-distinct(select(data.frame(sam),c(patient, disease,SampleType,tydzien_porodu_2)))
euc_dist_2 <- dist(t(vst_trans_count_tab[ , colnames(vst_trans_count_tab) %in% rownames(sam)]))

# and now making a sample info table with just the distinct samples
sample_info_tab_2 <- sample_info_tab[row.names(sample_info_tab) %in% rownames(sam), ]
sample_info_tab_2
# running betadisper on just these based on level of alteration as shown in the images above:
anova(betadisper(euc_dist_2, sample_info_tab_2$disease)) # 0.7

b<-adonis(formula = euc_dist_2~disease*SampleType*tydzien_porodu_2,
       data =  sample_info_tab_2,method = "bray",
       strata = sample_info_tab_2$SampleType, permutations = 9999) # 0.003
adonis(formula = euc_dist_2~disease*SampleType*tydzien_porodu_2,
       data =  sample_info_tab_2,method = "bray", permutations = 9999)

capture.output(b, file="euc_dist_2_adonis_disease_weekdeliv.txt")


###################################energ_bialko_proc_2########################################################################
tmp<-subset_samples(ASV_physeq_core, energ_bialko_proc_2!="NA")
tmp

deseq_counts <- phyloseq_to_deseq2(tmp, design = ~ disease)
deseq_counts
#deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~ disease) 
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
plot(euc_dend, ylab="VST Euc. dist.")

# generating and visualizing the PCoA with phyloseq
vst_physeq<-phyloseq(otu_table(vst_trans_count_tab, taxa_are_rows = T), sample_data(tmp))
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="energ_bialko_proc_2", shape="disease") + 
  geom_point(size=3) + labs(col="disease") + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + 
  ggtitle(paste("SampleType","PCoA"))+
  facet_wrap(~SampleType)+
  theme_minimal()  

#######BETADISPER and ADONIS  
#Test for homogeneity within groups
sample_info_tab<-data.frame(sample_data(tmp))
beta<-betadisper(euc_dist,group =  sample_info_tab$energ_bialko_proc_2)
permutest(beta)
TukeyHSD(beta)
capture.output(permutest(beta),file = "euc_vst_energ_bialko_proc_2_within_group_permutest_betadisper.tsv")
#+sample_info_tab$energ_bialko_proc_2
#Perumtational test for significance
library(vegan)
adonis(euc_dist~disease*energ_bialko_proc_2,strata =sample_info_tab$SampleType,data = sample_info_tab, permutations = 9999)

#############I DONT KNOW ABOUT THIS, BUT I THINK I NEED TO NORMALIZE FOR DUPLICATE SAMPLING##########
#############THIS NEXT PART IS MY ATTEMPT TO DO SO FOR THE EUC DIST VST COUNT TAB ADONIS#############

sam<-data.frame(sample_data(ASV_physeq_core))
sam<-distinct(select(data.frame(sam),c(patient, disease,SampleType,energ_bialko_proc_2)))
euc_dist_2 <- dist(t(vst_trans_count_tab[ , colnames(vst_trans_count_tab) %in% rownames(sam)]))

# and now making a sample info table with just the basalts
sample_info_tab_2 <- sample_info_tab[row.names(sample_info_tab) %in% rownames(sam), ]
sample_info_tab_2
# running betadisper on just these based on level of alteration as shown in the images above:
beta<-betadisper(euc_dist_2, sample_info_tab_2$energ_bialko_proc_2) # 0.7
permutest(beta)
capture.output(permutest(beta),file = "euc_vst_energ_bialko_proc_2_within_group_permutest_betadisper.tsv")


b<-adonis(formula = euc_dist_2~disease*SampleType*energ_bialko_proc_2,
          data =  sample_info_tab_2,method = "bray",
          strata = sample_info_tab_2$SampleType, permutations = 9999) # 0.003
b
capture.output(b, file="euc_dist_2_adonis_disease_energ_bialko_proc_2.txt")
