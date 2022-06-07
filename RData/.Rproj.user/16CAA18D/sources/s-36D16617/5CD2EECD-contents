library(knitr)
library(kableExtra)
library(gridExtra)
library(scales)
library(RColorBrewer)
library(tidyverse) #dplyr package included
library(stringr)
library(lme4)
library(nlme)
library(effects)
library(DESeq2)
library(qiime2R)
library(phyloseq)
library(DivNet)
library(vegan)
library(sjPlot)
library(DESeq2)
library(ggplot2)
library(microbiome)
library(vegan)
library(EnhancedVolcano)

############## CREATE DATAFRAMES, A PHYLOSEQ OBJECT, PRE-PROCESS DATA ######################## 
## REMOVE CHLOROPLAST AND MITOCHONDRIAL OTU'S using phyloseq if not done using QIIME [Current data: DONE IN QIIME WITH filter_taxa_from_otu_table.py]
ASV_physeq_core <- subset_taxa(ASV_physeq_core, Class != "c__Chloroplast")
ASV_physeq_core <- subset_taxa(ASV_physeq_core, Family != "f__mitochondria")
sam<-read.csv(file = "sample_coreJuly15.txt", header = T, row.names = 1, sep = "\t")
sam<-sample_data(sam)
sample_data(ASV_physeq_core)<-sam
sam<-as.data.frame(sam)
###########################DESeq2####################################################
# first we need to make a DESeq2 object
deseq_counts <- phyloseq_to_deseq2(ASV_physeq_core, design = ~ disease)

# NORMALIZE THE COUNTS USING VARIANCE STABILIZING TRANSFORMATION
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)

#MAKE A EUCLIDEAN DISTANCE MATRIX OF THE COUNT TAB
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
euc_dend <- as.dendrogram(euc_clust, hang=0.1)

#MAKE A PHYLOSEQ OBJECT WITH THE VST NORMALIZED DATA
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
vst_physeq <- phyloseq(vst_count_phy, sam)
length(sample_variables(ASV_physeq_core))
sample_info_tab$BMI_before_2
library(tidyverse)
tmp<-as.tibble(sam)
plot(tmp$Glu1)
hist(sample_info_tab$BMI_before_2)
# Use simulated abundance matrix
set.seed(711)
testOTU <- sample_data(ASV_physeq_core)
f1  <- filterfun_sample(topk(2))
wh1 <- genefilter_sample(testOTU, f1, A=2)
wh2 <- c(TRUE, TRUE, TRUE, FALSE, FALSE)
prune_taxa(wh1, testOTU)
prune_taxa(wh2, testOTU)

cont[1]
cont<-c('tydzien_porodu_2',	'newborn_weight',	'Glu1',	'Glu2',	'Glu3',	'Glu4',	'Glu5',	'Glu_mean',	'energ_tluszcz_proc_2',	'energ_bialko_proc_2',	'Bialko_ogolem_g_2',	'BMI_before_2',	'BMI_after_2',	'HbA1C_1trym_2',	'HbA1C_2trym_2',	'HBA1C_POROD_2',	'insulin_age',	'insulin_how_long')
#ADONIS
library("vegan")

# df= sample_data(ASV_physeq_core)
# OTU = t(otu_table(ASV_physeq_core))
# permutation<-adonis(OTU~disease+tydzien_porodu_2,data = as(df, "data.frame"), permutations = 999)
# dist<-vegdist(OTU)
# permutation
# anova(betadiver(x = dist, group = c(df$disease)))
# permustats(permanova)
# anova(OTU~age_binned+SampleType+Patient,data=df)
# permutation<-adonis(OTU~SampleType+Patient+Day*SampleType+Day*Patient,data=df, permutations = 999)

microbes.adonis <- adonis2(formula = euc_dist ~  tydzien_porodu_2,
                           data =  as(sample_data(p), "data.frame"),
                           permutations = 999)
summary(microbes.adonis)

microbes.adonis

mutate()
# microbedat.phyloseq.time <- phyloseq_to_deseq2(subset_samples(p,euc_dist ~ disease + val)


# # generating and visualizing the PCoA with phyloseq
# vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
# eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
# 
# plot_ordination(vst_physeq, vst_pcoa, color="disease_SampleType") + 
#   geom_point(size=3) + labs(col="disease") + 
#   coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle(paste("PCoA"))  
# #################################################################################
# d = Bray(physeq = ASV_physeq_core, parallel = T)#tax_glom(microbedat, taxrank = "Genus"))
# d.mds <- metaMDS(comm = euc_dist, zerodist=ignore)
# 
# scrs <- cbind(sample_info_tab, data.frame(MDS1 = d.mds$points[,1], MDS2 = d.mds$points[,2])) 
# sam
# cent <- cbind(sample_info_tab, data.frame(MDS1 = d.mds$points[,1], MDS2 = d.mds$points[,2])) %>% 
#   aggregate(cbind(MDS1, MDS2) ~ disease+
#               SampleType +
#               Glu1 + 
#               Glu2 +
#               Glu3 +
#               Glu4 +
#               BMI_before_2 +
#               BMI_after_2 +
#               HbA1C_1trym_2 + 
#               HbA1C_2trym_2 + 
#               HBA1C_POROD_2 + 
#               insulin_age + 
#               insulin_how_long+
#               Antibiotics_2 , 
#             data = ., FUN = mean) 
# sample_data(ASV_physeq_core)
# segs <- merge(scrs, setNames(cent,c("disease","SampleType","BMI_before_2","BMI_after_2",'oNMDS1','oNMDS2')),
#               by = c("disease","SampleType","BMI_before_2","BMI_after_2"), sort = FALSE)
# 
# ggplot(scrs, mapping = aes(x = MDS1, y = MDS2, shape =disease, color = BMI_before_2)) +
#   scale_colour_viridis_b() +
#   geom_segment(data = segs,
#                mapping = aes(xend = oNMDS1, yend = oNMDS2)) + # spiders
#   geom_point(data = cent, size = 3) +                         # centroids
#   geom_point() +                                              # sample scores
#   coord_fixed()                                               # same axis scaling
###########################################################################################

