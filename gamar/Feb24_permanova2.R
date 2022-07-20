library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library(microbiome)
library(microbiomeutilities)


count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")
tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
sample_info_tab <- read.table("10_Feb_2020_metadata_MG.csv", header=T, row.names=1,
                              check.names=F, sep=",")

# and setting the color column to be of type "character", which helps later
sample_info_tab # to take a peek

count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
sample_info_tab_phy <- sample_data(sample_info_tab)
ASV_physeq<- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)  # and now we can call the plot_richness() function on our phyloseq object
#to remove 2,4,5,6

ASV_physeq<-subset_samples(ASV_physeq, patient_disease!="2_T1D")
ASV_physeq<-subset_samples(ASV_physeq, patient_disease!="4_T1D")
ASV_physeq<-subset_samples(ASV_physeq, patient_disease!="5_T1D")
ASV_physeq<-subset_samples(ASV_physeq, patient_disease!="6_T1D")

ASV_physeq



ASV_physeq <- subset_taxa(ASV_physeq, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
  
  
  

ASV_physeq_trim<-prune_samples(sample_sums(ASV_physeq)>0, ASV_physeq)




ASV_physeq_trim
ASV_physeq_trim<- format_to_besthit(ASV_physeq_trim)
ASV_physeq_rel <-  microbiome::transform(ASV_physeq_trim, "compositional")
otu <- abundances(ASV_physeq_rel)

meta <- meta(ASV_physeq_rel)




library(vegan)
#> Loading required package: permute
#> Loading required package: lattice
#> This is vegan 2.5-3
#> 
#> Attaching package: 'vegan'
#> The following object is masked from 'package:microbiome':
#> 
#>     diversity

library('DESeq2')
# first we need to make a DESeq2 object



permanova <- adonis(t(otu) ~ disease_SampleType,
                    data = meta, permutations = 999, method = "bray")# P-value
permanova


#ear
ear<-subset_samples(ASV_physeq_core,SampleType=="Ear")
deseq_counts <- DESeqDataSetFromMatrix(data.frame(otu_table(ear)), colData = sample_data(ear), design = ~disease) 
meta<-meta(ear)
#deseq_counts <- estimateSizeFactors(deseq_counts, type="poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))
anova(betadisper(euc_dist, meta$disease))
adonis(euc_dist~meta$disease)
ear_permanova <- adonis(euc_dist ~ disease_SampleType, data = meta, permutations = 999, method = "bray")# P-value
ear_permanova
euc_clust<-hclust(euc_dist, method = "ward.D2")
plot(euc_clust, ylab="VST Eud dend")


#cervix
#ear
cervix<-subset_samples(ASV_physeq_core,SampleType=="Cervix")
deseq_counts <- DESeqDataSetFromMatrix(data.frame(otu_table(cervix)), colData = sample_data(cervix), design = ~disease) 
meta<-meta(cervix)
#deseq_counts <- estimateSizeFactors(deseq_counts, type="poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))
anova(betadisper(euc_dist, meta$disease))
cervix_permanova <- adonis(euc_dist ~ disease_SampleType, data = meta, permutations = 999, method = "bray")# P-value
cervix_permanova
euc_clust<-hclust(euc_dist, method = "ward.D2")
plot(euc_clust, ylab="VST Eud dend")


#stool

#vagina

#Anus


#introitus = vestibule


#SHINY
save.image(file = "ASV_physeq2.RData")
saveRDS(file = "ASV_physeq.rds", ASV_physeq_trim)
getwd()
#install.packages("shiny")
library(shiny)
shiny::runGitHub("shiny-phyloseq","joey711")


#Feb18


ASV_physeq

write.csv((otu_table(ASV_physeq)), "otu.tsv")
write.csv((tax_table(ASV_physeq)), "tax.tsv")
write.csv((sample_data(ASV_physeq)), "sample_data.tsv")

