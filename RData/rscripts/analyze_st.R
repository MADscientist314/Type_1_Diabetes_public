library(DESeq2)
library(ggplot2)
library(microbiome)
library(vegan)
library(EnhancedVolcano)
#setwd("~/Type1Diabetes/")
# #and lets make a new phyloseq object for each sampletype


infant<-subset_samples(ASV_physeq_core,Host=="Infant")
#infant<-tax_glom(infant,taxrank="Genus")
infant
filename<-ASV_physeq_core

  count_tab<-as.data.frame(otu_table(filename))
  tax_tab<-as.matrix(tax_table(filename))
  sample_info_tab<-data.frame(sample_data(filename))
  sample_info_tab
  
  dim(sample_info_tab)
  # first we need to make a DESeq2 object
  deseq_counts <- phyloseq_to_deseq2(filename,design = ~SampleType+disease)
  deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
  # now followed by the transformation function:
  deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
  # and here is pulling out our transformed table
  vst_trans_count_tab <- assay(deseq_counts_vst)
  # and calculating our Euclidean distance matrix
  euc_dist <- dist(t(vst_trans_count_tab))
  euc_clust <- hclust(euc_dist, method="ward.D2")
  # making our phyloseq object with t
  # making our phyloseq object with transformed table
  vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
  sample_info_tab_phy <- sample_data(sample_info_tab)
  vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)
  
  # generating and visualizing the PCoA with phyloseq
  vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
  eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
  
  plot_ordination(vst_physeq, vst_pcoa, color="disease") + 
    geom_point(size=3) + labs(col="disease") + scale_color_aaas()+theme_bw()+
    coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle(paste("neonatal PCoA"))  +
    facet_grid(SampleType~Antibiotics_2+Delivery)
  
  
  #######BETADISPER and ADONIS  
  #Test for homogeneity within groups
  
  #######BETADISPER and ADONIS  
  #Test for homogeneity within groups
  anova(betadisper(euc_dist, sample_info_tab$disease))
  #Perumtational test for significance
  adonis(data = sample_info_tab,formula = euc_dist~disease*SampleType) # 0.003
  # now converting our phyloseq object to a deseq object
  deseq <- phyloseq_to_deseq2(filename, ~ disease)
  # and running deseq standard analysis:
  deseq <- DESeq(deseq)
  # pulling out our results table, we specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about
  deseq_res <- results(deseq, alpha=0.01, contrast=c("disease", "Control", "T1D"))
  # we can get a glimpse at what this table currently holds with the summary command
  summary(deseq_res) 
  # this tells us out of ~1,800 ASVs, with adj-p < 0.01, there are X increased when comparing Control to T1D, and about X decreased
  # "decreased" in this case means at a lower count abundance in the Control than in the T1d, and "increased" means greater proportion in Control than in T1d
  # let's subset this table to only include these that pass our specified significance level
  sigtab_res_deseq<- deseq_res[which(deseq_res$padj < 0.01), ]
  
  # now we can see this table only contains those we consider significantly differentially abundant
  summary(sigtab_res_deseq) 
  
  # next let's stitch that together with these ASV's taxonomic annotations for a quick look at both together
  sigtab_deseq_with_tax <- cbind(as(sigtab_res_deseq, "data.frame"), as(tax_table(filename)[row.names(sigtab_res_deseq), ], "matrix"))
  
  # and now let's sort that table by the baseMean column
  sigtab_deseq_with_tax[order(sigtab_deseq_with_tax$baseMean, decreasing=T), ]
  # this puts sequence derived from an Actinobacteria at the top that was detected in ~10 log2fold greater abundance in the glassy basalts than in the more highly altered basalts
  write.table(sigtab_deseq_with_tax, file = paste(st,"DESEQ2_VST_results.tsv"), sep = "\t")
  
  
  res1 <- lfcShrink(dds = deseq,contrast = c('disease','Control','T1D'), res=deseq_res, type = 'normal')
  # ####BASIC#######
  EnhancedVolcano(res1,
                  lab = rownames(res1),
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  xlim = c(-8, 8),
                  title = paste(st,' Control versus Type 1 Diabetes'))
  ########Modified cutoffs for pvalue and log2FC
  EnhancedVolcano(res1,
                  lab = sigtab_deseq_with_tax$Genus,
                  x = 'log2FoldChange',
                  y = 'pvalue',
                  xlim = c(-8, 8),
                  title =  paste(st,' Control versus Type 1 Diabetes'),
                  pCutoff = 0.01,
                  FCcutoff = 1.5)
  
}

ASV_physeq_core
ASV_physeq_core_Introitus<-subset_samples(ASV_physeq_core, SampleType==" Introitus")
ASV_physeq_core_Vagina<-subset_samples(ASV_physeq_core, SampleType=="Vagina")
ASV_physeq_core_Cervix<-subset_samples(ASV_physeq_core, SampleType=="Cervix")
ASV_physeq_core_Stool<-subset_samples(ASV_physeq_core, SampleType=="Stool")
ASV_physeq_core_Anus<-subset_samples(ASV_physeq_core, SampleType=="Anus")
ASV_physeq_core_Ear<-subset_samples(ASV_physeq_core, SampleType=="Ear")

analyze_st(ASV_physeq_core_Introitus)
analyze_st(ASV_physeq_core_Vagina)
analyze_st(ASV_physeq_core_Cervix)
analyze_st(ASV_physeq_core_Anus)
analyze_st(ASV_physeq_core_Ear)
analyze_st(ASV_physeq_core_Stool)


getwd()