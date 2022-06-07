library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library(microbiome)
library(vegan)
library(EnhancedVolcano)
st<-as.character(meta(ASV_physeq_core_Introitus)$SampleType[1])
count_tab<-as.data.frame(otu_table(ASV_physeq_core_Introitus))
tax_tab<-as.matrix(tax_table(ASV_physeq_core_Introitus))
sample_info_tab<-as.data.frame(sample_data(ASV_physeq_core_Introitus))
# first we need to make a DESeq2 object
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~ disease) 
# we have to include the "colData" and "design" arguments because they are 
# required, as they are needed for further downstream processing by DESeq2, 
# but for our purposes of simply transforming the data right now, they don't 
# matter

#deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# NOTE: If you get this error here with your dataset: "Error in
# estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
# gene contains at least one zero, cannot compute log geometric means", that
# can be because the count table is sparse with many zeroes, which is common
# with marker-gene surveys. In that case you'd need to use a specific
# function first that is equipped to deal with that. You could run:
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# now followed by the transformation function:
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)
# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 
# but i like to change them to dendrograms for two reasons:
# 1) it's easier to color the dendrogram plot by groups
# 2) if wanted you can rotate clusters with the rotate() 
#    function of the dendextend package
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_info_tab$colors[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols

plot(euc_dend, ylab="VST Euc. dist.")
#dend_cols <- as.character(sample_info_tab$disease[order.dendrogram(euc_dend)])
#labels(euc_dend) <- dend_cols

plot(euc_dend, ylab="VST Euc. dist.")
# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="disease_SampleType") + 
  geom_point(size=3) + labs(col="disease") + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle(paste(st,"PCoA"))  


#######BETADISPER and ADONIS  
#Test for homogeneity within groups
anova(betadisper(euc_dist, sample_info_tab$disease))
#Perumtational test for significance
adonis(euc_dist~sample_info_tab$disease) # 0.003
# now converting our phyloseq object to a deseq object
deseq <- phyloseq_to_deseq2(ASV_physeq_core_Introitus, ~ disease)
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
sigtab_deseq_with_tax <- cbind(as(sigtab_res_deseq, "data.frame"), as(tax_table(ASV_physeq_core_Introitus)[row.names(sigtab_res_deseq), ], "matrix"))

# and now let's sort that table by the baseMean column
sigtab_deseq_with_tax[order(sigtab_deseq_with_tax$baseMean, decreasing=T), ]
# this puts sequence derived from an Actinobacteria at the top that was detected in ~10 log2fold greater abundance in the glassy basalts than in the more highly altered basalts
write.table(sigtab_deseq_with_tax, file = paste(st,"DESEQ2_VST_results.tsv"), sep = "\t")


res1 <- lfcShrink(dds = deseq,contrast = c('disease','Control','T1D'), res=deseq_res, type = 'normal')
# ####BASIC#######
EnhancedVolcano(res1,
                lab = paste(sigtab_deseq_with_tax$Genus,"_",rownames(sigtab_deseq_with_tax)),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-8, 8),
                title = paste(st,' Control versus Type 1 Diabetes'))
# ########Modified cutoffs for pvalue and log2FC
# EnhancedVolcano(res1,
#                 lab = paste(sigtab_deseq_with_tax$Genus,"_",rownames(sigtab_deseq_with_tax)),
#                 x = 'log2FoldChange',
#                 y = 'pvalue',
#                 xlim = c(-8, 8),
#                 title =  paste(st,' Control versus Type 1 Diabetes'),
#                 pCutoff = 0.01,
#                 FCcutoff = 1.5,
#                 pointSize = 3.0,
#                 labSize = 3.0)
