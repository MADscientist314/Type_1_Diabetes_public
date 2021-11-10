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
library(mosaic)

############## CREATE DATAFRAMES, A PHYLOSEQ OBJECT, PRE-PROCESS DATA ######################## 
## REMOVE CHLOROPLAST AND MITOCHONDRIAL OTU'S using phyloseq if not done using QIIME [Current data: DONE IN QIIME WITH filter_taxa_from_otu_table.py]
# and running deseq standard analysis:
ASV_physeq_core_genus<-tax_glom(physeq = ASV_physeq_core,taxrank = "Genus")
deseq_counts <- phyloseq_to_deseq2(ASV_physeq_core_genus, design = ~ disease)
deseq <- DESeq(deseq_counts)
# pulling out our results table, we specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res <- results(deseq, alpha=0.01, contrast=c("disease", "Control", "T1D"))
summary(deseq_res)
# let's subset this table to only include these that pass our specified significance level
sigtab_res_deseq <- deseq_res[which(deseq_res$padj < 0.01), ]
# now we can see this table only contains those we consider significantly differentially abundant
summary(sigtab_res_deseq)
# next let's stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_deseq_with_tax <- cbind(as(sigtab_res_deseq, "data.frame"), as(tax_table(ASV_physeq_core)[row.names(sigtab_res_deseq), ], "matrix"))
# and now let's sort that table by the baseMean column
sigtab_deseq_with_tax[order(sigtab_deseq_with_tax$baseMean, decreasing=T), ]
# pulling out our results table, we specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res <- results(deseq, alpha=0.01, contrast=c("disease", "Control", "T1D"))
summary(deseq_res)
# now we can see this table only contains those we consider significantly differentially abundant
capture.output(summary(sigtab_res_deseq),file = "vst_counts_deseq2_disease_results_padj_genus.tsv") 
# and now let's sort that table by the baseMean column
sigtab_deseq_with_tax[order(sigtab_deseq_with_tax$baseMean, decreasing=T), ]
# this puts sequence derived from an Actinobacteria at the top that was detected in ~10 log2fold greater abundance in the glassy basalts than in the more highly altered basalts
write.table(sigtab_deseq_with_tax, file = paste("DESEQ2_disease_VST_results.tsv"), sep = "\t")

res1 <- lfcShrink(dds = deseq,contrast = c('disease','Control','T1D'), res=deseq_res, type = 'normal')
# ####BASIC#######
sigtab_deseq_with_tax <- cbind(as(res1, "data.frame"), as(tax_table(ASV_physeq_core)[row.names(sigtab_res_deseq), ], "matrix"))
sigtab_deseq_with_tax

library(EnhancedVolcano)
EnhancedVolcano(sigtab_deseq_with_tax,
                lab = paste(sigtab_deseq_with_tax$Genus,"_",rownames(sigtab_deseq_with_tax)),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-8, 8),
                title = paste('DESEQ2 VST Control versus Type 1 Diabetes'))
# ########Modified cutoffs for pvalue and log2FC
EnhancedVolcano(sigtab_deseq_with_tax,
                lab = paste(sigtab_deseq_with_tax$Genus,"_",rownames(sigtab_deseq_with_tax)),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 5),
                title =  paste('DESEQ2 VST Control versus Type 1 Diabetes'),
                pCutoff = 0.01,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 3.0)

library(microbiome)
length(rownames(sigtab_deseq_with_tax))
tmp<-transform(Cervix, "compositional")
tmp<-prune_taxa(taxa = rownames(sigtab_deseq_with_tax), x = tmp)
tmp
taxa_names(tmp)<-sigtab_deseq_with_tax$Genus
length(taxa_names(tmp))
# for (z in (taxa_names(tmp))){
# a<-boxplot_abundance(d = tmp, x = "SampleType",y=z, violin = F, show.points = T, na.rm = T)+
#   scale_color_discrete("disease")
# print(a)
# }
tmp2<-as.data.frame(otu_table(tmp))
rownames(t(tmp2))==rownames(sample_data(tmp))
tmp3<-merge.data.frame(t(tmp2),sample_data(tmp))
tmp2<-psmelt(tmp)
tmp2$Abundance<-rownames(tmp2)
tmp2
tmp2$disease

ggplot(tmp2,mapping = aes(x = disease, y = Abundance, color=disease))+
  geom_boxplot(na.rm = T, outlier.color = "grey", outlier.fill = "grey",position = "dodge")
  #ylim(0,20)+
#  facet_wrap(facets = ~ OTU)
  #theme_minimal()+
  #theme(panel.margin=unit(.05, "lines"),
  #      panel.border = element_rect(color = "black", fill = NA, size = 0.1), 
  #      strip.background = element_rect(color = "black", size = 0.1))


