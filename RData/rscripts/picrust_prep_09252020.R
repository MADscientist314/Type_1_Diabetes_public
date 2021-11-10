library(phyloseq)
library(microbiome)
library(dplyr)
library(ggpubr)
library(ggsci)
library(DESeq2)
sam<-read.table(file = "sample_core_Sept16_2020.txt", header = T, sep = "\t", row.names = 1)
sam<-sample_data(sam)
getwd()
ASV_physeq_core
sample_data(ASV_physeq_core)<-sam
ASV_physeq_core<-subset_taxa(physeq = ASV_physeq_core, Kingdom!="Archaea")
#ASV_physeq_core_rel<-transform(ASV_physeq_core, "compositional")
ASV_physeq_core_genus<-tax_glom(ASV_physeq_core,taxrank = "Genus")
ASV_physeq_core_genus

#<-get_taxa_unique(ASV_physeq_core_genus, taxonomic.rank = "Genus")
#get_taxa_unique(ASV_physeq_core_genus, taxonomic.rank = "Genus")
taxa_sums(ASV_physeq_core)
sample_sums(ASV_physeq_core)
count_tab<-data.frame(otu_table(ASV_physeq_core_genus))
sample_info_tab<-data.frame(sample_data(ASV_physeq_core_genus))
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~disease) 

deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 
# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
tax_phy<-tax_table(ASV_physeq_core_genus)
vst_physeq <- phyloseq(vst_count_phy, tax_phy,sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="disease") + 
  geom_point(size=1) + labs(col="type") + 
#  geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
#  scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) + 
  theme(legend.position="right")

otu<-data.frame(otu_table(vst_physeq))
tax<-data.frame(tax_table(vst_physeq))
sam<-data.frame(sample_data(vst_physeq))

write.table(x = otu,file = "vst_counts.tsv",quote = F, sep = "\t")
write.table(x = tax,file = "vst_tax.tsv",quote = F, sep = "\t")
write.table(x = sam,file = "vst_metadata.tsv",quote = F, sep = "\t")
library(reltools)
get_taxa(physeq = vst_physeq,)
save_fasta(ps = ASV_physeq_core_genus, file = "vst_physeq.fasta", rank = "Genus")
write.fasta