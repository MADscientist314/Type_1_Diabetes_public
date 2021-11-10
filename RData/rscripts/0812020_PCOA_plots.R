
library(microbiome)
sam<-meta(ASV_physeq_core)
dim(sam)
sam$disease_SampleType<-paste(meta(ASV_physeq_core)$SampleType,"_",meta(ASV_physeq_core)$disease)
sample_data(ASV_physeq_core)<-sample_data(sam)
# plotPCA(deseq_counts_vst, intgroup="disease_SampleType" )+
#   # facet_wrap(facets = ~deseq_counts_vst$SampleType)+
#   stat_ellipse(aes(group = disease_SampleType), linetype = 2)+
#   scale_color_brewer(palette = "Paired")+
#   theme_classic()

ASV_physeq_core_clr<-transform(ASV_physeq_core, "clr")
#PCA via phyloseq
ord_clr <- phyloseq::ordinate(ASV_physeq_core, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))     
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(ASV_physeq_core, ord_clr, type="samples", color="disease_SampleType", shape="disease") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  # facet_wrap(facets = ~SampleType)+
  stat_ellipse(aes(group = disease_SampleType), linetype = 2)+
  scale_color_brewer(palette = "Paired")+
  theme_classic()

######################################################


# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

#adding a custom colum that has disease and sampletype
sam<-meta(vst_physeq)
sam$disease_SampleType<-paste(meta(vst_physeq)$SampleType,"_",meta(vst_physeq)$disease)
sample_data(vst_physeq)<-sample_data(sam)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(vst_physeq, vst_pcoa, color="disease_SampleType") +
  facet_grid(facets = disease_SampleType~Host,switch = "both", scales = "free")+
  # geom_point(size=1) + labs(col="SampleType") + 
  stat_ellipse(aes(group = disease_SampleType), linetype = 2)+
  # geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  #coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_brewer(palette = "Paired")+
  theme_bw()

#####################################################################################

a<-plotPCA(deseq_counts_vst, intgroup="energ_bialko_proc_2")+
  scale_color_viridis_c(direction = -1)+
  aes(shape=deseq_counts_vst$disease)+
  facet_wrap(facets = ~deseq_counts_vst$SampleType)+
  theme_bw()
a + geom_smooth(mapping = aes(fill=deseq_counts_vst$disease),
                alpha=0.1,
                method = "lm",
                inherit.aes = T,
                show.legend = T,
                na.rm = T)
library(Rtsne)

method <- "tsne"
trans <- "hellinger"
distance <- "euclidean"

# Distance matrix for samples
ps <- microbiome::transform(ASV_physeq_core, trans)

# Calculate sample similarities
dm <- vegdist(otu_table(ps), distance)

# Run TSNE
tsne_out <- Rtsne(dm, dims = 2) 
proj <- tsne_out$Y
rownames(proj) <- rownames(otu_table(ASV_physeq_core))
tsne_out()
library(ggplot2)
p <- plot_landscape(proj, legend = T, size = 1,col = "energ_bialko_proc_2") 
print(p)
                  
library(tsnemicrobiota)
library(ggforce)
tsne_res <- tsne_phyloseq(ASV_physeq_core, distance='bray',
                          perplexity = 50, verbose=0, rng_seed = 3901)
plot_tsne_phyloseq(ASV_physeq_core, tsne_res,color="disease",shape="disease",
                   title='t-SNE (Bray-Curtis)')+ 
  facet_grid(facets = meta(ASV_physeq_core)$tydzien_porodu_2~meta(ASV_physeq_core)$SampleType)+
  geom_mark_ellipse(aes(color = disease),expand = 0.01) +
  geom_jitter(size=2) + labs(col="SampleType") +
  geom_smooth(method = "lm",se = F, na.rm = T, linetype = 2)+
  scale_color_brewer(palette = "Set1")+
  theme_bw()
###########################################################################

ggplot(data =  meta(ASV_physeq_core))+
  geom_histogram(mapping = aes(x = tydzien_porodu_2, fill=disease),
                 bins = 50,binwidth = 0.5,
                 na.rm = T,
                 position = "dodge") +
  xlab(label = "Week of Delivery") +
  ylab(label = "Number of samples")+
  scale_fill_brewer(palette = "Set1")+
  theme_bw()


#######################################################################

plot_ordination(vst_physeq, vst_pcoa, color="disease") +
  facet_grid(facets = tydzien_porodu_2~SampleType, scales = "free")+
  geom_mark_ellipse(aes(fill = disease),alpha=0.07,expand = 0.01) +
  geom_jitter(size=1) + labs(col="SampleType") + 
  #geom_polygon(mapping = aes(group=disease, fill=disease),alpha=0.4,na.rm = T)+
  #stat_ellipse(aes(group = disease), linetype = 1)+
  # geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  #coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_brewer(palette = "Set1")+
  theme_bw()
