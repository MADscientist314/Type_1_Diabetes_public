library(microbiome)
ASV_physeq_core<-readRDS("ASV_physeq_core.RDS")
maternal<-subset_samples(physeq = ASV_physeq_core,Host=="Mother")
maternal
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
ord_clr <- phyloseq::ordinate(ASV_physeq_core_clr, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")
sapply(ord_clr$CA$eig[1:5], function(x) x / sum(ord_clr$CA$eig))     
clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(ASV_physeq_core, ord_clr, type="samples", color="disease", shape="SampleType") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  # facet_wrap(facets = ~SampleType)+
  #  stat_ellipse(aes(group = disease_SampleType), linetype = 2)+
  scale_color_manual(values = my_pal)+
  theme_classic()
######################################################
library(DESeq2)

count_tab<-as.data.frame(otu_table(ASV_physeq_core))
tax_tab<-as.matrix(tax_table(ASV_physeq_core))
sample_info_tab<-as.data.frame(sample_data(ASV_physeq_core))
# first we need to make a DESeq2 object
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~ 1) 
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
ds<-phyloseq_to_deseq2(ASV_physeq_core,design = ~1)
deseq_counts <- estimateSizeFactors(ds, type = "poscounts")
# now followed by the transformation function:
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts,blind = T)
library(parallel)
deseq_counts_vst
vst_trans_count_tab<-otu_table(object = deseq_counts_vst,taxa_are_rows = T)
deseq_counts <- phyloseq_to_deseq2(ASV_physeq_core_genus, design = ~ 1)
vst_trans_count_tab <- assay(deseq_counts_vst)
# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(ASV_physeq_core)
vst_physeq <- phyloseq(vst_count_phy, tax_table(ASV_physeq_core),sample_data(ASV_physeq_core))
#adding a custom colum that has disease and sampletype
sam<-meta(vst_physeq)
sam$disease_SampleType<-paste(meta(vst_physeq)$SampleType,"_",meta(vst_physeq)$disease)
sample_data(vst_physeq)<-sample_data(sam)
library(ggpubr)
library(ggsci)
library(scales)
show_col(pal_aaas(palette = "default",alpha = 1)(10))
pal_aaas(palette = "default",alpha = 1)(10)
my_pal<-c("#BB0021FF","#3B4992FF","#1B1919FF","#EE0000FF","#008B45FF","#631879FF","#BB0021FF","#5F559BFF","#A20056FF","#808180FF","#1B1919FF")
show_col(my_pal)
c("unifrac","wunifrac","dpcoa","jsd")
distances<-c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord")
methods<-c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
vst_physeq
# generating and visualizing the PCoA with phyloseq
vst_pcoa<- ordinate(vst_physeq, 
                    method="MDS", 
                    distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
ord<-plot_ordination(vst_physeq, vst_pcoa, color="disease",shape="SampleType") +
  scale_color_manual(values = my_pal)+
  facet_grid(facets = Host~.,margins = T)+
  stat_ellipse(aes(group = SampleType), linetype = 3) +
  #geom_point(size=1) + 
  #labs(col="SampleType") + 
  #geom_text(aes(label=rownames(sam), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  theme_pubr()
ord
#####################################################################################
ds<-phyloseq_to_deseq2(ASV_physeq_core,design = ~1)
deseq_counts <- estimateSizeFactors(ds, type = "poscounts")
# now followed by the transformation function:

deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
deseq_counts_vst

my_pal<-c("#1B1919FF","#BB0021FF",,"#E0000FF","#008B45FF","#631879FF","#BB0021FF","#5F559BFF","#A20056FF","#808180FF","#1B1919FF","#3B4992FF")
show_col(my_pal)
ntop = 500, returnData = FALSE)
p<-plotPCA(object = deseq_counts_vst,  
        intgroup="disease", ntop = 500, returnData = T)
p
a<-plotPCA(object = deseq_counts_vst, 
           intgroup="disease", ntop = 500, returnData = FALSE)+
  aes(shape=deseq_counts_vst$SampleType)+
  facet_grid(facets = deseq_counts_vst$Host~.,switch = "x")+
  geom_smooth(mapping = aes(color=deseq_counts_vst$disease,fill=deseq_counts_vst$disease),
              alpha=0.1,
              method = "lm",
              inherit.aes = T,
              show.legend = T,
              na.rm = T)+
  scale_fill_manual(values = my_pal)+
  scale_color_manual(values = my_pal)+
  theme_pubr(margin = T,border = T,legend = "right",base_size = 12,base_family = "mono")
a


ggscatter(data = p,
          x = "PC1",
          y = "PC2",
          color = "disease",
          shape = "SampleType",
          rug = T,palette = my_pal,
          fill = "disease",
          title = "VST transformed")


# using rlog transformed data:
dds <- makeExampleDESeqDataSet(betaSD=1)
dds
rld <- rlog(deseq_counts)
a<-plotPCA(object = rld, 
           intgroup="disease")+
  aes(shape=deseq_counts_vst$SampleType)+
  facet_grid(facets = deseq_counts_vst$Host~.)+
  geom_smooth(mapping = aes(fill=deseq_counts_vst$SampleType),
              alpha=0.1,
              method = "lm",
              inherit.aes = T,
              show.legend = T,
              na.rm = T)+
  scale_fill_manual(values = my_pal)+
  scale_color_manual(values = my_pal)+
  theme_pubclean()
# also possible to perform custom transformation:
ds<-phyloseq_to_deseq2(ASV_physeq_core,design = ~1)
deseq_counts <- estimateSizeFactors(ds, type = "poscounts")
dds <- estimateSizeFactors(deseq_counts)
# shifted log of normalized counts
se <- SummarizedExperiment(log(counts(deseq_counts, normalized=TRUE) + 1),colData=colData(deseq_counts))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method.
plotPCA( DESeqTransform(SummarizedExperiment = se),  intgroup="disease")+
  aes(shape=deseq_counts_vst$SampleType)+
  facet_grid(facets = deseq_counts_vst$Host~.)+
  # geom_smooth(mapping = aes(color=deseq_counts_vst$disease),
  #             alpha=0.1,
  #             method = "lm",
  #             inherit.aes = T,
  #             show.legend = T,
  #             na.rm = T)+
           scale_fill_manual(values = my_pal)+
           scale_color_manual(values = my_pal)+
           theme_pubclean()
library(Rtsne)
library(vegan)
method <- "tsne"
trans <- "none"
distance <- "bray"
# Distance matrix for samples
ps <- microbiome::transform(ASV_physeq_core, trans)
# Calculate sample similarities
dm <- vegdist(t(otu_table(vst_physeq)), distance)

# Run TSNE
nrow(matrix(dm)) - 1
tsne_out2 <- Rtsne(dm, 
                  dims = 2,
                  is_distance = T,
                  num_threads=48,
                  perplexity =  150 , 
                  theta = 0.5, 
                  check_duplicates = TRUE,
                  pca_center=T,
                  pca_scale=T,
                  pca = TRUE, 
                  partial_pca = FALSE, 
                  max_iter = 1000,
                  verbose = getOption("verbose", TRUE)) 

sam<-meta(dm)
sam$X<-tsne_out$Y[,1]
sam$Y<-tsne_out$Y[,2]

ggscatter(data = sam,
          x = "X",
          y = "Y",
          color = "disease",
          shape="SampleType",
          palette = my_pal,
          rug = T,
          add = "reg.line",
          #star.plot = T,
          add.params = list(size=0.8,color="disease",linetype=3),
          xlab = "tsne X",
          ylab = "tsne y",
          #ellipse = T,
          #ellipse.type = "euclid",
          title = "Beta diversity Bray Curtis tsne plot",
          ggtheme =   theme_pubr(margin = T,border = T,legend = "right",base_size = 12,base_family = "mono"))+
  #facet_wrap(~Host+SampleType,nrow = 3)
  facet_wrap(~Host,nrow = 3)




  
                      
  
  #labs(col="disease") +
  #geom_smooth(mapping = aes(group=disease,color=disease),inherit.aes = T,position = "identity",method = "lm",se = F, na.rm = T, linetype = 3)+
#  scale_fill_manual(values = my_pal)+
 # scale_color_manual(values = my_pal)+
  theme_pubclean(base_size = 12,base_family = "mono")

  theme_pubr(margin = T,border = F,legend = "right",base_size = 12,base_family = "mono")
###########################################################################

library(ggplot2)
p <- plot_landscape(proj, legend = T, size = 1,col = "disease",transformation = ) 
print(p)
library("devtools")
install_github("opisthokonta/tsnemicrobiota")
library("tsnemicrobiota")
library(ggforce)

tsne_res <- tsne_phyloseq(ASV_physeq_core, distance = "bray",
                          #precomputed_distance =dm,  
                          perplexity = 500, 
                          verbose=2, 
                          control = list(max_iter=10000, momentum_iter=c(250, 10000)),
                          rng_seed = 314)
plot_tsne_phyloseq(physeq = ps,
                   tsne_obj =  tsne_res,
                   color="disease",
                   shape="SampleType",
                   title='t-SNE (Bray-Curtis)')+ 
  facet_grid(facets = meta(ASV_physeq_core)$Host~.)+
  geom_mark_ellipse(aes(fill = disease),alpha=0.07,expand = 0.01) +
  #geom_jitter(size=2) + 
  #labs(col="disease") +
  #geom_smooth(method = "lm",se = F, na.rm = T, linetype = 3)+
  scale_fill_manual(values = my_pal)+
  scale_color_manual(values = my_pal)+
  theme_pubr(margin = T,border = T,legend = "right",base_size = 12,base_family = "mono")
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
