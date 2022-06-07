library(dplyr)
library(DESeq2)
library(microbiome)

ASV_physeq_core
#Lets agllomerate to different taxonomic levels, and transform compositionally
ASV_physeq_core<-tax_glom(physeq = ASV_physeq_core, taxrank = "Genus", NArm = T)
ASV_physeq_core

ASV_physeq_core<-microbiome::transform(ASV_physeq_core, 'compositional')
#Lets rename the taxa so the show up on the plot correctly
taxa_names(ASV_physeq_core)<-get_taxa_unique(ASV_physeq_core, taxonomic.rank = "Genus", errorIfNULL = T)

st<-as.character(meta(ASV_physeq_core)$SampleType[1])
dds <- phyloseq_to_deseq2(ASV_physeq_core, ~ disease)
geoMeans<-apply(counts(dds), 1, gm_mean)
dds<-estimateSizeFactors(dds, geoMeans=geoMeans)
dds <- DESeq(dds,test = "Wald",fitType = "local",parallel = T)
dds <- DESeq(dds, betaPrior=FALSE, parallel = T)
res1 <- results(dds,contrast = c('disease','Control','T1D'))
res1 <- lfcShrink(dds,contrast = c('disease','Control','T1D'), res=res1, type = 'normal')
# ####BASIC#######
# EnhancedVolcano(res1,
#                 lab = rownames(res1),


ASV_physeq_core_Introitus <- subset_samples(ASV_physeq_core_Introitus, disease != "NA")
ASV_physeq_core_Introitus <- prune_samples(sample_sums(kostic) > 0, ASV_physeq_core_Introitus)
ASV_physeq_core_Introitus
sample_data(ASV_physeq_core_Introitus)$disease
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


diagdds = phyloseq_to_deseq2(ASV_physeq_core_Introitus, ~ disease)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds,test = "Wald",sfType = 'poscounts',fitType = "local",parallel = T)
