#####Import the necessary libraries)
library('EnhancedVolcano')
library('magrittr')
library('DESeq2')
library('microbiome')
##################################
# Make a function to run DESEQ2 differential abundance analyses
# to look for statitistically significance taxa associated with
# each one of the sample types 
##################################################

analyze_volcano<-function(filename){

  #Lets agllomerate to different taxonomic levels, and transform compositionally
  filename<-tax_glom(physeq = filename, taxrank = "Genus", NArm = T)
  filename<-microbiome::transform(filename, 'compositional')
  #Lets rename the taxa so the show up on the plot correctly
  taxa_names(filename)<-get_taxa_unique(filename, taxonomic.rank = "Genus", errorIfNULL = T)
  
  st<-as.character(meta(filename)$SampleType[1])
  dds <- phyloseq_to_deseq2(filename, ~ disease)
  geoMeans<-apply(counts(dds), 1, gm_mean)
  dds<-estimateSizeFactors(dds, geoMeans=geoMeans)
  dds <- DESeq(dds,test = "Wald",fitType = "local",parallel = T)
  # #dds <- DESeq(dds, betaPrior=FALSE, parallel = T)
  # res1 <- results(dds,contrast = c('disease','Control','T1D'))
  # res1 <- lfcShrink(dds,contrast = c('disease','Control','T1D'), res=res1, type = 'normal')
  # ####BASIC#######
  # EnhancedVolcano(res1,
  #                 lab = rownames(res1),
  #                 x = 'log2FoldChange',
  #                 y = 'pvalue',
  #                 xlim = c(-8, 8),
  #                 title = paste(st,' Control versus Type 1 Diabetes'))
  # ########Modified cutoffs for pvalue and log2FC
  # EnhancedVolcano(res1,
  #                 lab = rownames(res1),
  #                 x = 'log2FoldChange',
  #                 y = 'pvalue',
  #                 xlim = c(-8, 8),
  #                 title =  paste(st,' Control versus Type 1 Diabetes'),
  #                 pCutoff = 10e-16,
  #                 FCcutoff = 1.5,
  #                 pointSize = 3.0,
  #                 labSize = 3.0)
  # 
  # #####################adjust color for point shading#####
  # EnhancedVolcano(res1,
  #                 lab = rownames(res1),
  #                 x = 'log2FoldChange',
  #                 y = 'pvalue',
  #                 xlim = c(-8, 8),
  #                 title =  paste(st,' Control versus Type 1 Diabetes'),
  #                 pCutoff = 10e-16,
  #                 FCcutoff = 1.5,
  #                 pointSize = 3.0,
  #                 labSize = 3.0,
  #                 col=c('black', 'black', 'black', 'red3'),
  #                 colAlpha = 1)
}




# #and lets make a new phyloseq object for each sampletype
 ASV_physeq_core_Introitus<-subset_samples(ASV_physeq_core, SampleType==" Introitus")
# ASV_physeq_core_Vagina<-subset_samples(ASV_physeq_core, SampleType=="Vagina")
# ASV_physeq_core_Cervix<-subset_samples(ASV_physeq_core, SampleType=="Cervix")
# ASV_physeq_core_Stool<-subset_samples(ASV_physeq_core, SampleType=="Stool")
# ASV_physeq_core_Anus<-subset_samples(ASV_physeq_core, SampleType=="Anus")
# ASV_physeq_core_Ear<-subset_samples(ASV_physeq_core, SampleType=="Ear")
# 

#now lests run each sampletype throught the previous function and identify 
#genus level taxa that are statistically significantly different between T1D and control
analyze_volcano(ASV_physeq_core_Introitus)
analyze_volcano(ASV_physeq_core_Vagina)
analyze_volcano(ASV_physeq_core_Cervix)
analyze_volcano(ASV_physeq_core_Stool)
analyze_volcano(ASV_physeq_core_Anus)
analyze_volcano(ASV_physeq_core_Ear)

