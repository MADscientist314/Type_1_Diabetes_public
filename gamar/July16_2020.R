sample_corejuly15<-read.csv("sample_coreJuly15.txt", header = T, row.names = 1, sep = "\t")
getwd()
setwd("Type1")
setwd("/home/gamar/Type_1_Diabetes_Project/paired")
sam<-sample_data(sample_corejuly15)
sample_data(ASV_physeq_core)<-sam
ASV_physeq_core
sample_data(ASV_physeq_core)


#CREATING A BOXPLOT FOR THE TAXA ASSOCIATED WITH EACH VARIABLE
#1. Subset the phyloseq_pbject by SampleType
ASV_physeq_core_genus<-tax_glom(ASV_physeq_core,taxrank = "Genus",NArm = T )
ASV_physeq_core_Introitus<-subset_samples(ASV_physeq_core_genus, SampleType==" Introitus")
ASV_physeq_core_Vagina<-subset_samples(ASV_physeq_core_genus, SampleType=="Vagina")
ASV_physeq_core_Cervix<-subset_samples(ASV_physeq_core_genus, SampleType=="Cervix")
ASV_physeq_core_Stool<-subset_samples(ASV_physeq_core_genus, SampleType=="Stool")
ASV_physeq_core_Anus<-subset_samples(ASV_physeq_core_genus, SampleType=="Anus")
ASV_physeq_core_Ear<-subset_samples(ASV_physeq_core_genus, SampleType=="Ear")

library(microbiome)
library(microbiomeutilities)

analyze_Maternal<-function(filename){
  ps.f2 <- format_to_besthit(filename)
  psf2.rel <- microbiome::transform(ps.f2, "compositional")
  otu <- abundances(psf2.rel)
  meta <- meta(psf2.rel)
  permanova <- adonis(t(otu) ~ disease*SGA_AGA_LGA,data = meta, permutations = 999, method = "bray")
  filename_type<-as.character(sample_data(psf2.rel)$SampleType[1])
  print(filename_type)
  print(permanova)
  coef <- coefficients(permanova)["disease1", ]
  #rownames(coefficients(permanova))
  top.coef <- coef[rev(order(abs(coef)))[1:3]]
  top.coef.df <- as.data.frame(top.coef)
  my_taxa <- c(rownames(top.coef.df))
  p <- plot_select_taxa(psf2.rel, my_taxa,variableA =  "disease", palette = "Set1", plot.type = "boxplot")
  #> An additonal column Sam_rep with sample names is created for reference purpose
  top.coef <- coef[rev(order(abs(coef)))[1:20]]
  par(mar = c(3, 14, 2, 1))
  b<-barplot(sort(top.coef),horiz = T, las = 1, main = paste("Top taxa", filename_type, sep = " "))
  #> An additonal column Sam_rep with sample names is created for reference purpose
  print(p)
}

#function

analyze_Infant<-function(filename){
  
  ps.f2 <- format_to_besthit(filename)
  
  psf2.rel <- microbiome::transform(ps.f2, "compositional")
  
  otu <- abundances(psf2.rel)
  
  meta <- meta(psf2.rel)
  
  permanova <- adonis(t(otu) ~ disease*Delivery*Antibiotics,data = meta, permutations = 999, method = "bray")
  
  filename_type<-as.character(sample_data(psf2.rel)$SampleType[1])
  
  print(filename_type)
  
  print(permanova)
  
  coef.deliv <- coefficients(permanova)["Delivery1", ]
  
  coef.disease <- coefficients(permanova)["disease1", ]
  
  rownames(coefficients(permanova))
  
  top.coef.disease <- coef[rev(order(abs(coef)))[1:3]]
  
  top.coef.disease.df <- as.data.frame(top.coef.disease)
  
  top.coef.deliv <- coef[rev(order(abs(coef)))[1:3]]
  
  top.coef.deliv.df <- as.data.frame(top.coef.deliv)
  
  my_taxa.deliv <- c(rownames(top.coef.deliv.df))
  
  my_taxa.disease <- c(rownames(top.coef.disease.df))
  
  p <- plot_select_taxa(psf2.rel, my_taxa.deliv,variableA =  "Delivery", palette = "Set2", plot.type = "boxplot")
  
  q <- plot_select_taxa(psf2.rel, my_taxa.disease,variableA =  "disease", palette = "Set1", plot.type = "boxplot")
  
  top.coef.disease <- coef[rev(order(abs(coef.disease)))[1:20]]
  
  par(mar = c(3, 14, 2, 1))
  
  b<-barplot(sort(top.coef.disease), horiz = T, las = 1, main = paste("Top taxa disease", filename_type, sep = " "))
  
  #> An additonal column Sam_rep with sample names is created for reference purpose
  
  top.coef.deliv <- coef[rev(order(abs(coef.deliv)))[1:20]]
  
  par(mar = c(3, 14, 2, 1))
  
  b<-barplot(sort(top.coef.deliv), horiz = T, las = 1, main = paste("Top taxa deliv", filename_type, sep = " "))
  
  
  
  print(p)
  
  print(q)
  
}


#MATERNAL SAMPLES
analyze_Maternal(ASV_physeq_core_Introitus)
analyze_Maternal(ASV_physeq_core_Vagina)
analyze_Maternal(ASV_physeq_core_Cervix)
analyze_Maternal(ASV_physeq_core_Anus)
#INFANT SAMPLES
analyze_Infant(ASV_physeq_core_Stool)
analyze_Infant(ASV_physeq_core_Ear)

library(vegan)

plot_taxa_boxplot(x = ASV_physeq_core_Ear, 
                                                     taxonomic.level = "Genus", 
                                                     top.otu = 3, 
                                                     VariableA = "Delivery", 
                                                     title = "Ear Canal Delivery")
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")

#individual taxa 
#differential abundance testing



#############################################################

#DESeq2: differential abundance testing for sequencing data

############################################################

#DESeq2 analysis can accommodate those particular assumptions about sequencing data.

# Example data

library(microbiome)

library(microbiomeutilities)

library(tibble)

# Start by converting phyloseq object to deseq2 format

library(DESeq2)

ASV_physeq_core<-format_to_besthit(ASV_physeq_core)

ds2 <- phyloseq_to_deseq2(ASV_physeq_core, ~ disease + SampleType + disease*SampleType)

# Run DESeq2 analysis (all taxa at once!)

dds <- DESeq(ds2, sfType = 'poscounts', parallel = TRUE)

# Investigate results

res <- results(dds)

deseq.results <- as.data.frame(res)

df <- deseq.results

df$taxon <- rownames(df)

df <- df %>% arrange(log2FoldChange, padj)

# Print the results; flitered and sorted by pvalue and effectsize

library(knitr)

df <- df %>% filter(pvalue < 0.05 & log2FoldChange > 1.5) %>%
  
  arrange(pvalue, log2FoldChange)

kable(df, digits = 5)

#For comparison purposes, assess significances and effect sizes based on Wilcoxon test.

d<-ASV_physeq_core

test.taxa <- taxa(d)

pvalue.wilcoxon <- c()

foldchange <- c()

for (taxa in test.taxa) {
  
  # Create a new data frame for each taxonomic group
  
  df <- data.frame(Abundance = abundances(d)[taxa,],
                   
                   Log10_Abundance = log10(1 + abundances(d)[taxa,]),
                   
                   Group = meta(d)$nationality)
  
  # Calculate pvalue and effect size (difference beween log means)       
  
  pvalue.wilcoxon[[taxa]] <- wilcox.test(Abundance ~ Group, data = df)$p.value
  
  foldchange[[taxa]] <- coef(lm(Log10_Abundance ~ Group, data = df))[[2]]
  
}

# Correct p-values for multiple testing

pvalue.wilcoxon.adjusted <- p.adjust(pvalue.wilcoxon)

par(mfrow = c(1,2))

plot(deseq.results$padj, pvalue.wilcoxon.adjusted,
     
     xlab = "DESeq2 adjusted p-value",
     
     ylab = "Wilcoxon adjusted p-value",
     
     main = "P-value comparison")

abline(v = 0.05, h = 0.05, lty = 2)

plot(deseq.results$log2FoldChange, foldchange, 
     
     xlab = "DESeq2",
     
     ylab = "Linear model",
     
     main = "Effect size comparison")

abline(0,1)

