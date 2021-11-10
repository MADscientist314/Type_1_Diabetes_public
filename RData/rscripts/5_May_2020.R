#CREATING A BOXPLOT FOR THE TAXA ASSOCIATED WITH EACH VARIABLE
#1. Subset the phyloseq_pbject by SampleType
ASV_physeq_core_genus<-tax_glom(ASV_physeq_core,taxrank = "Genus",NArm = T )
ASV_physeq_core_Introitus<-subset_samples(ASV_physeq_core_genus, SampleType==" Introitus")
ASV_physeq_core_Vagina<-subset_samples(ASV_physeq_core_genus, SampleType=="Vagina")
ASV_physeq_core_Cervix<-subset_samples(ASV_physeq_core_genus, SampleType=="Cervix")
ASV_physeq_core_Stool<-subset_samples(ASV_physeq_core_genus, SampleType=="Stool")
ASV_physeq_core_Anus<-subset_samples(ASV_physeq_core_genus, SampleType=="Anus")
ASV_physeq_core_Ear<-subset_samples(ASV_physeq_core_genus, SampleType=="Ear")

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

plot_taxa_boxplot(x = ASV_physeq_core_Ear, 
                  taxonomic.level = "Genus", 
                  top.otu = 3, 
                  VariableA = Delivery, 
                  title = "Ear Canal Delivery")
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
