library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library(microbiome)
library(microbiomeutilities)


count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")
tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
sample_info_tab <- read.table("10_Feb_2020_metadata_MG.csv", header=T, row.names=1,
                              check.names=F, sep=",")
# and setting the color column to be of type "character", which helps later
sample_info_tab # to take a peek

count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
sample_info_tab_phy <- sample_data(sample_info_tab)
ASV_physeq<- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)  # and now we can call the plot_richness() function on our phyloseq object
#to remove 2,4,5,6

ASV_physeq<-subset_samples(ASV_physeq, patient_disease!="2_T1D")
ASV_physeq<-subset_samples(ASV_physeq, patient_disease!="4_T1D")
ASV_physeq<-subset_samples(ASV_physeq, patient_disease!="5_T1D")
ASV_physeq<-subset_samples(ASV_physeq, patient_disease!="6_T1D")

ASV_physeq


ASV_physeq <- subset_taxa(ASV_physeq, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
  
  
  

ASV_physeq_trim<-prune_samples(sample_sums(ASV_physeq)>0, ASV_physeq)




ASV_physeq_trim
ASV_physeq_trim<- format_to_besthit(ASV_physeq_trim)
ASV_physeq_rel <-  microbiome::transform(ASV_physeq_trim, "compositional")
otu <- abundances(ASV_physeq_rel)

meta <- meta(ASV_physeq_rel)




library(vegan)
#> Loading required package: permute
#> Loading required package: lattice
#> This is vegan 2.5-3
#> 
#> Attaching package: 'vegan'
#> The following object is masked from 'package:microbiome':
#> 
#>     diversity
permanova <- adonis(t(otu) ~ disease_SampleType,
                    data = meta, permutations = 999, method = "bray")# P-value
permanova


#ear
ear<-subset_samples(ASV_physeq_trim,SampleType=="Ear")
ear
ear_rel<-transform(ear,transform = "compositional")
otu <- abundances(ear_rel)
meta <- meta(ear_rel)
ear_permanova <- adonis(t(otu) ~ disease_SampleType,
                        data = meta, permutations = 999, method = "bray")# P-value


ear_permanova

#cervix
Cervix<-subset_samples(ASV_physeq_trim,SampleType=="Cervix")
Cervix
Cervix_rel<-transform(Cervix,transform = "compositional")
otu <- abundances(Cervix_rel)
meta <- meta(Cervix_rel)
Cervix_permanova <- adonis(t(otu) ~ disease_SampleType,
                           data = meta, permutations = 999, method = "bray")# P-value
Cervix_permanova

#stool
Stool<-subset_samples(ASV_physeq_trim,SampleType=="Stool")
Stool
Stool_rel<-transform(Stool,transform = "compositional")
otu <- abundances(Stool_rel)
meta <- meta(Stool_rel)
Stool_permanova <- adonis(t(otu) ~ disease_SampleType,
                          data = meta, permutations = 999, method = "bray")# P-value
Stool_permanova

#vagina
Vagina<-subset_samples(ASV_physeq_trim,SampleType=="Vagina")
Vagina
Vagina_rel<-transform(Vagina,transform = "compositional")
otu <- abundances(Vagina_rel)

meta <- meta(Vagina_rel)
Vagina_permanova <- adonis(t(otu) ~ disease_SampleType,
                           data = meta, permutations = 999, method = "bray")# P-value
Vagina_permanova

#Anus
Anus<-subset_samples(ASV_physeq_trim,SampleType=="Anus")
Anus
Anus_rel<-transform(Anus,transform = "compositional")
otu <- abundances(Anus_rel)
meta <- meta(Anus_rel)
Anus_permanova <- adonis(t(otu) ~ disease_SampleType,
                         data = meta, permutations = 999, method = "bray")# P-value
Anus_permanova

#introitus = vestibule
Vestibule<-subset_samples(ASV_physeq_trim,SampleType=="Anus")
Vestibule
Vestibule_rel<-transform(Vestibule,transform = "compositional")
otu <- abundances(Vestibule_rel)

meta <- meta(Vestibule_rel)

Vestibule_permanova <- adonis(t(otu) ~ disease,
                              data = meta, permutations = 999, method = "bray")# P-value

Vestibule_permanova

#SHINY
save.image(file = "ASV_physeq2.RData")
saveRDS(file = "ASV_physeq.rds", ASV_physeq_trim)
getwd()
#install.packages("shiny")
library(shiny)
shiny::runGitHub("shiny-phyloseq","joey711")


#Feb18


ASV_physeq

write.csv((otu_table(ASV_physeq)), "otu.tsv")
write.csv((tax_table(ASV_physeq)), "tax.tsv")
write.csv((sample_data(ASV_physeq)), "sample_data.tsv")

