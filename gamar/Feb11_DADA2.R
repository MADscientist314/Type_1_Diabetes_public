#install.packages("devtools")
#library("devtools")
#devtools::install_github("benjjneb/dada2", ref="v1.14") # change the ref argument to get other versions
library('dada2') #load the library
getwd() #get the working directory
setwd('/home/gamar/Type_1_Diabetes_Project/paired')#set the working directory 
samples <- scan("list.good.txt", what="character")
samples
# one holding the file names of all the forward reads
forward_reads <- paste0(samples, ".fwd.paired.fastq")
# and one with the reverse
reverse_reads <- paste0(samples, ".rev.paired.fastq")
length(forward_reads)
# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads <- paste0(samples, ".175.fwd.filt.fastq")
filtered_reverse_reads <- paste0(samples, ".175.rev.filt.fastq")                         
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                              rm.phix=TRUE, minLen=150, trimLeft = c(10,10), trimRight = c(25,50),multithread = 48,compress = FALSE)
#path <-getwd()
class(filtered_out) # matrix
dim(filtered_out) # 20 2
plotQualityProfile(filtered_forward_reads[1:20], aggregate = T)
plotQualityProfile(filtered_reverse_reads[1:20], aggregate = T)
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=48) 
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=48) 
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
#filt_samples<-gsub(pattern = ".fwd.filt.fastq",replacement = "",x = filtered_forward_reads)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples
dada_forward <- dada(derep_forward, err=err_forward_reads, verbose = TRUE,pool=TRUE, multithread=48)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, verbose = TRUE, pool=TRUE, multithread=48) 

summary_tab
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,derep_reverse, trimOverhang=TRUE, minOverlap=5,verbose = T)# this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
#merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,derep_reverse, trimOverhang=TRUE, minOverlap=100, verbose = T )
class(merged_amplicons) # list
length(merged_amplicons) # 20 elements in this list, one for each of our samples
names(merged_amplicons) # the names() function gives us the name of each element of the list seqtab <- makeSequenceTable(merged_amplicons)




seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix
dim(seqtab) # 20 2521
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) # Howd many bimeras out of how many input sequences.
sum(seqtab.nochim)/sum(seqtab) # what is the abundance rate?  0.9931372 # good, we barely lost any in terms of abundance# set a little function
getN <- function(x) sum(getUniques(x))# making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
summary_tab
summary_tab
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz",verbose = T, multithread=TRUE, tryRC=T) # giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {asv_headers[i] <- paste(">ASV", i, sep="_")}# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)# tax table:
asv_tax <- summary_tabtaxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

asv_tab

library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
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
mother <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)  # and now we can call the plot_richness() function on our phyloseq object

mother


plot_bar(mother,x="Patient",y="Abundance",fill = "Phylum")
sample_variables(mother)

plot_bar(mother,x="Host",y="Abundance",fill = "Phylum")


plot_bar(mother,x="disease",y="RODZAJ_PORODU",fill = "Phylum",facet_grid = "SampleType")
plot_bar(mother_rel,x="disease",y="RODZAJ_PORODU",fill = "Phylum",facet_grid = "SampleType")
plot_bar(mother_rel,x="disease",y="sample_id",fill = "Phylum",facet_grid = "SampleType")
plot_bar(mother,x="disease",y="Abundance",fill = "Phylum")

plot_bar(mother,x="sample.id",y="RODZAJ_PORODU",fill = "Phylum")


#BiocManager::install("microbiome")
library(microbiome)

mother_rel<-transform(mother,transform="compositional")
plot_bar(mother_rel,x="Patient",y="Abundance",fill = "Phylum")


#blad w zapisie sample_id, nie powinno być sample.id
plot_bar(mother_rel,x="sample.id",y="Abundance",fill = "Phylum")

#materials in the infant
infant<-subset_samples(mother,Host=="Infant")
infant
infant_rel<-transform(infant,transform="compositional")
plot_bar(infant_rel,x="sample_id", "Abundance", fill = "Phylum",facet_grid = "SampleType")
#all materials from the infant together
plot_bar(infant_rel,x="sample_id", "Abundance", fill = "Phylum")

plot_bar(infant_rel,x="Patient", y="Abundance", fill = "Phylum")
plot_bar(infant_rel,x="Patient", y="Abundance", fill = "Class")
plot_bar(infant_rel,x="Patient", y="Abundance", fill = "Order")
plot_bar(infant_rel,x="Patient", y="Abundance", fill = "Family")
plot_bar(infant_rel,x="Patient", y="Abundance", fill = "Genus")
#podzial na 4 materiały u matki z uwzględnieniem abundance
plot_bar(infant_rel,x="SampleType", y="Abundance", fill = "Order")
#podzial na 4 materiały u matki z uwzględnieniem abundance
plot_bar(mother_rel,x="SampleType", y="Abundance", fill = "Order")



# to divide into the ear and the stool samples in infants: faced_grid
plot_bar(infant_rel,x="sample_id", y="Abundance", fill = "Phylum",facet_grid = "SampleType")

#materials in the mothers
mother<-subset_samples(mother,Host=="Mother")
mother
mother_rel<-transform(mother,transform = "compositional")
plot_bar(mother_rel,x="sample_id", "Abundance", fill ="Phylum")

# to divide into the m1, m2, 3 and o samples in mothers: faced_grid

#dlaczego tutaj mi nie wychodzila figure? BO PRZED "ABUNDANCE' NIE BYŁO y=
plot_bar(mother_rel,x="sample_id", "Abundance", fill = "Phylum",faced_grid = "SampleType")
#tutaj dodałam druga zmienna y= i poszlo, TO Z PODZIALEM NA MATERIALY U MATEK, BO MOTHER_REL
plot_bar(mother_rel,x="sample_id", y="Abundance", fill = "Phylum",facet_grid = "SampleType")
#tutaj zaczernilo
plot_bar(mother_rel,x="Host", y="SampleType", fill ="Phylum",facet_grid = "SampleType")

plot_bar(mother_rel,x="Disease", y="SampleType", fill ="Phylum",facet_grid = "SampleType")

mother
sample_variables(mother)
#sprawdzilam variables, zamiast Disease powinno byc disease:
plot_bar(mother_rel,x="disease", y="SampleType", fill ="Phylum",facet_grid = "SampleType")
#tutaj na osi X T1D vs. kontrola, na osi Y abundance a wszystko z podzialem na 4 materialy (SampleType) u matki:
plot_bar(mother_rel,x="disease", y="Abundance", fill ="Phylum",facet_grid = "SampleType")
#tutaj na osi X T1D vs. kontrola, na osi Y abundance a wszystko z podzialem na 2 materialy (SampleType) u dziecka:
plot_bar(infant_rel,x="disease", y="Abundance", fill ="Phylum",facet_grid = "SampleType")


# FEB 11 Microbiome COMPOSITION


library(microbiome)
library(dplyr)
data(dietswap)
transform <- microbiome::transform
pseq <- core(dietswap, detection = 50, prevalence = 50/100)

library(phyloseq)
pseq2 <-subset_samples(pseq, group == "DI" & nationality == "AFR" & timepoint.within.group == 1)
data(atlas1006)
pseq3 <- atlas1006 %>%
        subset_samples(DNA_extraction_method == "r") %>%
        aggregate_taxa(level = "Phylum") %>%
        microbiome::transform(transform = "compositional")


#composition heatmaps
p <- plot_composition(microbiome::transform(pseq, "compositional"),
                      plot.type = "heatmap",
                      sample.sort = "neatmap",
            otu.sort = "neatmap")
print(p)                      
#plot taxa prevalence
p0 <- subset_samples(atlas1006, DNA_extraction_method == "r")
p0 <- core(p0, detection =10, prevalence = 0)
plot_taxa_prevalence(p0, "Phylum", detection = 10)
print(p_prev)


install.packages("hrbrthemes")
install.packages("gcookbook")

themes_set(theme_bw(21))
p <-pseq3 %>%
  plot_composition(sample.sort = "Firmicutes", otu.sort = "abundance") +
  scale_fill_manual(values = default_colors("Phylum")[taxa(pseq3)])
print(p)

# Limit the analysis on core taxa and specific sample group
p <- plot_composition(pseq2,
                      taxonomic.level = "Genus",
                      sample.sort = "nationality",
                      x.label = "nationality") +
  guides(fill = guide_legend(ncol = 1)) +
  scale_y_percent() +
  labs(x = "Samples", y = "Relative abundance (%)",
       title = "Relative abundance data",
       subtitle = "Subtitle",
       caption = "Caption text.") + 
  theme_ipsum(grid="Y")
print(p) 


# Averaged by group
p <- plot_composition(pseq3,
                      average_by = "bmi_group", transform = "compositional")
print(p)
p <- NULL








library("RColorBrewer")
library('microbiomeutilities')
mother
plot_taxa_composition(x = mother, 
                      sample.sort = "SampleType", 
                      taxonomic.level = "Order", 
                      transform = "compositional", 
                      otu.sort = "abundance", 
                      palette = brewer.pal(12, "Paired"), 
                      x.label = "disease_SampleType", 
                      plot.type = "barplot", 
                      average_by = "disease_SampleType", 
                      verbose = TRUE)# + 
                      #scale_x_discrete("SampleType")
#instead of order will be Class
plot_taxa_composition(x = mother, 
                      sample.sort = "SampleType", 
                      taxonomic.level = "Class", 
                      transform = "compositional", 
                      otu.sort = "abundance", 
                      palette = brewer.pal(12, "Paired"), 
                      x.label = "disease_SampleType", 
                      plot.type = "barplot", 
                      average_by = "disease_SampleType", 
                      verbose = TRUE)# + 
#scale_x_discrete("SampleType")

#instead of order will be Family
plot_taxa_composition(x = mother, 
                      sample.sort = "SampleType", 
                      taxonomic.level = "Family", 
                      transform = "compositional", 
                      otu.sort = "abundance", 
                      palette = brewer.pal(12, "Paired"), 
                      x.label = "disease_SampleType", 
                      plot.type = "barplot", 
                      average_by = "disease_SampleType", 
                      verbose = TRUE)

infant

#instead of Order will be Class
plot_taxa_composition(x = infant, 
                      sample.sort = "SampleType", 
                      taxonomic.level = "Class", 
                      transform = "compositional", 
                      otu.sort = "abundance", 
                      palette = brewer.pal(12, "Paired"), 
                      x.label = "disease_SampleType", 
                      plot.type = "barplot", 
                      average_by = "disease_SampleType", 
                      verbose = TRUE)
#plot_taxa_boxplot

plot_taxa_boxplot(x = mother_rel,taxonomic.level = "Phylum", top.otu = 5, VariableA = "disease_SampleType",title = "title", color ="Paired" )

#Kingdom
plot_taxa_boxplot(x = mother_rel,taxonomic.level = "Kingdom", top.otu = 5, VariableA = "disease_SampleType",title = "title", color ="Paired" )


#12feb2020

#infant only

plot_taxa_boxplot(x = infant_rel,taxonomic.level = "Genus", top.otu = 25, VariableA = "disease_SampleType",title = "title", color ="Paired" )

infant_trim<-prune_taxa(taxa_sums(infant) > 10, infant)

taxa_names(infant_trim)
#dla ASV_20 tylko, zatem wybrany jakiś
plot_select_taxa(x = infant_trim, select.taxa = "ASV_20", variableA = "disease_SampleType", palette = brewer.pal(12, "Paired"), plot.type = "boxplot")
plot_select_taxa(x = infant_trim, select.taxa = "ASV_24", variableA = "disease_SampleType", palette = brewer.pal(12, "Paired"), plot.type = "boxplot")
mother
taxa_names(mother_trim)
plot_select_taxa(x = mother_trim, select.taxa = "ASV_956", variableA = "disease_SampleType", palette = brewer.pal(12, "Paired"), plot.type = "boxplot")
plot_select_taxa(x = mother_trim, select.taxa = "ASV_986", variableA = "disease_SampleType", palette = brewer.pal(12, "Paired"), plot.type = "boxplot")

#inny
plot_select_taxa(x = infant_trim, select.taxa = "ASV_149", variableA = "disease_SampleType", palette = brewer.pal(12, "Paired"), plot.type = "boxplot")

mother
A<-aggregate_taxa(mother, level = "Genus")
B<-tax_glom(mother, taxrank = "Genus",NArm = F)

plot_taxa_composition(x = A, 
                      sample.sort = "SampleType", 
                      taxonomic.level = "Order", 
                      transform = "compositional", 
                      otu.sort = "abundance", 
                      palette = brewer.pal(12, "Paired"), 
                      x.label = "disease_SampleType", 
                      plot.type = "barplot", 
                      average_by = "disease_SampleType", 
                      verbose = TRUE)# + 
#scale_x_discrete("SampleType")

#boxplot instead of barplot
plot_taxa_composition(x = A, 
                      sample.sort = "SampleType", 
                      taxonomic.level = "Order", 
                      transform = "compositional", 
                      otu.sort = "abundance", 
                      palette = brewer.pal(12, "Paired"), 
                      x.label = "disease_SampleType", 
                      plot.type = "boxplot", 
                      average_by = "disease_SampleType", 
                      verbose = TRUE)
#Class intead of Order, boxplot
plot_taxa_composition(x = A, 
                      sample.sort = "SampleType", 
                      taxonomic.level = "Class", 
                      transform = "compositional", 
                      otu.sort = "abundance", 
                      palette = brewer.pal(12, "Paired"), 
                      x.label = "disease_SampleType", 
                      plot.type = "errorbar", 
                      average_by = "disease_SampleType", 
                      verbose = TRUE)

#Class intead of Order, barplot
plot_taxa_composition(x = A, 
                      sample.sort = "SampleType", 
                      taxonomic.level = "Class", 
                      transform = "compositional", 
                      otu.sort = "abundance", 
                      palette = brewer.pal(12, "Paired"), 
                      x.label = "disease_SampleType", 
                      plot.type = "barplot", 
                      average_by = "disease_SampleType", 
                      verbose = TRUE)

plot_taxa_composition(x = B, 
                      sample.sort = "SampleType", 
                      taxonomic.level = "Class", 
                      transform = "compositional", 
                      otu.sort = "abundance", 
                      palette = brewer.pal(12, "Paired"), 
                      x.label = "disease_SampleType", 
                      plot.type = "barplot", 
                      average_by = "disease_SampleType", 
                      verbose = TRUE)# + 
#scale_x_discrete("SampleType")


plot_taxa_composition(x = mother, 
                      sample.sort = "SampleType", 
                      taxonomic.level = "Order", 
                      transform = "compositional", 
                      otu.sort = "abundance", 
                      palette = brewer.pal(12, "Paired"), 
                      x.label = "disease_SampleType", 
                      plot.type = "barplot", 
                      average_by = "disease_SampleType", 
                      verbose = TRUE)# + 
#scale_x_discrete("SampleType")
#czestosci prob materialow u matek - 4 podgrupy 
plot_frequencies(x = meta(mother), Groups = "disease", Factor = "SampleType")

mother_trim<-prune_taxa(taxa_sums(mother) > 10, mother)
mother_trim<-prune_samples(sample_sums(mother) > 10, mother)
mother_ord<-ordinate(physeq = mother_trim, method = "CCA", distance = "bray")
plot_ordination(physeq = mother_trim,ordination = mother_ord, color = "disease_SampleType")+ scale_color_brewer(palette = "Paired")
    +theme_minimal()
    +geom_polygon(alpha=0.1)+facet_wrap("SampleType")

#czestosci prob materialow u infants - 2 podgrupy 
plot_frequencies(x = meta(infant), Groups = "disease", Factor = "SampleType")


#przepisane powyzej na infant
infant_trim<-prune_taxa(taxa_sums(infant) > 10, infant)
infant_trim<-prune_samples(sample_sums(infant) > 10, infant)
infant_ord<-ordinate(physeq = infant_trim, method = "DCA", distance = "bray")
plot_ordination(physeq = infant_trim, ordination = infant_ord, color = "disease_SampleType")+ scale_color_brewer(palette = "Paired")
  +theme_minimal()+geom_polygon(alpha=0.1)+facet_wrap("SampleType")


#other configurations:

infant_trim<-prune_taxa(taxa_sums(infant) > 5, infant)
infant_trim<-prune_samples(sample_sums(infant) > 5, infant)
infant_ord<-ordinate(physeq = infant_trim, method = "MDS", distance = "bray")
plot_ordination(physeq = infant_trim, ordination = infant_ord, color = "disease_SampleType")+ scale_color_brewer(palette = "Paired")
    +theme_minimal()+geom_polygon(alpha=0.1)+facet_wrap("SampleType")
  
#Feb13 alfa-diversity

p <- plot_alpha_diversities(ASV_physeq_trim,
                            type = "evenness", #diversity, dominancem or eveness
                            index.val = "all",
                            plot.type = "stripchart",# options "stripchart", "boxplot", "violin"
                            variableA = "disease_SampleType", #Variable to be checked
                            palette = "Paired") # any palette in Rcolor brewer)
print(p)

plot_richness(ASV_physeq_trim)
plot_richness(ASV_physeq_trim, color="disease_SampleType", measures=c("Chao1", "Shannon"))

ASV_physeq_merge = merge_samples(ASV_physeq, "disease_SampleType")
p = plot_richness(ASV_physeq_merge, color="disease_SampleType", measures=c("Chao1", "Shannon"))
print(p)

#some STAT, 
#is there a significant difference between the disese and control?

ASV_physeq_trim<- format_to_besthit(ASV_physeq)
ASV_physeq_rel <-  microbiome::transform(ASV_physeq_trim, "compositional")
otu <- abundances(ASV_physeq_rel)
meta <- meta(ASV_physeq_rel)


#ear<-subset_samples(ASV_physeq_trim,SampleType=="Ear")
#ear
#ear_rel<-transform(ear,transform = "compositional")

#otu <- abundances(ear_rel)
#meta <- meta(ear_rel)


library(vegan)
#> Loading required package: permute
#> Loading required package: lattice
#> This is vegan 2.5-3
#> 
#> Attaching package: 'vegan'
#> The following object is masked from 'package:microbiome':
#> 
#>     diversity
permanova <- adonis(t(otu) ~ disease,
                    data = meta, permutations = 999, method = "bray")# P-value
permanova
print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])
#> [1] NAcoef <- coefficients(permanova)["disease_SampleType", ]
#top.coef <- coef[rev(order(abs(coef)))[1:3]]top.coef.df <- as.data.frame(top.coef)
#my_taxa <- c(rownames(top.coef.df))p <- plot_select_taxa(psf2.rel, my_taxa, "disease_SampleType", "Paired", plot.type = "boxplot")
#> An additonal column Sam_rep with sample names is created for reference purpose
print(p)



