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
filtered_forward_reads <- paste0(samples, ".fwd.filt.fastq")
filtered_reverse_reads <- paste0(samples, ".rev.filt.fastq")                         
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                              rm.phix=TRUE, minLen=175, truncLen=c(220,200),multithread = 48,compress = FALSE)
#path <-getwd()
class(filtered_out) # matrix
dim(filtered_out) # 20 2
plotQualityProfile(filtered_reverse_reads[17:20], aggregate = T)
plotQualityProfile(filtered_forward_reads[17:20], aggregate = T)
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


#merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,derep_reverse, trimOverhang=TRUE, minOverlap=25,verbose = T)# this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,derep_reverse, trimOverhang=TRUE, minOverlap=100, verbose = T )
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
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)  # and now we can call the plot_richness() function on our phyloseq object

ASV_physeq


plot_bar(ASV_physeq,x="Patient",y="Abundance",fill = "Phylum")
sample_variables(ASV_physeq)

plot_bar(ASV_physeq,x="Host",y="Abundance",fill = "Phylum")


plot_bar(ASV_physeq,x="disease",y="RODZAJ_PORODU",fill = "Phylum")

plot_bar(ASV_physeq,x="sample.id",y="RODZAJ_PORODU",fill = "Phylum")


#BiocManager::install("microbiome")
library(microbiome)

ASV_physeq_rel<-transform(ASV_physeq,transform="compositional")
plot_bar(ASV_physeq_rel,x="Patient",y="Abundance",fill = "Phylum")

plot_bar(ASV_physeq_rel,x="sample.id",y="Abundance",fill = "Phylum")
#materials in the infant
infant<-subset_samples(ASV_physeq,Host=="Infant")
infant
infant_rel<-transform(infant,transform="compositional")
plot_bar(infant_rel,x="sample_id", "Abundance", fill = "Phylum",facet_grid = "SampleType")
#all materials from the infant together
plot_bar(infant_rel,x="sample_id", "Abundance", fill = "Phylum")
# to divide into the ear and the stool samples in infants: faced_grid
plot_bar(infant_rel,x="sample_id", "Abundance", fill = "Phylum",facet_grid = "SampleType")

#materials in the mothers
mother<-subset_samples(ASV_physeq,Host=="Mother")
mother
mother_rel<-transform(mother,transform = "compositional")
plot_bar(mother_rel,x="sample_id", "Abundance", fill ="Phylum")                  
