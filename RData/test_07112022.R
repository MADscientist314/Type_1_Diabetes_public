library(dada2)
library(phyloseq)
library(microbiome)
library(tidyverse)
library(phangorn)
library(DECIPHER)
library(dada2)
library(phyloseq)
library(microbiome)
library(tidyverse)
setwd("k:/github/Type_1_Diabetes_public/gamar/")
getwd()
count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")[ , -c(1:4)]
dim(count_tab)
dim(tax_tab)
tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))
dim(tax_tab)
sample_info_tab <- read.table("k:/github/type1_diabetes_project/sample_core.csv", header=T, row.names=1,
                              check.names=F, sep=",")
dim(sample_info_tab)#[1] 527 159

ASV_physeq<-phyloseq(otu_table(count_tab,taxa_are_rows = T),
            tax_table(tax_tab),
            sample_data(sample_info_tab))

ASV_physeq_core<-core(ASV_physeq,detection = 0,prevalence = 3/100)
ASV_physeq_core #otu_table()   OTU Table:         [ 2779 taxa and 532 samples ]


write.table(taxa_names(ASV_physeq_core),"ASV_list.tsv",sep="\t",row.names = F,quote = F)

dim(count_tab)#[1] 1609  532
dim(tax_tab)#[1] 1609    6
dim(sample_info_tab)#[1] 527 159
colnames(count_tab)
sample_info_tab
#make a list of the taxa names in the ASV_physeq_core object
getwd()
# import the fasta file with all ASVs
dna<-readDNAStringSet(filepath = "k:/github/Type_1_Diabetes_public/gamar/ASVs.fa",
                      format = "fasta",
                      nrec = -1L,
                      use.names = T)
dna


# This propagates to the tip labels of the tree
alignment <-AlignSeqs(dna,
                      processors = 20,
                      verbose = T, 
                      anchor=NA,
                      iterations = 2,
                      refinements = 1,)

# 
# run it again but longer this time and look at the rana and protein levels
#source("alignment.R")
#
summary(alignment)
summary(alignment2)
saveRDS(alignment2,"alignment2.RDS")
alignment<-alignment2

library(phangorn)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
class(phang.align)
dm <- dist.ml(x = phang.align)
class(dm)

treeNJ <- NJ(dm) # Note, tip order != sequence order
class(treeNJ)

fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(object = fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR
treeNJ$tip.label
ASV_physeq_core
#Combine data into the existing phyloseq object
phy_tree(ASV_physeq)<-phy_tree(fitGTR$tree)


#optional troubleshooting code that renames the tree taxa to match the phylo obj 
# phytree<-phy_tree(fitGTR$tree)
# taxa_names(phytree)<-taxa_names(ASV_physeq)
# phy_tree(ASV_physeq)<-phytree
########## Construct phylogenetic tree##############

ASV_physeq

#save the origninal in case you mess it up
ASV_physeq
saveRDS(ASV_physeq_core,"ASV_physeq_core_with_tree_07112022.RDS")
saveRDS(ASV_physeq,"ASV_physeq_with_tree_07112022.RDS")

