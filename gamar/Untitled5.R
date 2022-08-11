library(dada2)
library(phyloseq)
library(microbiome)
library(tidyverse)
ASV_physeq_core
keep<-taxa_names(ASV_physeq_core)
library(phangorn)
library(DECIPHER)
library(dada2)
library(phyloseq)
library(microbiome)
library(tidyverse)

keep<-taxa_names(ASV_physeq_core)

# 
# 
# count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1,
#                         check.names=F, sep="\t")[ , -c(1:4)]
# 
# tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T,
#                                 row.names=1, check.names=F, sep="\t"))
# 
# sample_info_tab <- read.table("meta.tsv", header=T, row.names=1,
#                               check.names=F, sep="\t")
dna<-readDNAStringSet(filepath = "picrust2/ASVs.fa",
                      format = "fasta",
                      nrec = -1L,
                      use.names = T)
#extract only the reads in ASV_physeq_core
dna2<-dna[keep]

dna2
# This propagates to the tip labels of the tree
alignment <-AlignSeqs(dna2,
                      processors = 48,
                      verbose = T, 
                      anchor=NA,
                      iterations = 2,
                      refinements = 1,)

# 
# run it again but longer this time and look at the rana and protein levels
# alignment2<-AlignSeqs(DNAStringSet(seqs),
#                       processors = 20,
#                       verbose = T,
#                       useStructures = T,
#                       guideTree = alignment, 
#                       anchor=NA,
#                       iterations = 10,
#                       refinements = 5,)
# 
fa<-readDNAStringSet(filepath = "ASVs.fa",format = "fasta",nrec = -1L,use.names = T)
fa
ASV_physeq_core
library(phangorn)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")

class(phang.align)

dm <- dist.ml(x = phang.align)

class(dm)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(object = fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))


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
saveRDS(ASV_physeq,"ASV_physeq.RDS")

