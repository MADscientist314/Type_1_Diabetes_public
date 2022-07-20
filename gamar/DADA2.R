install.packages("devtools", dependencies = TRUE)
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.14") # change the ref argument to get other versions

library('dada2')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2",)
library('dada2')
getwd()
list.files()
## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("../samples", what="character")  # one holding the file names of all the forward reads
forward_reads <- paste0(samples, "_R1_001.paired.fastq.gz")
# and one with the reverse
reverse_reads <- paste0(samples, "_R2_001.paired.fastq.gz")  # and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads <- paste0(samples, "filtered_R1_001.paired.fastq.gz")
filtered_reverse_reads <- paste0(samples, "filtered_R2_001.paired.fastq.gz")

filtered_out <- filterAndTrim(fwd = forward_reads,filt =  filtered_forward_reads,rev = reverse_reads,filt.rev =  filtered_reverse_reads, maxEE=c(2,2),
                              rm.phix=TRUE, minLen=175,minQ = 25, truncLen=c(250,200), maxN = 0, multithread = 48, verbose = T)

filtered_out2<-data.frame(filtered_out)
data <- filtered_out2[filtered_out2$reads.out != 0,]
data



filtered_out<-as.matrix(data)

class(filtered_out) # matrix
dim(filtered_out) # 20 2filtered_out
plotQualityProfile(row.names(filtered_out))
head(filtered_out)
#plotQualityProfile(filtered_reverse_reads)
plotQualityProfile(filtered_reverse_reads[1:2])

filtered_forward_reads<-row.names(filtered_out2[filtered_out2$reads.out != 0,])

err_forward_reads <- learnErrors(filtered_forward_reads, multithread=48) 
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=48) 
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples
dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pooled", multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pooled", multithread=TRUE) merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                                                                                                                           derep_reverse, trimOverhang=TRUE, minOverlap=170)# this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
class(merged_amplicons) # list
length(merged_amplicons) # 20 elements in this list, one for each of our samples
names(merged_amplicons) # the names() function gives us the name of each element of the list seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix
dim(seqtab) # 20 2521seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) # Howd many bimeras out of how many input sequences.
# How many sequenced did we lose?
sum(seqtab.nochim)/sum(seqtab) # what is the abundance rate?  0.9931372 # good, we barely lost any in terms of abundance# set a little function
getN <- function(x) sum(getUniques(x))# making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))summary_tabtaxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa.gz", multithread=TRUE, tryRC=T) # giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)




