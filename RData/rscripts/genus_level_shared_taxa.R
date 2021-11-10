#this script is designed to identify shared taxa in both the mother and the infant
#load libraries and set working direactory
library(phyloseq)
library(microbiome)
library(tidyverse)
setwd("/media/jochum00/Jochum_Raid/T1d")
getwd()


cont<-as_tibble(read.table("picrust2/picrust2_out2/pathways_out/path_abun_contrib.tsv",sep = "\t",header = T,check.names = F))
want<-as_tibble(read.table("short_shared_pathways.txt",sep = "\t",header = T,check.names = F))
colnames(cont)<-c("sample","Function","taxon","taxon_abun","taxon_rel_abun","genome_function_count","taxon_function_abun","taxon_rel_function_abun")
selected<-cont%>%filter(Function %in% want$`BioCyc ID`)
selected

ASV_physeq_core<-readRDS("ASV_physeq_core.RDS")
tax<-as.matrix(tax_table(ASV_physeq_core))
tax<-as.data.frame(tax)
tax<-as_tibble(tax,rownames = "taxon")
sam<-as.data.frame(sample_data(ASV_physeq_core))
sam<-as_tibble(sam,rownames="sample")
sam<-sam%>%select(sample,patient,disease,Host,SampleType,Delivery)



colnames(want)<-c("Function","Shared_pathway")
want

long<-selected%>%pivot_longer(cols = c(taxon_abun,
                                       taxon_abun,
                                       taxon_function_abun,
                                       taxon_rel_abun,
                                       taxon_rel_function_abun,
                                       genome_function_count),
                              names_to = "abund_factor", 
                              values_to="count")
long_contributions<-full_join(long,tax)
long_contributions
sam
final<-full_join(long_contributions,sam)%>%filter(!is.na(Function))
final
want
final<-left_join(want,final)
final<-final%>%pivot_wider(names_from = abund_factor,values_from=count)
write.table(final,"taxonomic_contributions_to_predicted_functions.tsv",sep = "\t",quote = F,row.names = F)





unique(final$abund_factor)
[1] "taxon_abun"              "taxon_function_abun"     "taxon_rel_abun"          "taxon_rel_function_abun"
[5] "genome_function_count"   NA                       

####################################################################
taxon_abun<-final%>%filter(abund_factor=="taxon_abun")%>%pivot_wider(names_from = abund_factor,values_from=count)
write.table(taxon_abun,"tax_abund_contrib_to_predicted_functions.tsv",sep = "\t",quote = F,row.names = F)
####################################################################
taxon_function_abun<-final%>%filter(abund_factor=="taxon_function_abun")%>%pivot_wider(names_from = abund_factor,values_from=count)
write.table(taxon_function_abun,"tax_function_abund_contrib_to_predicted_functions.tsv",sep = "\t",quote = F,row.names = F)

taxon_abun



as.data.frame(final)
write.table(final,"taxonomic_contributions_to_predicted_functions.tsv",sep = "\t",quote = F,row.names = F)

getwd()
