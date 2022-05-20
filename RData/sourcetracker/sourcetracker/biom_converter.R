install.packages("devtools") # if not already installed
library(devtools)
devtools::install_github("https://github.com/joey711/biom/")

BiocManager::install("biom")
library(biom)
library(biomformat)
library(phyloseq)
library(metagenomeSeq)

data = biom_data(ASV_physeq_core_genus)
smd = sample_metadata(phylo_object)
omd = observation_metadata(phylo_object)

biom_object = make_biom(data, smd, omd)

metagenomeSeq_object = biom2MRexperiment(biom_object)
a<-phyloseq_to_metagenomeSeq(ASV_physeq_core_genus)
b<-metagenomeSeq::MRexperiment2biom(a)
b
write_biom(x = b,biom_file = "ASV_physeq_core_genus.biom")

sam<-data.frame(sample_data(ASV_physeq_core_genus))%>%
  select(c(index,SourceSink,patient_disease,Host))
colnames(sam)<-c("SampleID","SourceSink","patient_disease","Host")
sam

write.table(x = sam,file = "sam.tsv",sep = "\t",row.names = F,quote = F)
df<-read.table("test/mixing_proportions.txt",header = T,sep = "\t",row.names = 1,check.names = F)
df<-as_tibble(t(df),rownames = "sampleID")
df
sam<-data.frame(sample_data(ASV_physeq_core))
sam$sampleID<-rownames(sam)
df2<-inner_join(df,sam)
df2

df3<-df2%>%select(disease,SampleType,Anus,Cervix,Introitus,Vagina,Unknown)#%>%pivot_longer(cols = -disease,names_to = "Source",values_to = "source_val" )

df3<-df3%>%pivot_longer(cols = -c(disease,SampleType),names_to = "Source",values_to = "source_val" )
# df3
# df3<-df2%>%select(disease,Anus,Cervix,Introitus,Vagina,Unknown)%>%group_by(disease)%>%summarise_all(.funs = mean)
# 
# df3<-df3%>%pivot_longer(cols = -disease,names_to = "SampleType",values_to = "source" )
library(ggpubr)
my_comparisons<-list(c("T1D","Control"))
ggboxplot(data = df3,x = "disease",y = "source_val",color = "disease",palette = "npg",add = "jitter",ggtheme = theme_bw())+
  yscale("log10",T)+
  facet_grid(SampleType~Source,scales = "free")+
  stat_compare_means(comparisons = my_comparisons)
#ok now that that is finished lets re-run this on a pt specfic loop

#############################################################################
#############################################################################
#############################################################################
ps<-psmelt(ASV_physeq_core_genus)
pt<-NULL
loop<-function(pt){
  print(paste0("LOG: running ", pt))
  ps<-ps%>%filter(patient_disease%in%pt)%>%filter(Abundance>1)
  ps
  otu<-ps%>%select(OTU,Abundance,Sample)%>%pivot_wider(names_from = Sample,
                                                       values_from = Abundance,
                                                       #values_fn = mean,
                                                       values_fill = 0)
  otu
  physeq_otu<-data.frame(otu%>%select(-OTU),check.names = F)
  physeq_otu
  rownames(physeq_otu)<-otu$OTU
  physeq_otu<-otu_table(physeq_otu,taxa_are_rows = T)
  physeq_otu
  tax<-ps%>%select(OTU,Kingdom,Phylum,Class,Order,Family,Genus)%>%distinct_all()
  physeq_tax<-as.matrix(tax%>%select(-OTU))
  rownames(physeq_tax)<-tax$OTU
  physeq_tax<-tax_table(physeq_tax)
  physeq_tax
  physeq_sam<-data.frame(ps)%>%select(Sample,patient_disease,index,Host,SourceSink)%>%distinct_all()
  rownames(physeq_sam)<-physeq_sam$Sample
  physeq_sam<-sample_data(physeq_sam)
  physeq_sam
  physeq_otu
  #put it together
  physeq<-phyloseq(physeq_otu,physeq_tax,physeq_sam)
  physeq
    
  a<-phyloseq_to_metagenomeSeq(physeq)
  b<-metagenomeSeq::MRexperiment2biom(a)
  
  write_biom(x = b,biom_file =paste0("patient_disease/",pt,"/otu.biom"))
  
  sam<-data.frame(sample_data(ASV_physeq_core_genus))%>%
    select(c(index,SourceSink,patient_disease,Host))
  colnames(sam)<-c("SampleID","SourceSink","patient_disease","Host")
  sam
  
  write.table(x = sam,file = paste0("./patient_disease/",pt,"/sam.tsv"),sep = "\t",row.names = F,quote = F)
  # ps<-ps%>%filter(patient_disease%in%pt)%>%filter(Abundance>1)
  # ps<-ps%>%mutate(SampleType=gsub(" Introitus","Introitus",SampleType))
  # ps<-as_tibble(ps)%>%mutate(Sample=SampleType)
  # 
  # otu<-ps%>%select(OTU,Abundance,Sample)%>%pivot_wider(names_from = Sample,
  #                                                      values_from = Abundance,
  #                                                      values_fn = mean,
  #                                                      values_fill = 0)
  # otu
  # 
  # sam<-ps %>%select(Sample,Host)%>%distinct_all()%>%
  #   mutate(SourceSink=gsub("Infant","sink",Host),
  #          SourceSink=gsub("Mother","source",SourceSink),
  #          SourceSink=gsub("Introitus","source",SourceSink),
  #          SourceSink=gsub("Cervix","source",SourceSink),
  #          SourceSink=gsub("Anus","source",SourceSink),
  #          SourceSink=gsub("Vagina","source",SourceSink))
  # #sam
  # print(paste0("LOG: printing otu ", pt))
  # write.table(otu,paste0("patient_disease/",pt,"/otu.tsv"),sep = "\t",row.names = F,quote = F)
  # print(paste0("LOG: printing sam.tsv ", pt))
  # write.table(sam,paste0("./patient_disease/",pt,"/sam.tsv"),sep = "\t",row.names = F,quote = F)
}

sapply(l$patient_disease,function(X)loop(pt = X))
