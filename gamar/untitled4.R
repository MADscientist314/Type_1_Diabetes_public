library(biom)
library(biomformat)
library(phyloseq)
library(metagenomeSeq)
setwd("/media/jochum00/Jochum_Raid/T1d/sourcetracker/sourcetracker")

setwd("/media/jochum00/Jochum_Raid/T1d/sourcetracker/sourcetracker")

ASV_physeq_core2<-prune_samples(x = ASV_physeq_core,sample_sums(ASV_physeq_core)>0)
library(microbiome)
ASV_physeq_core2<-transform(ASV_physeq_core,"compositional")
s<-data.frame(sample_data(ASV_physeq_core2))%>%
  mutate(SampleType=gsub(" Introitus","Introitus",SampleType))
s

t<-s%>%filter(Host=="Mother")%>%mutate(Host=paste0(SampleType))
t<-t%>%mutate(SampleID=rownames(t))
t
t
u<-s%>%filter(Host=="Infant")
u<-u%>%mutate(SampleID=rownames(u))
u
v<-full_join(t,u)
v

v$patient_disease
v<-v %>%select(SampleID,SampleType,patient_disease,Host)%>%distinct_all()%>%
  mutate(SourceSink=gsub("Infant","sink",Host),
         SourceSink=gsub("Introitus","source",SourceSink),
         SourceSink=gsub("Cervix","source",SourceSink),
         SourceSink=gsub("Anus","source",SourceSink),
         SourceSink=gsub("Vagina","source",SourceSink))

length(v$SampleID)
length(sample_names(ASV_physeq_core2))
rownames(v)<-v$SampleID
v
sample_names(ASV_physeq_core2)
sample_data(ASV_physeq_core2)<-sample_data(v)
unique(v$Host)

sam<-data.frame(v)#sample_data(ASV_physeq_core_genus))

getwd()
sam<-sam%>%select(index,Stample)
write.table(sam,"sam.tsv",sep = "\t",row.names = F,quote = F)

sam


#convert the phyloseq object to a metagenomeseq object
a<-phyloseq_to_metagenomeSeq(ASV_physeq_core2)
#conver the metagenomeseq object to a biom format
b<-metagenomeSeq::MRexperiment2biom(a)
#write the biom object to a file
write_biom(x = b,biom_file =paste0("otu.biom"))

sample_names(ASV_physeq_core2)
