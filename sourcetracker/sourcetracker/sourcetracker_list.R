tax_table(vst_physeq)<-tax_table(ASV_physeq_core)

library(microbiome)
library(phyloseq)
library(tidyverse)
getwd()
setwd("/media/jochum00/Jochum_Raid/T1d/sourcetracker/sourcetracker")
ASV_physeq_core_genus<-tax_glom(ASV_physeq_core,taxrank="Genus",NArm=T)
ASV_physeq_core_genus<-prune_samples(x = ASV_physeq_core_genus,sample_sums(ASV_physeq_core_genus)>0)

s<-data.frame(sample_data(ASV_physeq_core_genus))
s$Host
t<-s%>%filter(Host=="Mother")%>%mutate(Host=paste0(SampleType))%>%mutate(index=rownames(t))
t$Host
t
u<-s%>%filter(Host=="Infant")%>%mutate(Host=paste0(SampleType))%>%mutate(index=rownames(u))
u
v<-full_join(t,u)

rownames(v)<-v$index
sample_data(ASV_physeq_core_genus)<-sample_data(v)


l<-s%>%select(patient_disease)%>%distinct_all()
l
#write.table(l,"sample_manifest.tsv",row.names = F,quote = F)
#lets subset just one of the dyads for now
ps$Sample
ps<-psmelt(ASV_physeq_core_genus)
loop<-function(pt){
  print(paste0("LOG: running ", pt))
  ps<-ps%>%filter(patient_disease%in%pt)%>%filter(Abundance>1)
  ps<-ps%>%mutate(SampleType=gsub(" Introitus","Introitus",SampleType))
  ps<-as_tibble(ps)%>%mutate(Sample=SampleType)
  
  otu<-ps%>%select(OTU,Abundance,Sample)%>%pivot_wider(names_from = Sample,values_from = Abundance,values_fn = mean, values_fill = 0)
  otu
  
  sam<-ps %>%select(Sample,Host)%>%distinct_all()%>%
    mutate(SourceSink=gsub("Infant","sink",Host),
           SourceSink=gsub("Mother","source",SourceSink))
  #sam
  print(paste0("LOG: printing otu ", pt))
  write.table(otu,paste0("patient_disease/",pt,"/otu.tsv"),sep = "\t",row.names = F,quote = F)
  print(paste0("LOG: printing sam.tsv ", pt))
  write.table(sam,paste0("./patient_disease/",pt,"/sam.tsv"),sep = "\t",row.names = F,quote = F)
}

#shell script
#(base) jochum00@aagaardlab3:/media/jochum00/Jochum_Raid/T1d/sourcetracker/sourcetracker$ docker run -it -v /media/jochum00/:/media/jochum00 -w $PWD --name sourcetracker2 yangkl0909/sourcetracker2:0.1 bash
#for n in $(cat sample_manifest.tsv); do sourcetracker2 gibbs -i ./patient_disease/$n/otu.tsv -m ./patient_disease/$n/sam.tsv -o patient_disease/$n/ --jobs 48 --alpha1 0.001 --alpha2 0.001 --beta 10 --source_rarefaction_depth 1000 --sink_rarefaction_depth 1000 --restarts 10 --draws_per_restart 1 --burnin 100 --delay 10 --cluster_start_delay 10 --source_sink_column SourceSink --source_column_value source --sink_column_value sink --source_category_column Host; done


ps<-ps%>%filter(Abundance>1)
ps<-ps%>%mutate(SampleType=gsub(" Introitus","Introitus",SampleType))
ps<-as_tibble(ps)%>%mutate(Sample=SampleType)

otu<-ps%>%select(OTU,Abundance,Sample)%>%pivot_wider(names_from = Sample,values_from = Abundance,values_fn = mean, values_fill = 0)
otu

sam<-ps %>%select(Sample,Host)%>%distinct_all()%>%
  mutate(SourceSink=gsub("Infant","sink",Host),
         SourceSink=gsub("Mother","source",SourceSink))
#sam
write.table(otu,paste0("patient_disease/otu.tsv"),sep = "\t",row.names = F,quote = F)
write.table(sam,paste0("patient_disease/sam.tsv"),sep = "\t",row.names = F,quote = F)


sam$SourceSink


sapply(l$patient_disease,function(X)loop(pt = X))

tmp<-subset_samples(physeq = ASV_physeq_core_genus,patient_disease=="20_T1D")
tmp
tmp<-prune_samples(sample_sums(tmp)>=1, tmp)
tmp
tmp<-prune_taxa(taxa_sums(tmp) > 0, tmp) 
tmp# [ 134 taxa and 6 samples ]
ps<-psmelt(tmp)
#ps_backup<-ps
ps<-ps%>%mutate(SampleType=gsub(" Introitus","Introitus",SampleType))
ps<-as_tibble(ps)%>%mutate(Sample=SampleType)



otu<-ps%>%select(OTU,Abundance,Sample)%>%pivot_wider(names_from = Sample,values_from = Abundance)
otu
sam<-ps %>%select(Sample,Host)%>%distinct_all()%>%
  mutate(SourceSink=gsub("Infant","sink",Host),
         SourceSink=gsub("Mother","source",SourceSink))
sam

write.table(otu,"otu.tsv",sep = "\t",row.names = F,quote = F)
write.table(sam,"sam.tsv",sep = "\t",row.names = F,quote = F)


sample






otu
sam









ps
#convert the samples

dat<-ps%>%select(OTU,id,Host,SampleType,Abundance)%>%distinct_all()
dat
dat<-dat%>%pivot_wider(names_from = OTU,values_from = Abundance, values_fill = 0)
dat
#make a maktris
mat_dat<-data.frame(dat)%>%select(-c(id,Host,SampleType))
rownames(mat_dat)<-dat$id
mat_dat


tax<-ps%>%select(c(OTU,Kingdom,Phylum,Class,Order,Family,Genus))%>%distinct_all()
tax

tax
sam<-ps%>%select(-c(OTU,Abundance,Kingdom,Phylum,Class,Order,Family,Genus))%>%distinct_all()
sam<-sam%>%select(id,patient,disease,patient_disease,Host,SampleType,disease_SampleType,disease_Host_SampleType,Delivery)
sam<-data.frame(sam)
rownames(sam)<-sam$id
sam

#call the infant and mother source and sink

sam<-sam%>%mutate(SourceSink=gsub("Infant","sink",Host),
                  SourceSink=gsub("Mother","source",SourceSink))
sam$SourceSink
# load sample metadata
common.sample.ids <- intersect(rownames(sam), rownames(mat_dat))
otus<- mat_dat[common.sample.ids,]#%>%select(-c(Host,SampleType))
otus
metadata <- sam[common.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

otus


#otus<-ceiling(otus)
# extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='source')
otus[train.ix,]
metadata$SourceSink
test.ix<- which(metadata$SourceSink=='sink')
test.ix
envs<-metadata$SampleType
envs
if(is.element('SampleType',colnames(metadata))) 
desc<- metadata$SampleType
desc
class(otus[train.ix,])
class(envs[train.ix])

# load SourceTracker package
source('src/SourceTracker.r')
ps


# tune the alpha values using cross-validation (this is slow!)
#tune.results <- tune.st(otus = otus[train.ix,],envs = envs[train.ix],rarefaction_depth = 1000,verbosity = T,individual.samples = T,ntrials = 25)
# alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001
data.frame(otus[train.ix,])
# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix],rarefaction_depth = NULL)
st
ts <- sourcetracker(otus[test.ix,], envs[test.ix],rarefaction_depth = NULL)

# Estimate source proportions in test data
results<-predict.sourcetracker(stobj = st,test = ts$train)#otus[test.ix,])#, alpha1=alpha1, alpha2=alpha2,verbosity = T,rarefaction_depth = NULL, full.results = T)
results
# Estimate leave-one-out source proportions in training data 
results.train <- predict(st, alpha1=alpha1, alpha2=alpha2,verbosity = T,full.results =T,rarefaction_depth = NULL)

# plot results
labels <- sprintf('%s %s', envs,desc)
plot(results, labels[test.ix], type='pie')

# other plotting functions
plot(results, labels[test.ix], type='bar')
plot(results, labels[test.ix], type='dist')
plot(results.train, labels[train.ix], type='pie')
plot(results.train, labels[train.ix], type='bar')
plot(results.train, labels[train.ix], type='dist')

# plot results with legend
plot(results, labels[test.ix], type='pie', include.legend=TRUE, env.colors=c('#47697E','#5B7444','#CC6666','#79BEDB','#885588'))







list<-unique(sam$patient_disease)
list
#make a function that iteratively does this per patient
subme<-function(pt){
  print (paste0("LOG: pruning ",pt))
  res<-prune_samples(x = ASV_physeq_core_genus,samples= sample_data(ASV_physeq_core_genus)$patient_disease %in% pt)
  res<-prune_taxa(taxa_sums(res) > 0, res)
  res<-prune_samples(sample_sums(res) > 0,res)
  res
  print (paste0("LOG: psmelting ",pt))
  ps<-psmelt(res)
  
  print (paste0("LOG: data wranglinge ",pt))
  ps<-ps%>%mutate(id=paste0("s_",Sample))%>%select(-Sample)%>%arrange(OTU)%>%distinct_all()
  dat<-ps%>%select(OTU,id,Abundance)%>%distinct_all()
  dat<-dat%>%pivot_wider(names_from = OTU,values_from = Abundance, values_fill = 0)
  mat_dat<-data.frame(dat)%>%select(-id)
  rownames(mat_dat)<-dat$id
  mat_dat
  tax<-ps%>%select(c(OTU,Kingdom,Phylum,Class,Order,Family,Genus))%>%distinct_all()
  sam<-ps%>%select(-c(OTU,Abundance,Kingdom,Phylum,Class,Order,Family,Genus))%>%distinct_all()
  sam<-sam%>%select(id,patient,disease,patient_disease,Host,SampleType,disease_SampleType,disease_Host_SampleType,Delivery)
  rownames(sam)<-sam$id
  #call the infant and mother source and sink
  sam<-sam%>%mutate(SourceSink=gsub("Infant","sink",Host),
                    SourceSink=gsub("Mother","source",SourceSink))
  sam$SourceSink
  # load sample metadata
  common.sample.ids <- intersect(rownames(sam), rownames(mat_dat))
  otus<- mat_dat[common.sample.ids,]
  metadata <- sam[common.sample.ids,]
  # double-check that the mapping file and otu table
  # had overlapping samples
  if(length(common.sample.ids) <= 1) {
    message <- paste(sprintf('Error: there are %d sample ids in common '),
                     'between the metadata file and data table')
    stop(message)
  }
  #otus<-ceiling(otus)
  print (paste0("LOG: extracting source sinks from ",pt))
  
  # extract the source environments and source/sink indices
  train.ix <- which(metadata$SourceSink=='source')
  train.ix
  # metadata$SourceSink
  test.ix<- which(metadata$SourceSink=='sink')
  # test.ix
  # envs<-metadata$SampleType
  # envs
  if(is.element('SampleType',colnames(metadata))) 
  desc<- metadata$SampleType
  desc
  #otus[train.ix,]
  #envs[train.ix]
  
  print (paste0("LOG: sourcetrackering ",pt))
  
  # # load SourceTracker package
  source('src/SourceTracker.r')
  # envs[train.ix]
  # # tune the alpha values using cross-validation (this is slow!)
  #tune.results <- tune.st(otus = otus[train.ix,],envs = envs[train.ix],rarefaction_depth = 0,verbosity = T)
  #alpha1 <- tune.results$best.alpha1
  #alpha2 <- tune.results$best.alpha2
  # # note: to skip tuning, run this instead:
  alpha1 <- alpha2 <- 0.001
  # 
  # # train SourceTracker object on training data
  st <- sourcetracker(otus[train.ix,], envs[train.ix],rarefaction_depth = NULL)
  st
  # # Estimate source proportions in test data
  results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2,verbosity = T,rarefaction_depth = NULL, full.results = T)
  # 
  # # Estimate leave-one-out source proportions in training data 
  results.train <- predict(st, alpha1=alpha1, alpha2=alpha2,verbosity = T,full.results =T,rarefaction_depth = NULL)
  # 
  # # plot results
}
lapply(list, subme)
