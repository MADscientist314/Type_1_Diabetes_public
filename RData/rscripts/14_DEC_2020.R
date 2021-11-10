#this script is designed to identify shared taxa in both the mother and the infant

#load libraries and set working direactory
library(phyloseq)
library(microbiome)
library(tidyverse)
setwd("/media/jochum00/Jochum_Raid/T1d")
getwd()

#import the phyloseq object
ASV_physeq_core<-(readRDS("ASV_physeq_core.RDS"))
ASV_physeq_core # [ 2387 taxa and 527 samples ]
#AGGLOMERATE TO GENUS LEVEL
ASV_physeq_core<-tax_glom(ASV_physeq_core,taxrank = "Genus")
ASV_physeq_core #[ 228 taxa and 527 samples ]
#change the names from "ASV_XXX" to the genus name
taxa_names(ASV_physeq_core)<-get_taxa_unique(ASV_physeq_core,taxonomic.rank = "Genus")
#melt the phyloseq object into a tibble
ASV_physeq_core_comp<-microbiome::transform(transform = "compositional",x = ASV_physeq_core)


################################################################################
#### TRANSMISSION ANALYSIS
#identify shared taxa in both the mother and the infant BY PATIENT ID
################################################################################
#convert the physeq obj to a dataframe
# df<-as_tibble(psmelt(ASV_physeq_core))
# df
# df<-df%>%filter(!is.na(disease))
# unique(df$disease)
# #make a sheet that shows the status of the sample
# Deliv<-df%>%select(patient,disease,Host,SampleType,Delivery)
# status<-Deliv%>%select(patient,disease,Delivery)%>%filter(!is.na(Delivery))%>%distinct()%>%arrange(patient)
# status
# #merge the patient and disease columns to a unique idetnifier
# df$subject<-paste0(df$patient,"_",df$disease,"_",df$Host)#,"_",df$SampleType,"_",df$Delivery)
# #filter out emtpy rouws
# df2<-df %>%select(OTU,Abundance,subject)%>%arrange(subject)%>%filter(Abundance>0)
# df2
# #sum up the the Genus Counts and pivot them into a count table
# df3<-df2%>%
#   group_by(OTU,subject)%>%
#   summarise(count=sum(Abundance))%>%
#   arrange(OTU,subject)%>%arrange(subject)
# df3
# 
# ##OK HERE IS THE MONEY MAKER STEP
# # Remove TAXA that are not shared by both subejects (maybe I should of pivoted wider earlier)
# #remove samples that dont have both infant and mother samples remaining
# df4<-df3%>%filter(!duplicated(subject))%>%arrange(subject)
# df%>%pivot_wider(names_from = OTU,values_from=count)
# df3
# 
# ##OK HERE IS THE MONEY MAKER STEP
# # Remove columns that are not shared by both subejects (maybe I should of pivoted wider earlier)
# #df4<-df3%>%distinct(OTU,subject)%>%arrange(subject,Host)
# #df4
# 


################TEST CODE############################
#convert the physeq obj to a dataframe
df<-as_tibble(psmelt(ASV_physeq_core))

#make a sheet that shows the status of the sample
Deliv<-df%>%select(patient,disease,Host,SampleType,Delivery)
status<-Deliv%>%select(patient,disease,Delivery)%>%
  filter(!is.na(Delivery))%>%
  distinct()%>%arrange(patient)
status$subject<-paste0(status$patient,"_",status$disease)


#merge the patient and disease columns to a unique idetnifier
df$subject<-paste0(df$patient,"_",df$disease)#,"_",df$SampleType,"_",df$Delivery)
#filter out emtpy rouws
df2<-df %>%select(OTU,Abundance,subject,Host)%>%
  arrange(subject)%>%filter(Abundance>0)
df2
#sum up the the Genus Counts and pivot them into a count table
df3<-df2%>%
  group_by(OTU,subject,Host)%>%
  summarise(count=sum(Abundance))%>%
  arrange(OTU,subject,Host)%>%arrange(subject)
df3 # A tibble: 14,329 x 4

##OK HERE IS THE MONEY MAKER STEP
# Remove TAXA that are not shared by both subejects (maybe I should of pivoted wider earlier)
#remove samples that dont have both infant and mother samples remaining
df4<-df3%>%filter(!duplicated(subject))%>%
  arrange(subject)
df4 # A tibble: 10,695 x 4

#OMG I think it did it!
#lets pivot wider and find out
df5<-df4%>%arrange(subject,Host)#%>%pivot_wider(names_from = OTU,values_from=count)
status
df5
df6<-full_join(status,df5, by = "subject")%>%
  arrange(subject,Host)%>%
  pivot_wider(names_from = OTU,values_from=count)
df6
df6$subject<-NULL
df6<-df6%>%filter(!is.na(disease))
df6$patient<-as.character(df6$patient)
unique(df6$disease)

# get rid of the lines that arent duplicated (ie: no mother, not infant)
# 
# #all this is doing is getting rid of the patient column and I dont know why!
# 
# df7<-df6#%>%filter(!duplicated(subject))%>%arrange(subject) 
# df6$subject
# df7<-df7%>%filter(duplicated(subject)==TRUE)%>%arrange(subject)
# 
# df7$subject
# #df7<-df7%>%filter(!duplicated(subject))%>%arrange(subject)
# T1D_47<-df7%>%filter(subject=="47_T1D")
# T1D_47
# df7
#ok, lets pivot our taxa back into rows so we can filter out the NA's
df8<-df6%>%
  pivot_longer(cols=c(-patient,-disease,-Delivery,-Host),
               names_to="Genus",values_to="count")%>%
  filter(!is.na(count))
df8$fill<-NULL
unique(df8$disease)
#alright lets plot it out
library(ggpubr)
a<-ggballoonplot(data = df8,x = "patient",y = "Genus",color = "Host",size = "count",facet.by = "disease")
a


#lets try it out with just one patient
T1D_47<-df8%>%filter(patient=="47")%>%filter(disease=="T1D")%>%pivot_wider(names_from = Genus,values_from=count)
T1D_47

#well I am stopping here, but I was thinking maybe we need a vector with 
# taxa by ????? that can be used to filter out taxa that arent shared amounts mothers and infants
# also maybe make a line plot with the mother on the left and the infant on the right and connection
#https://rpkgs.datanovia.com/ggpubr/reference/ggline.html

