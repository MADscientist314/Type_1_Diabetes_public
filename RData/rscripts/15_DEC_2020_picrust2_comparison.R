#this script is designed to identify shared taxa in both the mother and the infant
#load libraries and set working direactory
library(phyloseq)
library(microbiome)
library(tidyverse)
setwd("/media/jochum00/Jochum_Raid/T1d")
getwd()


#############preprocessing####################
ASV_physeq_core<-(readRDS("ASV_physeq_core.RDS"))
ASV_physeq_core
#agglomerate to a genus level
ASV_physeq_core_genus<-tax_glom(ASV_physeq_core, taxrank = "Genus")
#Get only Bacteria
ASV_physeq_core_genus<-subset_taxa(Kingdom=="Bacteria")
# do a prevelance and detection cutoff at 1%
ASV_physeq_core_genus<-core(x = ASV_physeq_core_genus,detection = 1/100,prevalence = 1/100)
# ###Export your phyloseq into tsv files for biom conversion
# otu<-write.table(x = data.frame(otu_table(ASV_physeq_core_genus)),file = "ASV_physeq_core_genus_counts.tsv")
# tax<-write.table(x = data.frame(tax_table(ASV_physeq_core_genus)),file = "ASV_physeq_core_genus_taxonomy.tsv")
# sam<-write.table(x = data.frame(sampled_data(ASV_physeq_core_genus)),file = "ASV_physeq_core_genus_metadata.tsv")


otu<-read.table("picrust2/ASV_physeq_core_genus_counts.tsv",sep = "\t",header = T,row.names = 1, check.names = F)
tax<-read.table("picrust2/ASV_physeq_core_genus_taxonomy.tsv",sep = "\t",header = T,row.names = 1, check.names = F)
sam<-read.table("picrust2/ASV_physeq_core_genus_metadata.tsv",sep = "\t",header = T,row.names = 1, check.names = F)
otu

paths<-read.table("picrust2/picrust2_out_pipeline/pathways_out/path_abun_unstrat.tsv",sep = "\t", header = T,row.names = 1,check.names = F)
paths<-otu_table(paths,taxa_are_rows = T)
paths_tax<-read.table("picrust2/picrust2_out_pipeline/pathways_out/path_abun_unstrat_tax.txt",header = T, sep = "\t", row.names =1, quote = "", comment.char = "")
paths_tax<-tax_table(as.matrix(paths_tax),errorIfNULL = T)
sam<-read.table("picrust2/ASV_physeq_core_genus_metadata.tsv",sep = "\t",header = T,row.names = 1, check.names = F)
sam<-sample_data(sam)


path_physeq<-phyloseq(paths,paths_tax,sam)
path_physeq_core<-core(x =path_physeq,detection = 1/100,prevalence = 1/100)
path_physeq_core

#change the names from "path_XXX" to the genus name
path_physeq_core<-tax_glom(path_physeq_core,taxrank = "description")
path_physeq_core
taxa_names(path_physeq_core)<-get_taxa_unique(physeq = path_physeq_core,"description")

#melt the phyloseq object into a tibble
#get the relabunds
path_physeq_core_comp<-microbiome::transform(transform = "compositional",x = path_physeq_core)

length(unique(df$subject))
################TEST CODE############################
#convert the physeq obj to a dataframe
df<-as_tibble(psmelt(path_physeq_core))
df# A tibble: 173,383 x 163





#make a sheet that shows the status of the sample
Deliv<-df%>%select(patient,disease,Host,SampleType,Delivery)
status<-Deliv%>%select(patient,disease,Delivery)%>%filter(!is.na(Delivery))%>%distinct()%>%arrange(patient)
status$subject<-paste0(status$patient,"_",status$disease)

df
#merge the patient and disease columns to a unique idetnifier
df$subject<-paste0(df$patient,"_",df$disease)#,"_",df$SampleType,"_",df$Delivery)
#filter out emtpy rouws
df2<-df %>%select(OTU,Abundance,subject,Host,disease)%>%
  arrange(subject)%>%
  filter(Abundance>0)%>%
  arrange(OTU)
df2
#Get a tibble of moths and infants separeately
#then sum up the counts that came from each sampletype 
#so we just get 1 otu per row per sample
mother<-df %>%
  select(OTU,Abundance,subject,Host)%>%
  arrange(subject)%>%filter(Abundance>0)%>%
  filter(Host=="Mother")%>%
  group_by(OTU,subject)%>%
  summarise(Abundance=sum(Abundance))

mother_OTU<-mother %>% count(OTU)
mother_OTU


infant<-df %>%
  select(OTU,Abundance,subject,Host)%>%
  arrange(subject)%>%
  filter(Abundance>0)%>%
  filter(Host=="Infant")%>%
  group_by(OTU,subject)%>%
  summarise(Abundance=sum(Abundance))



#join the matching cells from each list to get a list of samples that 
#have both mothers and infants
df3<-inner_join(mother, infant, by = c("OTU", "subject"),
                suffix=c("_Mother","_Infant"))%>%
  select(OTU,subject,Abundance_Mother,Abundance_Infant)

#extra tibble witht the total count just in case
df3_optional<-df3%>%
  mutate(Abundance_total=Abundance_Mother+Abundance_Infant)



####YAY YOU DID IT!!!!
df3
status
dim(df2)
dim(df3)
min(df3$Abundance_Mother)
min(df3$Abundance_Infant)
tail(df3)

#double check the transformation by getting a list of the unique lists and comparing lengths

mother_list<-unique(mother$subject)
infant_list<-unique(infant$subject)
tot_list<-unique(df3$subject)
tot_list
length(tot_list)
length(mother_list)
length(infant_list)
df3

#do some data wrangling 
#combined the abundances into a single column
df4<-df3%>%
  pivot_longer(cols = c(Abundance_Mother,Abundance_Infant),
               names_to="Host",values_to="Abundance")

#lets fix the names
df4$Host<-gsub(x = df4$Host,pattern = "Abundance_",replacement = "")
df4
#ok lets data wrangle the Delivery back in
df5<-full_join(status,df4, by = "subject")%>%
  arrange(subject,Host)

###WHY DOES PATTIENT 42 CONTROL VAGINA NOT HAVE A HOST OR ABUNDANCE?
#PROBABLY BC IT WAS REMOVED FROM THE SAMPLE SET
df6<-df5%>%filter(is.na(Host))
# A tibble: 1 x 7
#patient disease Delivery subject    OTU   Host  Abundance
#<int> <chr>   <chr>    <chr>      <chr> <chr>     <dbl>
#  1      42 Control Vaginal  42_Control NA    NA           NA

df5#47,253 x 7
#OK WELL LETS get rid of the NA Host samples
df5<-df5%>%filter(!is.na(Host))
df5# 47,252 x 7

library(mosaic)

#set this function, idk why
shlab <- function(lbl_brk){
  sub("^[a-z]+\\.","",lbl_brk) # removes the starts of strings as a. or ab.
}


summary<-as_tibble(tally(~ OTU+Delivery+disease+Host, 
                         data =df5, 
                         format="data.frame"))%>%
  filter(Freq>0)%>%
  arrange(desc(Freq))

summary2<-as_tibble(tally(~ OTU+Delivery+disease+Host, 
                         data =df5, 
                         format="percent"))%>%
  filter(n>0)%>%
  arrange(desc(n))

a<-31/92
a
summary#1,309 x 5
summary2
tail(summary)
tail(summary2)
sum2<-summary%>%pivot_wider(names_from = Host,values_from=Freq)
sum2#1,309 x 5
sum2<-summary%>%pivot_wider(names_from = Host,values_from=Freq)%>%filter(Infant>0)%>%filter(Mother>0)
tail(sum2)
#alright lets plot it out
#make a like of the order OTU names so ggpubr doesnt jack it up alphabetically

ggballoonplot(Stack_Overflow_DummyData, x = "Region_prefix", y = "Species_prefix", size = "Frequency", size.range = c(1, 9), fill = "Dist") +
  scale_x_discrete(labels = shlab) +
  scale_y_discrete(labels = shlab) +
  scale_fill_manual(values = c("green", "blue", "red", "black", "white")) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  theme_set(theme_gray() + theme(legend.key=element_blank())) +     # Sets Grey Theme and removes grey background from legend panel
  theme(axis.ticks.y = element_text(size = 5))+
# Add Frequency Values Next to the circles

library(ggpubr)
a<-ggballoonplot(data = summary,x = "Host",
                 y = "OTU",
                 color = "disease",
                 size = "Freq",
                 size.range = c(1, 2))+
  scale_y_discrete(labels = shlab) +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  # Sets Grey Theme and removes grey background from legend panel
  theme(axis.text.y = element_text(size = 5))+facet_grid(~Delivery)
  
a

df5
df6
######ok well we need to scale this down some more
summary<-as_tibble(tally(OTU~Delivery+disease,data =df5, 
                           format="data.frame"))%>%
  filter(Freq>0)%>%
  arrange(desc(Freq))

summary$Freq<-summary$Freq/2 #divides by 2 bc infant and mom

summary

summary2<-as_tibble(tally(OTU~Delivery+disease,data =df5, 
                         format="percent"))%>%
  filter(n>0)%>%
  arrange(desc(n))
summary#1,309 x 5
tail(summary)

summary2$n<-summary2$n*100
head(summary2)
tail(summary2)



colnames(summary)<-c("OTU","Delivery","disease","num_mother_infant_samples")
colnames(summary2)<-c("OTU","Delivery","disease","perc_samples")

summary
tail(summary)
gghistogram(summary,
            x = "num_mother_infant_samples",
            y = "..count..",
            fill = "disease",
            facet.by = "Delivery",
            title = "Histogram showing number of mother infant shared",
            alpha = 0.1) 
gghistogram(summary2,
            x = "perc_samples",
            y = "..count..",
            fill = "Delivery",
            facet.by = "OTU",
            title = "Histogram showing number of mother infant shared",
            alpha = 0.1) 

#alright lets plot it out
  #make a like of the order OTU names so ggpubr doesnt jack it up alphabetically
  
  ggballoonplot(Stack_Overflow_DummyData, x = "Region_prefix", y = "Species_prefix", size = "Frequency", size.range = c(1, 9), fill = "Dist") +
    scale_x_discrete(labels = shlab) +
    scale_y_discrete(labels = shlab) +
    scale_fill_manual(values = c("green", "blue", "red", "black", "white")) +
    guides(fill = guide_legend(override.aes = list(size=8))) +
    theme_set(theme_gray() + theme(legend.key=element_blank())) +     # Sets Grey Theme and removes grey background from legend panel
    theme(axis.ticks.y = element_text(size = 5))+
    # Add Frequency Values Next to the circles

hist(summary2$perc_samples)
summary3<-summary2%>%filter(perc_samples>30)
library(ggpubr)
  a<-ggballoonplot(data = summary3,x = "Delivery",
                   y = "OTU",
                   color = "disease",
                   size = "perc_samples",
                   size.range = c(1, 2))+
    scale_y_discrete(labels = shlab) +
    guides(fill = guide_legend(override.aes = list(size=8))) +
    # Sets Grey Theme and removes grey background from legend panel
    theme(axis.text.y = element_text(size = 5))#+facet_grid(~Delivery)

df2
sam2<-as_tibble(sam)
sam2<-sam2%>%filter(!is.na(Delivery))
nsamples<-as_tibble(tally(~ Delivery+disease+Host, 
                         data =sam2, 
                         format="data.frame"))
nsamples

my_tbl <- tally(disease~OTU,
                data = df5)
mosaicplot(my_tbl, color = TRUE)
my_tbl
df5
prop.test(disease~OTU,
          success = "T1D",
          data = df5)



mod <- lm(Abundance ~ OTU*disease*Delivery,
          data = df5)
anova(mod)
df6<-df5
df6$OTU<-as.factor(df6$OTU)
mod <- lm(Abundance*OTU~disease*Delivery,
          data = df6)
anova(mod)
#Which differences are significant?
#mplot(TukeyHSD(mod))
#  Logistic regression
df6$OTU<-as.factor(df6$OTU)
df6$Csection<-df6$Delivery=="Csection"
df6$Vaginal<-df6$Delivery=="Vaginal"
df6$T1D<-df6$disease=="T1D"
df6$Control<-df6$disease=="Control"

summary<-as_tibble(tally(OTU~Delivery+disease,data =df5, 
                         format="data.frame"))%>%
  filter(Freq>0)%>%
  arrange(desc(Freq))

summary$Freq<-summary$Freq/2 #divides by 2 bc infant and mom

#install.packages("ggVennDiagram")
library(ggVennDiagram)

genes
set.seed(20190708)
x <- list(A=sample(genes,300),B=sample(genes,525),C=sample(genes,440),D=sample(genes,350))
x<-list(summary)
x
# four dimension venn plot
ggVennDiagram(x)

x<-tally()
x


logit_mod <- glm(T1D*Csection~OTU*Abundance,
                   data = df6,
                   family = binomial)
msummary(logit_mod)
#Odds ratios and confidence intervals
library(broom)
logit_mod_conf<-tidy(logit_mod, conf.int = TRUE,
       exponentiate = TRUE)


