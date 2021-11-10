getwd()
library(dplyr)
library(tidyverse)
library(mosaic)
sam<-read.table(file = "sample_core_Sept16_2020.txt", header = T, sep = "\t", row.names = 1)
sam
sam2<-tibble(sam)%>%filter(Host=="Infant")
sam2<-tibble(sam2)%>%group_by(patient,disease,newborn_weight,SGA_AGA_LGA_2)%>%distinct(patient)
sam3<-filter(sam2, !is.na(SGA_AGA_LGA_2))
#############################################################
##################################Delivery###################
sam2<-tibble(sam)%>%group_by(patient,disease)%>%distinct(Delivery)

sam2<-filter(sam2, !is.na(Delivery))
Delivery_perc<-tally(Delivery~disease,data =sam2,format = "percent")
Delivery.table<-tally(Delivery~disease, data=sam2)
class(Delivery.table)
Delivery_perc
chisq.test(Delivery_perc)

Deliv<-read.table(file = "Deliv.txt", header = T, sep = "\t", row.names=1)
Deliv
chisq.test(x = Deliv$Cases, y = Deliv$Control)
msummary(chisq(x = Delivery.table))
## Hollander & Wolfe (1973), 116.
## Mucociliary efficiency from the rate of removal of dust in normal
##  subjects, subjects with obstructive airway disease, and subjects
##  with asbestosis.
Caesarian<-c(33,20)
Vaginal<-c(12,17)
Vacuum<-c(5,5)
# with asbestosis
kruskal.test(list(Caesarian,Vaginal,Vacuum))
chisq.test(Caesarian)
chisq.test(Vacuum)

## Formula interface.
require(graphics)
boxplot(Ozone ~ Month, data = airquality)
kruskal.test(Ozone ~ Month, data = airquality)
chisq(Deliv)


#############################################################
#######################Antibiotics_2########################
sam2<-tibble(sam)%>%group_by(patient,disease)%>%distinct(Antibiotics_2)
sam2<-filter(sam2, !is.na(Antibiotics_2))
Antibiotics_2_perc<-tally(Antibiotics_2~disease,data =sam2,format = "percent")
Antibiotics_2_table <- tally(Antibiotics_2~disease, data=sam2)
Antibiotics_2_table
Antibiotics_2_perc
chisq.test(Antibiotics_2_perc)
#############################################################
library(mosaic)
library(stats)
library(tidyverse)
sam2<-tibble(sam)%>%filter(Host=="Infant")
sam2<-tibble(sam2)%>%group_by(patient,disease,newborn_weight,SGA_AGA_LGA_2)%>%distinct(patient)
sam3<-filter(sam2, !is.na(SGA_AGA_LGA_2))
sam3$SGA_AGA_LGA_2<- factor(sam3$SGA_AGA_LGA_2,levels = c("SGA","AGA","LGA"))
SGA_AGA_LGA_2_perc<-tally(SGA_AGA_LGA_2~disease,data =sam3,format = "percent")
SGA_AGA_LGA_2_table <- tally(SGA_AGA_LGA_2~disease, data=sam3)
SGA_AGA_LGA_2_table
SGA_AGA_LGA_2_perc
wilcox.test(SGA_AGA_LGA_2~as.factor(sam2$disease), data = sam2)
chisq.test(SGA_AGA_LGA_2_table)
SGA_AGA_LGA_2_perc
mcnemar.test(SGA_AGA_LGA_2_perc)

mantelhaen.test(SGA_AGA_LGA_2_perc)
SGA_perc<-tally(SGA_AGA_LGA_2=="AGA"~disease,data =sam2,format = "percent")
SGA_perc

GA_box<-ggplot(sam2, mapping = aes(y = newborn_weight,x=SGA_AGA_LGA_2,color=SGA_AGA_LGA_2))+
  geom_boxplot()+geom_jitter()+theme_bw()
  #geom_histogram(mapping = aes(color=SGA_AGA_LGA_2),position = "dodge",bins = 10)
aov(SGA_AGA_LGA_2~disease,sam3)
GA_box
GA_hist<-ggplot(sam2,  mapping = aes(x = newborn_weight))+
  geom_point()+
  #geom_dotplot(mapping = aes(color=SGA_AGA_LGA_2,fill=SGA_AGA_LGA_2))+
  geom_label(mapping = aes(y=newborn_weight),label=pt)+theme_bw()
GA_hist
#mapping = aes(x = newborn_weight,y=stat(),)
####################OLD WAY##################################
#############################################################
# sam2<-tibble(sam)%>%group_by(patient,disease)%>%distinct(Antibiotics_2)
# sam2<-filter(sam2, !is.na(Antibiotics_2))
# Antibiotics_2_perc<-tally(Antibiotics_2=="Abx"~disease,data =sam2,format = "percent")
# Antibiotics_2_table <- tally(Antibiotics_2=="Abx"~disease, data=sam2)
# Antibiotics_2_table
# Antibiotics_2_perc
# chisq.test(Antibiotics_2_perc)
# ##################################Delivery###################
# sam2<-tibble(sam)%>%group_by(patient,disease)%>%distinct(Delivery)
# sam2<-filter(sam2, !is.na(Delivery))
# Delivery_perc<-tally(Delivery=="Vaginal"~disease,data =sam2,format = "percent")
# Delivery.table <- tally(Delivery=="Vaginal"~disease, data=sam2)
# Delivery.table
# Delivery_perc
# chisq.test(Delivery_perc)