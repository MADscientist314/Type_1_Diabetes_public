sam<-sample_data(ASV_physeq_core)

library(dplyr)
library(broom)
library(mosaic)
library(tidyverse)
library(magrittr)
sam
categorical<-c('disease',	
               'Host',	
               'SampleType',	
               'Delivery',	
               'Antibiotics',
               'Antibiotics_2',	
               'antibiotic_csec_prophylaxis',	
               'antibiotic_S_agalactiae',	
               'antibiotic_Premature_rupture ')

# one categorical variable
#Counts by category
sam<-summarize(sam, .groups = patient)
sam
tally(~ disease, data =sam)
#Percentages by category
tally(~ disease,data = sam, format= "percent")
#Bar graph of percentages
gf_percents(Delivery ~ Antibiotics_2,data = q1, 
            fill = NULL,color = NULL)

#Exact test
result1 <- binom.test(~ disease =="disease", 
                      data = sam)
#Approximate test (large samples)
result2 <- prop.test(~ disease =="disease", 
                     data = sam, alternative = "less",
                     p = 0.4)
#Extract confidence intervals and p values
confint(result1)
pval(result2)


#################################################
#########2 catergorical variables#################
q <- sam %>% summarise(patient)
sam<-data.frame(sam)
q <- sam %>% filter(!is.na(sam$Delivery))
q
#Contingency table with margins

tally(~disease+Delivery,data =q , margins = TRUE)
#Percentages by column
tally(~ disease | Delivery,data =sam,format = "percent")
#Mosaic plot
my_tbl<-tally(disease ~ Delivery,data =sam)

mosaicplot(my_tbl , color = TRUE)

#Test for proportions (approximate)
prop.test(Delivery~disease,
          success = "Vaginal",
          data =q)
