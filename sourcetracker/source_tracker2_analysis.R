library(tidyverse)
setwd("k:/github/Type_1_Diabetes_public/sourcetracker/")
df<-as_tibble(read.table("k:/github/Type_1_Diabetes_public/sourcetracker/sourcetracker_results.txt",header = T,sep = "\t"))
df
df2<-df%>%
  select(Study,SampleType,Delivery,Antibiotics,
         Anus,Cervix,Introitus,Vagina,Unknown)#,
         #X2_Anus,X2_Cervix,X2_Introitus,X2_Vagina,X2_Unknown)

df2
df2
df2<-df2%>%pivot_longer(cols =-c(Study,SampleType,Delivery,Antibiotics),names_to = "Source",values_to = "contribution")

df2<-df2%>%rename(disease=Study)%>%rename(Sink=SampleType)


library(ggpubr)
my_comparisons <- list(c("Control","T1D"))

ggboxplot(data = df2,
          x = "disease",
          y = "contribution",
          color="disease",
          add="jitter",
          palette = "npg",
          title = "sourcetracker2 derived maternal source contibutions to neonatal microbiomes")+
  facet_grid(Sink+Delivery~Source)+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox")+
  yscale("log10")+
  ylim(c(0,1.1))

df3<-df2%>%group_by(disease,Sink,Delivery,Source)%>%
  summarize(mean_contribution=mean(contribution),
            sd_contribution=sd(contribution))
df2
write.table(df3,"k:/github/Type_1_Diabetes_public/sourcetracker/contribution_overview.tsv",
            sep = "\t")
df2<-df[,17:21]

df2<-df2%>%mutate(across(where(is.character),factor))


library(mosaic)
l<-lm(contribution~disease*Delivery*Source*Sink,df2)
msummary(anova(l))


library(broom)
l2<-broom::tidy(l)
write.table(l2,file = "lm.tsv",sep = "\t")
msummary(TukeyHSD(l))