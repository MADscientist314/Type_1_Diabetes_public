library(tidyverse)


#############preprocessing####################
getwd()
#lets take a look at the outputs
load("FEAST_res_T1D_list.RData")
load("FEAST_res_control_list.RData")

#lets write the raw lists out to a file for marzena
sink("Feast_res_control_list.tsv")
print(FEAST_res_control_list)
sink()

sink("Feast_res_T1d_list.tsv")
print(FEAST_res_T1D_list)
sink()


#############preprocessing####################

#lets take a look at the outputs
FEAST_res_control_list
FEAST_res_T1D_list
getwd()

#lets write the raw lists out to a file for marzena
sink("Feast_res_control_list.tsv")
print(FEAST_res_control_list)
sink()

sink("Feast_res_T1d_list.tsv")
print(FEAST_res_T1D_list)
FEAST_res_T1D_list
sink()


length(FEAST_res_T1D_list)
#lets concatenate the FEAST outputs into a single mat for each analysis
control_mat<-tibble()
T1d_mat<-tibble()
control_meta<-tibble()
T1d_meta<-tibble()

#unlist the T1d list
for (i in 1:length(FEAST_res_T1D_list))
{
  k=FEAST_res_T1D_list[i]
  #get the matrices
  if (i %% 2 ==1) 
    T1d_mat<-bind_rows(k,T1d_mat, .id = NULL)
  #get the meta lists
  if (i %% 2 ==0) 
    print(j)
  T1d_meta <- bind_rows(k,T1d_meta, .id = NULL)
}
for (i in 1:length(FEAST_res_control_list))
{
  j=FEAST_res_control_list[i]
  #get the matrices
  if (i %% 2 ==1) 
    control_mat <- bind_rows(j,control_mat, .id = NULL)
  #get the meta lists
  if (i %% 2 ==0) 
    print(j)
  control_meta <- bind_rows(j,control_meta, .id = NULL)
}

######################################################################################
#this section of code takes the unlisted matrix output and generates a tibble with
#the outputs of disease,patient,SampleType,Source,contribution, and contribution to the sink

#add a disease column and rename the rownames as the Sinks
T1d_mat$disease<-"T1D"
T1d_mat$Sink<-rownames(T1d_mat)
#strip the patient names from the sink and tidy up the contribution column
T1d<-T1d_mat%>%
  pivot_longer(cols = -c(disease,Sink),names_to="Source",values_to="contribution")%>%
  filter(contribution!="NA")%>%
  separate(col=Sink,into = c("patient","SampleType"),sep = "_",remove = F)
#fix the pt nomenclature
T1d$patient<-sub(pattern="c.*",replacement ="", x = T1d$patient)
T1d
#Strip the SampleType from the Sink name, and regex the NAs to "Unknown"
T1d<-T1d%>%
  separate(col =Source,sep ="_",into = c("pt","Source"))%>%
  select(-c(pt,Sink))%>%
  mutate_at(c(4),~replace(., is.na(.), "Unknown"))
T1d

######################################################################################
#do the same thing for  the control samples
control_mat$disease<-"control"
control_mat$Sink<-rownames(control_mat)
control<-control_mat%>%
  pivot_longer(cols = -c(disease,Sink),names_to="Source",values_to="contribution")%>%
  filter(contribution!="NA")%>%
  separate(col=Sink,into = c("patient","SampleType"),sep = "_",remove = F)
control
#fix the pt nomenclature
control$patient<-sub(pattern="k.*",replacement ="", x = control$patient)
control
control<-control%>%
  separate(col =Source,sep ="_",into = c("pt","Source"))%>%
  select(-c(pt,Sink))%>%
  mutate_at(c(4),~replace(., is.na(.), "Unknown"))
control
######################################################################################
#ok lets merge them together
df<-full_join(T1d,control)
df
#and lets write this to a file for supplementary purposes
write.table(df,"Feast_res.tsv",sep = "\t",row.names = F)
#and lets make tomse figures and tables
library(ggpubr)
mypal<-c("#90CBF2","#B0E7A7")
my_comparisons <- list( c("control", "T1D"))
ggviolin(data = df,
         add = "jitter",
         draw_quantiles = 0.5,
         ggtheme = theme_pubr(),
         x ="disease",
         y = "contribution",
         fill = "disease",
         color="grey25",
         alpha = 1,
         #add = "jitter",
         palette = mypal)+
  facet_wrap(facets =SampleType~Source,nrow = 2,scales = "free")+
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test",
                     hide.ns = F,
                     na.rm = T,
                     label ="p.format",
                     label.x.npc = "top")
library(mosaic)
library(ggsci)
gf_boxplot(contribution~Source|SampleType,fill = ~disease,data = df)+
  scale_fill_manual(values = mypal)+
  stat_compare_means(method = "t.test",label ="p.format",label.y = 1.02)+
  theme_bw()

######################################################################################
#STATS ANALYSIS
######################################################################################

favstats(~mean_contribution | disease+SampleType+Source,
         data = df)

#Compute summary statistics by group
df2<-df%>%group_by(disease, SampleType,Source) %>%
  summarize(mean_contribution =
              mean(contribution))
gf_barh(data =df2,group = disease,color = Source)

gf_barsh(~contribution | disease+SampleType,
         data = df,color=~Source,position = position_dodgev())+scale_color_npg()
gf_percents(~contribution|disease,
            data = df,color=~Source,position = position_dodgev())+scale_color_npg()



gf_qq(~ contribution | disease,
      data = df) %>%
  gf_qqline() %>%
  gf_labs(x = "Normal quantile",
          y = "contribution to SampleType")
#Two-sample t-test and confidence interval
result <- t_test(contribution ~ disease,
                 data = df)
result # view results
confint(result)
pval(result)
#More than two levels (Analysis of variance)
#Numeric and graphic summaries

mod<-lm(contribution ~ Source*disease+SampleType,
        data = df)
lm(contribution ~ Source*disease+SampleType,
   data = df) %>%anova()
#mplot(anova(mod))

#Which differences are significant?
tuk<-TukeyHSD(mod)
mplot(object = tuk,level=0.95,system = "lattice")

save.image(file = "FEAST_analysis.RData")

