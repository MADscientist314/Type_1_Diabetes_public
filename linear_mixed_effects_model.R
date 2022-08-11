library(tidyverse)
library(ggpubr)
df
df<-as_tibble(read.table("sourcetracker_results2.txt",header = T,sep = "\t"))%>%
  rename(patient=Description)%>%
  rename(disease=Study)%>%
  rename(Sink=SampleType)%>%
  mutate(across(where(is.character),factor))


df
df2<-as_tibble(read.table("sourcetracker_results2.txt",header = T,sep = "\t"))%>%
  rename(patient=Description)%>%
  rename(disease=Study)%>%
  rename(Sink=SampleType)%>%
  select(SampleID,patient,disease,Sink,Delivery,Antibiotics,
         Anus,Cervix,Introitus,Vagina,Unknown)%>%
  pivot_longer(cols =-c(SampleID,patient,disease,Sink,Delivery,Antibiotics),
               names_to = "Source",
               values_to = "contribution")%>%mutate(across(where(is.character),factor))

df%>%pivot_wider(names_from = Sink,values_from=contribution)

df3<-df2%>%#select(-SampleID)%>%
  group_by(SampleID,patient,disease,Delivery,Antibiotics,Source,Sink)%>%
  summarise(contribution=mean(contribution))%>%
  pivot_wider(names_from =Sink,values_from=contribution)%>%ungroup()
df3


df2
library(lme4)
stool_lme<-glmer(Stool ~ disease + Delivery + Antibiotics+(1 | patient), data = df3, family = binomial )
ear_lme<-glmer(Ear ~ disease + Delivery + Antibiotics+(1 | patient), data = df3, family = binomial )
source_lme<-glmer(Source~Ear+Stool+ disease + Delivery + Antibiotics+(1 | SampleID), data = df3, family = binomial )

summary(stool_lme)
summary(ear_lme)
df2

a<-glmer(contribution ~ disease + Sink + Source + Antibiotics+(1|patient), data = df2, family = binomial )

b<-glm(Stool ~ disease + Delivery + Antibiotics, data = df3, family = binomial )
b<-glm(disease ~Stool+Ear, data = df3, family = binomial )
summary(b)
summary(aov(b))
qb<-glm(Stool ~ disease + Delivery + Antibiotics, data = df3, family = quasibinomial )
g<-glm(Stool ~ disease + Delivery + Antibiotics, data = df3, family = gaussian )
q<-glm(Stool ~ disease + Delivery + Antibiotics, data = df3, family = quasi )

## GLMM with individual-level variability (accounting for overdispersion)
## For this data set the model is the same as one allowing for a period:herd
## interaction, which the plot indicates could be needed.
df$disease_SampleType

cbpp$obs <- 1:nrow(cbpp)
#gm <- lm(Anus*Cervix~disease_SampleType,data=df)
df

t1d_lm <- lm(cbind(Anus,Unknown,Cervix,Introitus,Vagina)~disease*Sink*Delivery+Antibiotics,data=df)
t1d_lm
t1d_lm_tidy<-broom::tidy(t1d_lm)%>%filter(p.value<0.05)%>%
  filter(term!="(Intercept)")%>%
  mutate(term=gsub("Sink","",term),
         term=gsub("disease","",term),
         term=gsub("Delivery","",term),)

aov_t1d_lm<-aov(t1d_lm)
summary(aov_t1d_lm)
TukeyHSD(aov_t1d_lm)


t1d_lm_aov<-broom::tidy(aov_t1d_lm)
t1d_lm_aov

df
t1d_glm<-glm(Anus+Unknown+Cervix+Introitus+Vagina~disease_SampleType*Delivery,data=df)
t1d_anova<-aov(Anus+Unknown+Cervix+Introitus+Vagina~disease+Sink+Delivery+Antibiotics,data=df)
t1d_anova<-aov(t1d_glm)
summary(t1d_anova)
t1d_manova <- lm(cbind(Anus,Unknown,Cervix,Introitus,Vagina)~disease*Sink*Delivery,data=df)


broom::tidy(t1d_manova)%>%filter(p.value<0.05)

msummary(t1d_glm)


TukeyHSD(t1d_anova)
broom::tidy(t1d_manova)%>%filter(p.value<0.05)

t1d_glm

analysis <- Anova(t1d_lm, idata = df, idesign = ~SampleID)
print(analysis)
t1d_lm
gm
s<-broom::tidy(gm)
s
t<-broom::tidy(TukeyHSD(gm))
t%>%filter(adj.p.value<0.05)
msummary(gm$contrasts)
library(broom)
library(MASS)
fit <- lda(factor(disease_SampleType) ~  Anus+Unknown+Cervix+Introitus+Vagina,data=df)



gm2 <- glmer(contribution~disease*Delivery*Antibiotics*Source*Sink* +(1 | SampleID),family = binomial,data=df2,)
summary(gm2)
#gm2 <- glmer(disease~contribution*Source*Sink* +(1 | SampleID),family = binomial,data=df2,)

bin<-tally(disease~Delivery,"data.frame",data=df)
bin
binom.test(175,365)
binom.test(190,365)
binom.test(265,430)
binom.test(165,430)
binom.test(175,440)

175+265

175+190
265+165
msummary(gm2)
print(gm2)
summary(a)
aov(a)
One-way repeated measures ANOVA
require(foreign)
require(car)
model <- glmer(disease_SampleType~Anus+(1|patient), data = df,family = binomial)
summary(model)
analysis <- Anova(gm, idata = df2, idesign = ~patient)
print(analysis)