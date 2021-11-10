library(mosaic)
library(dplyr)
library(ggpubr)
library(ggsci)
colnames(sam)
############AGE##############################
sam<-read.table(file = "sample_core_Sept16_2020.txt", header = T, sep = "\t", row.names = 1)
sam_mat<-as_tibble(sam)%>%filter(Host=="Mother")%>%distinct(disease,Host,age_mother,patient)
t_test(age_mother~disease, data = sam_mat)
wilcox.test(age_mother~disease, data = sam_mat)
gghistogram(data = sam_mat, x = "age_mother", y = "..count..", color = "disease")
ggdensity(data = sam_mat, x = "age_mother", y = "..count..", color = "disease")
favstats(age_mother~disease, data = sam_mat)
my_comparisons<-list(c("T1D", "Control"))
ggviolin(sam_mat, x = "disease", y = "age_mother", fill = "disease",
         palette = "npg",alpha=0.3,
         add = "jitter", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label.y=49)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "t.test",label.y = 50) 

############AGE##############################
sam_mat<-as_tibble(sam)%>%filter(Host=="Mother")%>%distinct(disease,Host,age_at_delivery_full_years,patient)
t_test(age_at_delivery_full_years~disease, data = sam_mat)
wilcox.test(age_at_delivery_full_years~disease, data = sam_mat)
gghistogram(data = sam_mat, x = "age_at_delivery_full_years", y = "..count..", color = "disease", fill="disease",alpha=0.4)
ggdensity(data = sam_mat, x = "age_at_delivery_full_years", y = "..density..", color = "disease", fill="disease", alpha = 0.4)
favstats(age_at_delivery_full_years~disease, data = sam_mat)
my_comparisons<-list(c("T1D", "Control"))
ggviolin(sam_mat, x = "disease", y = "age_at_delivery_full_years", fill = "disease",
         palette = "npg",alpha=0.3,
         add = "jitter", add.params = list(fill = "white"))+
  stat_compare_means(method = "wilcox.test",comparisons = my_comparisons, label.y=47)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox.test",label.y = 50) 

############Ethnicity##############################
sam_mat<-as_tibble(sam)%>%filter(Host=="Mother")%>%distinct(disease,Host,Ethnicity,patient)
Ethnicity<-table(tally(sam_mat$disease),tally(sam_mat$Ethnicity))

my_tbl<-tally(~disease,data =sam_mat)
my_tbl
chisq.test(tally(~disease,data =sam_mat))

############Delivery##############################
props<-c(50,42)/92
props
Caesarian<-c(33,20)
Vaginal<-c(12,17)
Vacuum<-c(5,5)
chisq.test(c(33,20),p = c(50,42)/92)
chisq.test(c(12,17),p =c(50,42)/92)
chisq.test(c(5,5),p = c(50,42)/92)

chisq.test(c(7,5),p= c(50,42)/92)
chisq.test(c(1,2),p= c(50,42)/92)
chisq.test(c(42,35),p= c(50,42)/92)

#############Gestation#########################
sam_mat<-as_tibble(sam)%>%filter(Host=="Mother")%>%distinct(disease,Host,Gestation_week,patient)
sam_mat
t_test(Gestation_week~disease, data = sam_mat)
wilcox.test(x = sam_mat$Gestation_week, y = tally(sam_mat$disease), data = sam_mat)
wilcox.test(Gestation_week~disease, data = sam_mat)

gghistogram(data = sam_mat, x = "Gestation_week", y = "..count..", color = "disease", fill="disease",alpha=0.4)
ggdensity(data = sam_mat, x = "Gestation_week", y = "..density..", color = "disease", fill="disease", alpha = 0.4)
favstats(Gestation_week~disease, data = sam_mat)
my_comparisons<-list(c("T1D", "Control"))
ggviolin(sam_mat, x = "disease", y = "Gestation_week", fill = "disease",
         palette = "npg",alpha=0.3,
         add = "jitter", add.params = list(fill = "white"))+
  stat_compare_means(method = "wilcox.test",comparisons = my_comparisons, label.y=47)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox.test",label.y = 50) 
###############Parity######################
chisq.test(c(28,14),p= c(50,42)/92)
chisq.test(c(22,27),p= c(50,42)/92)
chisq.test(c(0,1),p= c(50,42)/92)

########Misscarriages#########
chisq.test(c(7,8),p= c(50,42)/92)
chisq.test(c(41,33),p= c(50,42)/92)
chisq.test(c(2,1),p= c(50,42)/92)

#####Pre-pregnancy BMI#########
sam$pre_pregnancy_BMI
sam_mat<-as_tibble(sam)%>%filter(Host=="Mother")%>%distinct(disease,Host,pre_pregnancy_BMI,patient)
sam_mat
t_test(pre_pregnancy_BMI~disease, data = sam_mat)
wilcox.test(pre_pregnancy_BMI~disease, data = sam_mat)

gghistogram(data = sam_mat, x = "pre_pregnancy_BMI", y = "..count..", color = "disease", fill="disease",alpha=0.4)
ggdensity(data = sam_mat, x = "pre_pregnancy_BMI", y = "..density..", color = "disease", fill="disease", alpha = 0.4)
favstats(pre_pregnancy_BMI~disease, data = sam_mat)
my_comparisons<-list(c("T1D", "Control"))
ggviolin(sam_mat, x = "disease", y = "pre_pregnancy_BMI", fill = "disease",
         palette = "npg",alpha=0.3,
         add = "jitter", add.params = list(fill = "white"))+
  stat_compare_means(method = "wilcox.test",comparisons = my_comparisons, label.y=47)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox.test",label.y = 50)

###########Pre-pregnancy category#########
chisq.test(c(1,1),p= c(50,42)/92)
chisq.test(c(35,33),p= c(50,42)/92)
chisq.test(c(13,7),p= c(50,42)/92)
chisq.test(c(1,0),p= c(50,42)/92)
chisq.test(c(0,1),p= c(50,42)/92)



#####Pre-pregnancy BMI#########
sam$pre_delivery_BMI
sam_mat<-as_tibble(sam)%>%filter(Host=="Mother")%>%distinct(disease,Host,pre_delivery_BMI,patient)
sam_mat
t_test(pre_delivery_BMI~disease, data = sam_mat)
wilcox.test(pre_delivery_BMI~disease, data = sam_mat)

gghistogram(data = sam_mat, x = "pre_delivery_BMI", y = "..count..", color = "disease", fill="disease",alpha=0.4)
ggdensity(data = sam_mat, x = "pre_delivery_BMI", y = "..density..", color = "disease", fill="disease", alpha = 0.4)
favstats(pre_delivery_BMI~disease, data = sam_mat)
my_comparisons<-list(c("T1D", "Control"))
ggviolin(sam_mat, x = "disease", y = "pre_delivery_BMI", fill = "disease",
         palette = "npg",alpha=0.3,
         add = "jitter", add.params = list(fill = "white"))+
  stat_compare_means(method = "wilcox.test",comparisons = my_comparisons, label.y=47)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox.test",label.y = 50)

###########Weight###########
sam_mat<-as_tibble(sam)%>%filter(Host=="Mother")%>%distinct(disease,Host,weight_gain,patient)
sam_mat
t_test(weight_gain~disease, data = sam_mat)
wilcox.test(weight_gain~disease, data = sam_mat)

gghistogram(data = sam_mat, x = "weight_gain", y = "..count..", color = "disease", fill="disease",alpha=0.4)
ggdensity(data = sam_mat, x = "weight_gain", y = "..density..", color = "disease", fill="disease", alpha = 0.4)
favstats(weight_gain~disease, data = sam_mat)
my_comparisons<-list(c("T1D", "Control"))
ggviolin(sam_mat, x = "disease", y = "weight_gain", fill = "disease",
         palette = "npg",alpha=0.3,
         add = "jitter", add.params = list(fill = "white"))+
  stat_compare_means(method = "wilcox.test",comparisons = my_comparisons, label.y=47)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox.test",label.y = 50)


###########Weight gain###########
chisq.test(c(32,26),p= c(50,42)/92)
chisq.test(c(18,15),p= c(50,42)/92)
chisq.test(c(0,1),p= c(50,42)/92)

############Sat fat########
sam$energy_SFA
sam_mat<-as_tibble(sam)%>%filter(Host=="Mother")%>%filter(!is.na(energy_SFA))%>%distinct(disease,Host,energy_SFA,patient)
sam_mat
t_test(energy_SFA~disease, data = sam_mat)
wilcox.test(energy_SFA~disease, data = sam_mat)

gghistogram(data = sam_mat, x = "energy_SFA", y = "..count..", color = "disease", fill="disease",alpha=0.4)
ggdensity(data = sam_mat, x = "energy_SFA", y = "..density..", color = "disease", fill="disease", alpha = 0.4)
favstats(energy_SFA~disease, data = sam_mat)
my_comparisons<-list(c("T1D", "Control"))
ggviolin(sam_mat, x = "disease", y = "energy_SFA", fill = "disease",
         palette = "npg",alpha=0.3,
         add = "jitter", add.params = list(fill = "white"))+
  stat_compare_means(method = "wilcox.test",comparisons = my_comparisons, label.y=47)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox.test",label.y = 50)
############fat########
sam$energy_fat
sam_mat<-as_tibble(sam)%>%filter(Host=="Mother")%>%filter(!is.na(energy_fat))%>%distinct(disease,Host,energy_fat,patient)
sam_mat
t_test(energy_fat~disease, data = sam_mat)
wilcox.test(energy_fat~disease, data = sam_mat)

gghistogram(data = sam_mat, x = "energy_fat", y = "..count..", color = "disease", fill="disease",alpha=0.4)
ggdensity(data = sam_mat, x = "energy_fat", y = "..density..", color = "disease", fill="disease", alpha = 0.4)
favstats(energy_fat~disease, data = sam_mat)
my_comparisons<-list(c("T1D", "Control"))
ggviolin(sam_mat, x = "disease", y = "energy_fat", fill = "disease",
         palette = "npg",alpha=0.3,
         add = "jitter", add.params = list(fill = "white"))+
  stat_compare_means(method = "wilcox.test",comparisons = my_comparisons, label.y=47)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox.test",label.y = 50)

############carb########
sam$energy_carbohydrates
sam_mat<-as_tibble(sam)%>%filter(Host=="Mother")%>%filter(!is.na(energy_carbohydrates))%>%distinct(disease,Host,energy_carbohydrates,patient)
sam_mat
t_test(energy_carbohydrates~disease, data = sam_mat)
wilcox.test(energy_carbohydrates~disease, data = sam_mat)

gghistogram(data = sam_mat, x = "energy_carbohydrates", y = "..count..", color = "disease", fill="disease",alpha=0.4)
ggdensity(data = sam_mat, x = "energy_carbohydrates", y = "..density..", color = "disease", fill="disease", alpha = 0.4)
favstats(energy_carbohydrates~disease, data = sam_mat)
my_comparisons<-list(c("T1D", "Control"))
ggviolin(sam_mat, x = "disease", y = "energy_carbohydrates", fill = "disease",
         palette = "npg",alpha=0.3,
         add = "jitter", add.params = list(fill = "white"))+
  stat_compare_means(method = "wilcox.test",comparisons = my_comparisons, label.y=47)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox.test",label.y = 50)


############protein########
sam$energy_protein
sam_mat<-as_tibble(sam)%>%filter(Host=="Mother")%>%filter(!is.na(energy_protein))%>%distinct(disease,Host,energy_protein,patient)
sam_mat
t_test(energy_protein~disease, data = sam_mat)
wilcox.test(energy_protein~disease, data = sam_mat)

gghistogram(data = sam_mat, x = "energy_protein", y = "..count..", color = "disease", fill="disease",alpha=0.4)
ggdensity(data = sam_mat, x = "energy_protein", y = "..density..", color = "disease", fill="disease", alpha = 0.4)
favstats(energy_protein~disease, data = sam_mat)
my_comparisons<-list(c("T1D", "Control"))
ggviolin(sam_mat, x = "disease", y = "energy_protein", fill = "disease",
         palette = "npg",alpha=0.3,
         add = "jitter", add.params = list(fill = "white"))+
  stat_compare_means(method = "wilcox.test",comparisons = my_comparisons, label.y=47)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox.test",label.y = 50)

#####h1baC#######
chisq.test(c(50,0),p= c(50,42)/92)
chisq.test(c(0,42),p= c(50,42)/92)


##########Age at Diagnosis############
# 
# sam$age_T1D_diagnosis
# sam_mat<-as_tibble(sam)%>%filter(Host=="Mother")%>%filter(!is.na(age_T1D_diagnosis))%>%distinct(disease,Host,age_T1D_diagnosis,patient)
# sam_mat
# t_test(~ age_T1D_diagnosis, data = sam_mat, mu=13.88)
# wilcox.test(age_T1D_diagnosis~disease, data = sam_mat)
# 
# gghistogram(data = sam_mat, x = "age_T1D_diagnosis", y = "..count..", color = "disease", fill="disease",alpha=0.4)
# ggdensity(data = sam_mat, x = "age_T1D_diagnosis", y = "..density..", color = "disease", fill="disease", alpha = 0.4)
# favstats(age_T1D_diagnosis~disease, data = sam_mat)
# my_comparisons<-list(c("T1D", "Control"))
# ggviolin(sam_mat, x = "disease", y = "age_T1D_diagnosis", fill = "disease",
#          palette = "npg",alpha=0.3,
#          add = "jitter", add.params = list(fill = "white"))+
#   stat_compare_means(method = "wilcox.test",comparisons = my_comparisons, label.y=47)+ # Add pairwise comparisons p-value
#   stat_compare_means(method = "wilcox.test",label.y = 50)


############birthweight########
sam$birth_weight
sam_mat<-as_tibble(sam)%>%filter(Host=="Infant")%>%filter(!is.na(birth_weight))%>%distinct(disease,Host,birth_weight,patient)
sam_mat
t_test(birth_weight~disease, data = sam_mat)
wilcox.test(birth_weight~disease, data = sam_mat)

gghistogram(data = sam_mat, x = "birth_weight", y = "..count..", color = "disease", fill="disease",alpha=0.4)
ggdensity(data = sam_mat, x = "birth_weight", y = "..density..", color = "disease", fill="disease", alpha = 0.4)
favstats(birth_weight~disease, data = sam_mat)
my_comparisons<-list(c("T1D", "Control"))
ggviolin(sam_mat, x = "disease", y = "birth_weight", fill = "disease",
         palette = "npg",alpha=0.3,
         add = "jitter", add.params = list(fill = "white"))+
  stat_compare_means(method = "wilcox.test",comparisons = my_comparisons, label.y=47)+ # Add pairwise comparisons p-value
  stat_compare_means(method = "wilcox.test",label.y = 50)

