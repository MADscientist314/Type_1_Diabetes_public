
################################################################################
# missing 1st trimester A1c value imputations
# objective: impute missing values for 1st trimester hemologbin A1C
# and make a normal distribution of values for the controls
################################################################################

library(caret)
mom<-as_tibble(sam)%>%
  select(patient,
         disease,
         Host,
         tydzien_porodu_2,
         newborn_weight,
         energ_tluszcz_proc_2,
         energ_bialko_proc_2,
         HbA1C_1trym_2,
         HbA1C_2trym_2,
         HBA1C_POROD_2,
         BMI_after_2,
         BMI_before_2,
         insulin_age,
         insulin_how_long)%>%
  filter(Host=="Mother")%>%
  distinct_all()
mom
mom<-mom%>%mutate(across(where(is.character),factor))%>%
  mutate(patient=as.factor(patient))
mom_T1d<-mom%>%filter(disease!="T1D")%>%select(patient,HbA1C_1trym_2)%>%distinct_all

favstats(mom_T1d$HbA1C_1trym_2)
knn_preproc<-preProcess(mom_T1d,method = c("knnImpute",
                                           #"bagImpute",
                                           #"medianImpute",
                                           #"zv",
                                           "nzv"#,
                                           #"conditionalX"
))
center_preproc<-preProcess(mom_T1d,method = c("center","scale"),verbose = T,outcome = "disease")

preproc<-predict(center_preproc,mom_T1d_not_missing)
mom_T1d_train<-preproc%>%slice_sample(prop = 0.9,weight_by = patient)
mom_T1d_train
mom_T1d_test<-preproc%>%filter(!patient%in%mom_T1d_train$patient)

mom_T1d_test

knn_preproc<-preProcess(mom_T1d_train,method = c("knnImpute"),verbose = T,outcome = "disease")

mom_T1d_test2<-predict(object = knn_preproc,newdata = mom_T1d_test)

# bag_preproc<-preProcess(mom_T1d,method = c("bagImpute"))
# median_preproc<-preProcess(mom_T1d,method = c("medianImpute"))
# 

# mom_T1d_knn<-predict(object = preproc,newdata = mom_T1d_test)

postResample(obs = mom_T1d_test$HbA1C_1trym_2,pred = mom_T1d_test2$HbA1C_1trym_2)

mom_T1d_bag<-predict(object = preproc,newdata = mom_T1d)
mom_T1d_median<-predict(object = preproc,newdata = mom_T1d)


mom_T1d<-mom_T1d%>%
  mutate(knn=mom_T1d_knn$HbA1C_1trym_2,
         bag=mom_T1d_bag$HbA1C_1trym_2,
         median=mom_T1d_median$HbA1C_1trym_2)

model <- lm(HbA1C_1trym_2~knn+bag+median, 
            data = mom_T1d%>%select(where(is.numeric)))
msummary(model)
#Prediction
lm_fun <- makeFun(model)
mom_T1d_missing<-mom_T1d%>%filter(is.na(HbA1C_1trym_2))
mom_T1d_not_missing<-mom_T1d%>%filter(!is.na(HbA1C_1trym_2))%>%
  mutate(imputed="actual")

mom_T1d_missing<-mom_T1d_missing%>%
  #  group_by(patient,disease)%>%
  rowwise()%>%
  mutate(HbA1C_1trym_2=lm_fun(knn = knn,bag = bag,median = median))%>%
  mutate(imputed="imputed")
summary(knn_preproc)

mom_T1d_2<-full_join(mom_T1d_missing,mom_T1d_not_missing)%>%
  mutate(imputed=factor(imputed))
show_col(pal_aaas()(10))
aaas<-pal_aaas()(10)
my_pal<-c(aaas[6],aaas[10])

p<-ggscatter(data = mom_T1d_2,
             x = "HbA1C_1trym_2",
             y = 0,repel = T,#add.params = list(max.overlaps=10),
             group = "patient",color = "imputed",palette = my_pal)+
  ggrepel::geom_label_repel(label=mom_T1d_2$patient,mapping = aes(color=imputed),
                            max.overlaps = 100)

ggpar(p = p,main = "Distibution of 1st trimester HBA1C values",
      caption = "missing values were imputed using a 
      k nearest neighbors (knn) algorithm (RMSE=0.032,Rsquared=1.0,MAE=0.02)",
      legend = "right")

################################################################################
