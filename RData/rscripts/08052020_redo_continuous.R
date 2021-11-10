library(phyloseq)
# sam<-sample_data(ASV_physeq_core)

library(dplyr)
library(broom)
library(mosaic)
library(tidyverse)
mother_comp_list <- c("newborn_weight", "tydzien_porodu_2","energ_tluszcz_proc_2","energ_bialko_proc_2","Bialko_ogolem_g_2","BMI_before_2",	"BMI_after_2")
mother_T1d_only_list<-c( "HbA1C_1trym_2", "HbA1C_2trym_2", "HBA1C_POROD_2", "insulin_age",	"insulin_how_long")
newborn_only<-c("Glu1",	"Glu2",	"Glu3",	"Glu4","Glu5")
q <-sam %>% group_by(patient) %>%
  filter(!is.na(insulin_how_long))
gf_histogram(~ insulin_how_long,data =q , binwidth = 5,  center = 2.5)
###########BMI_before_2###############
##########tidy the data#########
q <-sam %>% group_by(patient) %>%
  filter(!is.na(BMI_before_2))
######LOG REGRESSION
#Logistic regression
logit_mod <- glm(q$BMI_before_2 ~ q$disease,
                 data = q)

capture.output(msummary(logit_mod), file = "BMI_before_2_log_regression_summary.txt")
#Odds ratios and confidence intervals
capture.output(tidy(logit_mod, conf.int = TRUE,exponentiate = TRUE), file = "BMI_before_2_OR_and_CI.txt")
#Odds ratios and confidence intervals
#Numeric summaries
capture.output(favstats(~ BMI_before_2 | disease, data = q),file = "BMI_before_2_fav_stats.txt")
#Graphic summaries
gf_qq(~BMI_before_2 | disease,
      data = q) %>%
  gf_qqline() %>%
  gf_labs(x = "Normal quantile",
          y = "BMI_before_2")#+
  ggsave(filename = "BMI_before_2_Normal_quantile.png", device = "png", height = 4, width = 4,units = "in" )

  #Two-sample t-test and confidence interval
capture.output(t_test(BMI_before_2 ~ disease,
                 data = q),file = "BMI_before_2_2_sample_t_test_w_CI.txt")

##############################################
#More than two levels (Analysis of variance)
############################################
#Numeric and graphic summaries
favstats(disease ~ BMI_before_2, data = q)
gf_boxplot(BMI_before_2 ~ disease,
           data = q)+
  ggsave(filename = "BMI_before_2_boxplot.png", device = "png", height = 4, width = 4,units = "in" )

#Fit and summarize model
mod <- lm(BMI_before_2 ~ disease,
          data = q)


a<-anova(mod)
a
capture.output(msummary(mod),file = "BMI_before_2_disease_anova_summary.txt")

#Which differences are significant?
capture.output(TukeyHSD(mod), file="BMI_before_2Tukey_HSD.txt")




##################newborn_weight##################
##########tidy the data#########

sam<-data.frame(sam)
q <-distinct(select(data.frame(sam),c(patient, disease, newborn_weight)))


######LOG REGRESSION
#Logistic regression
logit_mod <- glm(newborn_weight ~ disease,
                 data = q)
msummary(logit_mod)

capture.output(msummary(logit_mod), file = "newborn_weight_log_regression_summary.txt")
#Odds ratios and confidence intervals
capture.output(tidy(logit_mod, conf.int = TRUE,exponentiate = TRUE), file = "newborn_weight_OR_and_CI.txt")
#Odds ratios and confidence intervals
#Numeric summaries
capture.output(favstats(~ newborn_weight | disease, data = q),file = "newborn_weight_fav_stats.txt")
#Graphic summaries
gf_qq(~newborn_weight | disease,
      data = q) %>%
  gf_qqline() %>%
  gf_labs(x = "Normal quantile",
          y = "newborn_weight")+
  ggsave(filename = "newborn_weight_Normal_quantile.png", device = "png", height = 4, width = 4,units = "in" )
#Two-sample t-test and confidence interval
capture.output(t_test(newborn_weight ~ disease,
                      data = q),file = "newborn_weight_2_sample_t_test_w_CI.txt")

##############################################
#More than two levels (Analysis of variance)
############################################
#Numeric and graphic summaries
favstats(newborn_weight ~disease, data = q)
gf_boxplot(newborn_weight ~ disease,
           data = q)+
  ggsave(filename = "newborn_weight_boxplot.png", device = "png", height = 4, width = 4,units = "in" )

#Fit and summarize model
mod <- lm(newborn_weight ~ disease,
          data = q)


a<-anova(mod)
a

capture.output(msummary(mod),file = "newborn_weight_disease_anova_summary.txt")

#Which differences are significant?
capture.output(TukeyHSD(mod), file="newborn_weight_Tukey_HSD.txt")



###########tydzien_porodu_2###############
##########tidy the data#########
sam<-data.frame(sam)
q <-distinct(select(data.frame(sam),c(patient, disease, tydzien_porodu_2)))
dim(q)
# q <-sam %>% group_by(patient) %>%
#   filter(!is.na(tydzien_porodu_2))
######LOG REGRESSION
#Logistic regression
logit_mod <- glm(q$tydzien_porodu_2 ~ q$disease,
                 data = q)
msummary(logit_mod)
capture.output(msummary(logit_mod), file = "tydzien_porodu_2_log_regression_summary.txt")
#Odds ratios and confidence intervals
capture.output(tidy(logit_mod, conf.int = TRUE,exponentiate = TRUE), file = "tydzien_porodu_2_OR_and_CI.txt")
#Odds ratios and confidence intervals
#Numeric summaries
capture.output(favstats(~ tydzien_porodu_2 | disease, data = q),file = "tydzien_porodu_2_fav_stats.txt")
#Graphic summaries
gf_qq(~tydzien_porodu_2 | disease,
      data = q) %>%
  gf_qqline() %>%
  gf_labs(x = "Normal quantile",
          y = "tydzien_porodu_2")+
  ggsave(filename = "tydzien_porodu_2_Normal_quantile.png", device = "png", height = 4, width = 4,units = "in" )
#Two-sample t-test and confidence interval
capture.output(t_test(tydzien_porodu_2 ~ disease,
                      data = q),file = "tydzien_porodu_2_2_sample_t_test_w_CI.txt")

##############################################
#More than two levels (Analysis of variance)
############################################
#Numeric and graphic summaries
favstats(disease ~ tydzien_porodu_2, data = q)
gf_boxplot(tydzien_porodu_2 ~ disease,
           data = q)+
  ggsave(filename = "tydzien_porodu_2_boxplot.png", device = "png", height = 4, width = 4,units = "in" )

#Fit and summarize model
mod <- lm(tydzien_porodu_2 ~ disease,
          data = q)


a<-anova(mod)
a
capture.output(msummary(mod),file = "tydzien_porodu_2_disease_anova_summary.txt")

#Which differences are significant?
capture.output(TukeyHSD(mod), file="tydzien_porodu_2Tukey_HSD.txt")


##########energ_tluszcz_proc_2###############
##########tidy the data#########
q <-sam %>% group_by(patient) %>%
  filter(!is.na(energ_tluszcz_proc_2))
######LOG REGRESSION
#Logistic regression
logit_mod <- glm(q$energ_tluszcz_proc_2 ~ q$disease,
                 data = q)

capture.output(msummary(logit_mod), file = "energ_tluszcz_proc_2_log_regression_summary.txt")
#Odds ratios and confidence intervals
capture.output(tidy(logit_mod, conf.int = TRUE,exponentiate = TRUE), file = "energ_tluszcz_proc_2_OR_and_CI.txt")
#Odds ratios and confidence intervals
#Numeric summaries
capture.output(favstats(~ energ_tluszcz_proc_2 | disease, data = q),file = "energ_tluszcz_proc_2_fav_stats.txt")
#Graphic summaries
gf_qq(~energ_tluszcz_proc_2 | disease,
      data = q) %>%
  gf_qqline() %>%
  gf_labs(x = "Normal quantile",
          y = "energ_tluszcz_proc_2")+
  ggsave(filename = "energ_tluszcz_proc_2_Normal_quantile.png", device = "png", height = 4, width = 4,units = "in" )
#Two-sample t-test and confidence interval
capture.output(t_test(energ_tluszcz_proc_2 ~ disease,
                      data = q),file = "energ_tluszcz_proc_2_2_sample_t_test_w_CI.txt")

##############################################
#More than two levels (Analysis of variance)
############################################
#Numeric and graphic summaries
favstats(disease ~ energ_tluszcz_proc_2, data = q)
gf_boxplot(energ_tluszcz_proc_2 ~ disease,
           data = q)+
  ggsave(filename = "energ_tluszcz_proc_2_boxplot.png", device = "png", height = 4, width = 4,units = "in" )

#Fit and summarize model
mod <- lm(energ_tluszcz_proc_2 ~ disease,
          data = q)


a<-anova(mod)
a
capture.output(msummary(mod),file = "energ_tluszcz_proc_2_disease_anova_summary.txt")

#Which differences are significant?
capture.output(TukeyHSD(mod), file="energ_tluszcz_proc_2Tukey_HSD.txt")










###########energ_bialko_proc_2###############
##########tidy the data#########
q <-sam %>% group_by(patient) %>%
  filter(!is.na(energ_bialko_proc_2))
######LOG REGRESSION
#Logistic regression
logit_mod <- glm(q$energ_bialko_proc_2 ~ q$disease,
                 data = q)

capture.output(msummary(logit_mod), file = "energ_bialko_proc_2_log_regression_summary.txt")
#Odds ratios and confidence intervals
capture.output(tidy(logit_mod, conf.int = TRUE,exponentiate = TRUE), file = "energ_bialko_proc_2_OR_and_CI.txt")
#Odds ratios and confidence intervals
#Numeric summaries
capture.output(favstats(~ energ_bialko_proc_2 | disease, data = q),file = "energ_bialko_proc_2_fav_stats.txt")
#Graphic summaries
gf_qq(~energ_bialko_proc_2 | disease,
      data = q) %>%
  gf_qqline() %>%
  gf_labs(x = "Normal quantile",
          y = "energ_bialko_proc_2")+
  ggsave(filename = "energ_bialko_proc_2_Normal_quantile.png", device = "png", height = 4, width = 4,units = "in" )
#Two-sample t-test and confidence interval
capture.output(t_test(energ_bialko_proc_2 ~ disease,
                      data = q),file = "energ_bialko_proc_2_2_sample_t_test_w_CI.txt")

##############################################
#More than two levels (Analysis of variance)
############################################
#Numeric and graphic summaries
favstats(disease ~ energ_bialko_proc_2, data = q)
gf_boxplot(energ_bialko_proc_2 ~ disease,
           data = q)+
  ggsave(filename = "energ_bialko_proc_2_boxplot.png", device = "png", height = 4, width = 4,units = "in" )

#Fit and summarize model
mod <- lm(energ_bialko_proc_2 ~ disease,
          data = q)


a<-anova(mod)
a
capture.output(msummary(mod),file = "energ_bialko_proc_2_disease_anova_summary.txt")

#Which differences are significant?
capture.output(TukeyHSD(mod), file="energ_bialko_proc_2Tukey_HSD.txt")



###########Bialko_ogolem_g_2###############
##########tidy the data#########
q <-sam %>% group_by(patient) %>%
  filter(!is.na(Bialko_ogolem_g_2))
######LOG REGRESSION
#Logistic regression
logit_mod <- glm(q$Bialko_ogolem_g_2 ~ q$disease,
                 data = q)

capture.output(msummary(logit_mod), file = "Bialko_ogolem_g_2_log_regression_summary.txt")
#Odds ratios and confidence intervals
capture.output(tidy(logit_mod, conf.int = TRUE,exponentiate = TRUE), file = "Bialko_ogolem_g_2_OR_and_CI.txt")
#Odds ratios and confidence intervals
#Numeric summaries
capture.output(favstats(~ Bialko_ogolem_g_2 | disease, data = q),file = "Bialko_ogolem_g_2_fav_stats.txt")
#Graphic summaries
gf_qq(~Bialko_ogolem_g_2 | disease,
      data = q) %>%
  gf_qqline() %>%
  gf_labs(x = "Normal quantile",
          y = "Bialko_ogolem_g_2")+
  ggsave(filename = "Bialko_ogolem_g_2_Normal_quantile.png", device = "png", height = 4, width = 4,units = "in" )
#Two-sample t-test and confidence interval
capture.output(t_test(Bialko_ogolem_g_2 ~ disease,
                      data = q),file = "Bialko_ogolem_g_2_2_sample_t_test_w_CI.txt")

##############################################
#More than two levels (Analysis of variance)
############################################
#Numeric and graphic summaries
favstats(disease ~ Bialko_ogolem_g_2, data = q)
gf_boxplot(Bialko_ogolem_g_2 ~ disease,
           data = q)+
  ggsave(filename = "Bialko_ogolem_g_2_boxplot.png", device = "png", height = 4, width = 4,units = "in" )

#Fit and summarize model
mod <- lm(Bialko_ogolem_g_2 ~ disease,
          data = q)


a<-anova(mod)
a
capture.output(msummary(mod),file = "Bialko_ogolem_g_2_disease_anova_summary.txt")

#Which differences are significant?
capture.output(TukeyHSD(mod), file="Bialko_ogolem_g_2Tukey_HSD.txt")





###########BMI_after_2###############
##########tidy the data#########
q <-sam %>% group_by(patient) %>%
  filter(!is.na(BMI_after_2))
######LOG REGRESSION
#Logistic regression
logit_mod <- glm(q$BMI_after_2 ~ q$disease,
                 data = q)

capture.output(msummary(logit_mod), file = "BMI_after_2_log_regression_summary.txt")
#Odds ratios and confidence intervals
capture.output(tidy(logit_mod, conf.int = TRUE,exponentiate = TRUE), file = "BMI_after_2_OR_and_CI.txt")
#Odds ratios and confidence intervals
#Numeric summaries
capture.output(favstats(~ BMI_after_2 | disease, data = q),file = "BMI_after_2_fav_stats.txt")
#Graphic summaries
gf_qq(~BMI_after_2 | disease,
      data = q) %>%
  gf_qqline() %>%
  gf_labs(x = "Normal quantile",
          y = "BMI_after_2")+
  ggsave(filename = "BMI_after_2_Normal_quantile.png", device = "png", height = 4, width = 4,units = "in" )
#Two-sample t-test and confidence interval
capture.output(t_test(BMI_after_2 ~ disease,
                      data = q),file = "BMI_after_2_2_sample_t_test_w_CI.txt")

##############################################
#More than two levels (Analysis of variance)
############################################
#Numeric and graphic summaries
favstats(disease ~ BMI_after_2, data = q)
gf_boxplot(BMI_after_2 ~ disease,
           data = q)+
  ggsave(filename = "BMI_after_2_boxplot.png", device = "png", height = 4, width = 4,units = "in" )

#Fit and summarize model
mod <- lm(BMI_after_2 ~ disease,
          data = q)


a<-anova(mod)
a
capture.output(msummary(mod),file = "BMI_after_2_disease_anova_summary.txt")

#Which differences are significant?
capture.output(TukeyHSD(mod), file="BMI_after_2Tukey_HSD.txt")



favstats(disease ~ newborn_weight + tydzien_porodu_2 + energ_tluszcz_proc_2 + energ_bialko_proc_2 + Bialko_ogolem_g_2 + BMI_before_2 + BMI_after_2, data = q)



