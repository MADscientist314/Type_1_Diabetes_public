library(microbiome)
library(magrittr)
library(dplyr)
library(mosaic)
sam<-read.csv(file = "sample_coreJuly15.txt", header = T, sep = "\t", row.names = 1)
sam
##################################newborn_weight###################
sam2<-tibble(sam)%>%group_by(patient,disease)%>%distinct(newborn_weight)
#Two-sample t-test and confidence interval
capture.output(t_test(newborn_weight ~ disease,data = sam2),append = F, file = "08152020_analysis/newborn_weight_vs_disease.txt")
capture.output(favstats(x = newborn_weight~disease,data = sam2),append=T,file = "08152020_analysis/newborn_weight_vs_disease.txt")
tally(x = ~ disease, data = sam2)
write.table(sam2,file = "08152020_analysis/pt_disease_newborn_weight.tsv", sep = "\t")
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
# Plot
phist <- gghistogram(
  sam2, x = "newborn_weight", 
  add = "mean", rug = TRUE,
  fill = "disease", palette = c("#61D04F","#2297E6")
)
#> Warning: Using `bins = 30` by default. Pick better value with the argument
#> `bins`.

# Density plot
pdensity <- ggdensity(
  sam2, x = "newborn_weight", 
  color= "disease", palette = c("#61D04F","#2297E6"),
  alpha = 0, xlab = 
) +
  theme_half_open(11, rel_small = 1) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    position = "right"
  )  +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = grid::unit(0, "pt")
  )

# Aligning histogram and density plots
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
#######################################################################
##################################tydzien_porodu_2###################
sam2<-tibble(sam)%>%group_by(patient,disease)%>%distinct(tydzien_porodu_2)
#Two-sample t-test and confidence interval
capture.output(t_test(tydzien_porodu_2 ~ disease,data = sam2),append = F, file = "08152020_analysis/tydzien_porodu_2_vs_disease.txt")
capture.output(favstats(x = tydzien_porodu_2~disease,data = sam2),append=T,file = "08152020_analysis/tydzien_porodu_2_vs_disease.txt")
tally(x = ~ disease, data = sam2)
write.table(sam2,file = "08152020_analysis/pt_disease_tydzien_porodu_2.tsv", sep = "\t")
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
# Plot
phist <- gghistogram(
  sam2, x = "tydzien_porodu_2",facet.by = "disease", 
  add = "mean", rug = TRUE,xlab = "Gestation (weeks)",
  fill = "disease", palette = c("#61D04F","#2297E6")
)
#> Warning: Using `bins = 30` by default. Pick better value with the argument
#> `bins`.

# Density plot
pdensity <- ggdensity(
  sam2, x = "tydzien_porodu_2", 
  color= "disease", palette = c("#61D04F","#2297E6"),facet.by = "disease",
  alpha = 0.1, xlab = "Gestation (weeks)",fill = "disease")+
  theme_half_open(11, rel_small = 1) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    position = "right"
  )  +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = grid::unit(0, "pt")
  )

# Aligning histogram and density plots
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
################################################################
##################################BMI_before_2###################
sam2<-tibble(sam)%>%group_by(patient,disease)%>%distinct(BMI_before_2)
#Two-sample t-test and confidence interval
capture.output(t_test(BMI_before_2 ~ disease,data = sam2),append = F, file = "08152020_analysis/BMI_before_2_vs_disease.txt")
capture.output(favstats(x = BMI_before_2~disease,data = sam2),append=T,file = "08152020_analysis/BMI_before_2_vs_disease.txt")

write.table(sam2,file = "08152020_analysis/pt_disease_BMI_before_2.tsv", sep = "\t")
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
# Plot
phist <- gghistogram(
  sam2, x = "BMI_before_2",facet.by = "disease", 
  add = "mean", rug = TRUE,xlab = "Gestation (weeks)",
  fill = "disease", palette = c("#61D04F","#2297E6")
)
#> Warning: Using `bins = 30` by default. Pick better value with the argument
#> `bins`.

# Density plot
pdensity <- ggdensity(
  sam2, x = "BMI_before_2", 
  color= "disease", palette = c("#61D04F","#2297E6"),facet.by = "disease",
  alpha = 0.1, xlab = "Gestation (weeks)",fill = "disease")+
  theme_half_open(11, rel_small = 1) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    position = "right"
  )  +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = grid::unit(0, "pt")
  )

# Aligning histogram and density plots
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
################################################################

################################################################
##################################energ_bialko_proc_2###################
sam2<-tibble(sam)%>%group_by(patient,disease)%>%distinct(energ_bialko_proc_2)
#Two-sample t-test and confidence interval
capture.output(t_test(energ_bialko_proc_2 ~ disease,data = sam2),append = F, file = "08152020_analysis/energ_bialko_proc_2_vs_disease.txt")
capture.output(favstats(x = energ_bialko_proc_2~disease,data = sam2),append=T,file = "08152020_analysis/energ_bialko_proc_2_vs_disease.txt")

write.table(sam2,file = "08152020_analysis/pt_disease_energ_bialko_proc_2.tsv", sep = "\t")
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(cowplot))
# Plot
phist <- gghistogram(
  sam2, x = "energ_bialko_proc_2",facet.by = "disease", 
  add = "mean", rug = TRUE,xlab = "Gestation (weeks)",
  fill = "disease", palette = c("#61D04F","#2297E6")
)
#> Warning: Using `bins = 30` by default. Pick better value with the argument
#> `bins`.

# Density plot
pdensity <- ggdensity(
  sam2, x = "energ_bialko_proc_2", 
  color= "disease", palette = c("#61D04F","#2297E6"),facet.by = "disease",
  alpha = 0.1, xlab = "Gestation (weeks)",fill = "disease")+
  theme_half_open(11, rel_small = 1) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    position = "right"
  )  +
  theme(
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks = element_blank(),
    axis.ticks.length = grid::unit(0, "pt")
  )

# Aligning histogram and density plots
aligned_plots <- align_plots(phist, pdensity, align="hv", axis="tblr")
ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])
################################################################
sam$energ_bialko_proc_2<-as.numeric(sam$Delivery)
sample_data(ASV_physeq_core)<-sample_data(sam)
a<-subset_samples(ASV_physeq_core, Delivery!="NA")
plot_landscape(method = "PCoA", 
               distance = "bray",
               transformation ='compositional',
               col = "disease",size=2,
               x = a,shading = F)+
  facet_grid(facets = meta(a)$SampleType~meta(a)$Delivery)+scale_color_manual(values = c("#61D04F","#2297E6"))+theme_bw()








##################################Delivery###################
sam2<-tibble(sam)%>%group_by(patient,disease)%>%distinct(Delivery)
sam2<-filter(sam2, !is.na(Delivery))
Delivery_perc<-tally(Delivery~disease,data =sam2,format = "percent")
Delivery.table <- tally(Delivery~disease, data=sam2)
Delivery_perc
Delivery.table
fisher.test(Delivery.table)
capture.output(fisher.test(Delivery.table),append = F, file = "08152020_analysis/Delivery_vs_disease.txt")
capture.output(Delivery.table,append=T,file = "08152020_analysis/Delivery_vs_disease.txt")
capture.output(Delivery_perc,append=T,file = "08152020_analysis/Delivery_vs_disease.txt")

write.table(sam2,file = "08152020_analysis/pt_disease_Delivery.tsv", sep = "\t")
##########################################################################################
##################################Antibiotics_2###################
sam2<-tibble(sam)%>%group_by(patient,disease)%>%distinct(Antibiotics_2)
sam2<-filter(sam2, !is.na(Antibiotics_2))
Antibiotics_2_perc<-tally(Antibiotics_2~disease,data =sam2,format = "percent")
Antibiotics_2.table <- tally(Antibiotics_2~disease, data=sam2)
Antibiotics_2_perc
Antibiotics_2.table
fisher.test(Antibiotics_2.table)
capture.output(fisher.test(Antibiotics_2.table),append = F, file = "08152020_analysis/Antibiotics_2_vs_disease.txt")
capture.output(Antibiotics_2.table,append=T,file = "08152020_analysis/Antibiotics_2_vs_disease.txt")
capture.output(Antibiotics_2_perc,append=T,file = "08152020_analysis/Antibiotics_2_vs_disease.txt")

write.table(sam2,file = "08152020_analysis/pt_disease_Antibiotics_2.tsv", sep = "\t")
##########################################################################################
##################################abx_purpose###################
sam2<-tibble(sam)%>%group_by(patient,disease)%>%distinct(abx_purpose)
sam2<-filter(sam2, !is.na(abx_purpose))
abx_purpose_perc<-tally(abx_purpose~disease,data =sam2,format = "percent")
abx_purpose.table <- tally(abx_purpose~disease, data=sam2)
abx_purpose_perc
abx_purpose.table
a<-chisq.test(abx_purpose.table)
a$observed
a$expected
a$residuals
a$stdres
a$parameter
capture.output(fisher.test(abx_purpose.table),append = F, file = "08152020_analysis/abx_purpose_vs_disease.txt")
capture.output(abx_purpose.table,append=T,file = "08152020_analysis/abx_purpose_vs_disease.txt")
capture.output(abx_purpose_perc,append=T,file = "08152020_analysis/abx_purpose_vs_disease.txt")

write.table(sam2,file = "08152020_analysis/pt_disease_abx_purpose.tsv", sep = "\t")
##########################################################################################