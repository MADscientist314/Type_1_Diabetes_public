library(dplyr)
library(ggpubr)
# Add p-values comparing groups
# Specify the comparisons you want
# Load data
data("ToothGrowth")
df <- ToothGrowth
head(df, 4)
#>    len supp dose
#> 1  4.2   VC  0.5
#> 2 11.5   VC  0.5
#> 3  7.3   VC  0.5
#> 4  5.8   VC  0.5

# Box plots with jittered points
# :::::::::::::::::::::::::::::::::::::::::::::::::::
# Change outline colors by groups: dose
# Use custom color palette
# Add jitter points and change the shape by groups
sam<-read.csv("vst_adonis_08102020_BACKUP/sample_coreAug7.txt", header =T, row.names = 1, sep = "\t")
sam<-data.frame(sam)
sam<-distinct(select(data.frame(sam),c(patient, disease,SampleType,tydzien_porodu_2)))

p <- ggboxplot(sam, x = "disease", y = "tydzien_porodu_2",
               color = "disease", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter", shape = NULL)
p
my_comparisons <- list( c("T1D", "Control") )
p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)+facet_wrap(~SampleType)  
# Violin plots with box plots inside
# :::::::::::::::::::::::::::::::::::::::::::::::::::
# Change fill color by groups: dose
# add boxplot with white fill color
ggviolin(sam, x = "disease", y = "tydzien_porodu_2", fill = "disease",
         palette = c("#00AFBB", "#E7B800"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add significance levels
  stat_compare_means(label.y = 50)   
dim(sam)
