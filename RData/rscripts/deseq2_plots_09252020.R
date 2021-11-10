library(phyloseq)
library(microbiome)
ASV_physeq_core
sam<-read.table(file = "sample_core_Sept16_2020.txt", header = T, sep = "\t", row.names = 1)
sam$age_mother<-as.factor(sam$age_mother)
sam$age_at_delivery_full_years<-as.factor(sam$age_at_delivery_full_years)
sam$Gestation_week<-as.factor(sam$Gestation_week)
sam$disease<-as.factor(sam$disease)
sam$SampleType<-as.factor(sam$SampleType)
sample_data(ASV_physeq_core)<-sam
ASV_physeq_core_genus<-tax_glom(physeq = ASV_physeq_core,taxrank = "Genus")
ASV_physeq_core_genus_rel<-microbiome::transform(ASV_physeq_core, "compositional")
deseq_res<-read.table(file = "DESEQ2_disease_VST_results.tsv", header = T, sep = "\t", row.names = 1)
deseq_res$Genus
taxa_names(ASV_physeq_core_genus)
ASV_physeq_core_genus_deseq2<-prune_taxa(x = ASV_physeq_core_genus, deseq_res$Genus)
ASV_physeq_core_genus_rel_deseq2
saveRDS(object = ASV_physeq_core_genus_rel_deseq2, file = "ASV_physeq_core_genus_rel_deseq2.rds")
save.image(file = "ASV_physeq_core_genus_rel_deseq2.rda")

top11<-psmelt(ASV_physeq_core_genus_rel_deseq2)
top11<-as_tibble(top11)                   
                   
                 
ASV_physeq_core_genus_rel_deseq2<-subset_samples(physeq = ASV_physeq_core_genus_rel_deseq2, sample_sums(ASV_physeq_core_genus_rel_deseq2)>0)
ASV_physeq_core_genus_rel_deseq2<-subset_taxa(physeq = ASV_physeq_core_genus_rel_deseq2, taxa_sums(ASV_physeq_core_genus_rel_deseq2)>0)
plot_composition(ASV_physeq_core_genus_rel_deseq2, sample.sort = "neatmap",otu.sort = "abundance",group_by = "disease")+
#  facet_grid(facets = ~"SampleType")+
  scale_fill_brewer("Paired")
plot_heatmap(physeq = ASV_physeq_core_genus_rel_deseq2, 
             method = "PCA", 
             distance = "bray", 
             sample.label = "disease", 
             taxa.label = "Genus")

microbiome::heat()

install.packages("shiny")
shiny::runGitHub("shiny-phyloseq","joey711")

# The mtcars dataset:
data <- as.matrix(top11$Abundance)
data
# Default Heatmap
heatmap(x = data)

a<-saveRDS(object = ASV_physeq_core_genus_rel_deseq2, file = "ASV_physeq_core_genus_rel_deseq2.rds")


d1 <- top11$disease
d2 <- top11$Abpeerj32$microbes[, 1:10]
cc <- associate(d1, d2, method='pearson') 
p <- heatmap(cc, 'X1', 'X2', 'Correlation', star='p.adj')
p
top11<-psmelt(ASV_physeq_core_genus_rel_deseq2)

library(dplyr)                  

top11<-as_tibble(top11)
top11<-group_by(top11)
top11$Abundancepi<-pi*(top11$Abundance*top11$Abundance)*10000000000
top11$Abundance
top11$Abundancepi
top1<-top11%>%group_by(disease,Genus,Abundance)
top1
library(ggpubr)
library(ggsci)
my_comparisons<-list(c("Control","T1D"))

top1<-top11 %>%
  group_by(disease,SampleType,Genus) %>%
  summarise(avg = mean(Abundance))
top1
top11$Host<-factor(x = top11$Host,levels = c("Mother","Infant"))

ggboxplot(data = top11, x = "disease",
          fill = "disease",
          y = "Abundance",
          add="jitter",
          facet.by = c("Genus"))+
  yscale(.scale = "log10",.format = T)+
  rotate_x_text()+
  scale_color_aaas()+
  scale_fill_aaas()+
  stat_compare_means(comparisons = my_comparisons,label.y.npc = 0.9, hide.ns = T)

ggballoonplot(data = top1,size.range = c(10,25),#c(min(top1$log10avg),max(25),
              x = "disease",
              y = "Genus",
              size = "log10avg",
              color = NULL,
              fill="avg",
              ggtheme = theme_bw())+
  scale_fill_gradient(
    low = "darkred",
    high = "lightblue",
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill")

p<-ggviolin(data = top11,
         x = "disease",
         y = "Abundance",
         color = "disease",
         add="point",
#         facet.by = "SampleType",
         ggtheme = theme_bw())+
  scale_color_aaas()+
  facet_grid(SampleType~Genus)+
  yscale(.scale = "log10",.format = T)+rotate_x_text()


p+stat_compare_means(comparisons = my_comparisons,label.y.npc = 10, hide.ns = T)
p+stat_compare_means(hide.ns = T,comparisons = my_comparisons)


ggballoonplot(data = top11,
              size.range = c(3,12),
              x = "Abundance",
              y = "disease",
              size = "Abundance",
              color = "disease",
              fill="Abundance",
              palette="aaas",
              position=position_dodge(),
              ggtheme = theme_bw())+
  scale_fill_continuous(type = "viridis")+xscale("log10",.format = T)+
  facet_grid(facets =SampleType~Genus)+stat_compare_means()

ggballoonplot(data = top1,
              size.range = c(3,12),
              x = "disease",
              y = "avg",
              size = "avg",
              color = "disease",
              fill="Genus",
              palette="Paired",
              position=position_dodge(),
              ggtheme = theme_bw())+yscale("log10",.format = T)

top1
top1$SampleType<-factor(x = top1$SampleType, levels = c(" Introitus","Vagina","Cervix","Anus","Stool","Ear"), )
ggline(data = top1,
       linetype = 3,
       group = "Genus",
       point.size = 2,
       point.color = "Genus",
       numeric.x.axis = T,
       plot_type = "l",
       x = "avg",
       y = "disease",
       palette="Paired",
       add = "point", 
       add.params = list(size="avg",color="Genus"),
       ggtheme = theme_classic2())+
  facet_grid(facets =SampleType~Genus)+
  xscale("log10", .format = T)

library(microbiome) # Load libraries
library(phyloseq)
library(dplyr)
library(reshape2)
library(knitr)
# Z transformed abundance data
taxa_names(ASV_physeq_core_genus)<-get_taxa_unique(ASV_physeq_core_genus,taxonomic.rank = "Genus")





pseqz <-microbiome::transform(ASV_physeq_core_genus_deseq2, transform = "hellinger")


# Plot the abundances heatmap
# Plot the abundances heatmap
dfm <- melt(abundances(pseqz))
dfm <- psmelt(pseqz)
heatmap(otu_table(ASV_physeq_core_genus_rel_deseq2), )
require(graphics); require(grDevices)
x  <- as.matrix(mtcars)
rc <- rainbow(nrow(x), start = 0, end = .3)
cc <- rainbow(ncol(x), start = 0, end = .3)
hv <- heatmap(x, col = cm.colors(256), scale = "column",
              RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
              xlab = "specification variables", ylab =  "Car Models",
              main = "heatmap(<Mtcars data>, ..., scale = \"column\")")
utils::str(hv) # the two re-ordering index vectors

## no column dendrogram (nor reordering) at all:
heatmap(x, Colv = NA, col = cm.colors(256), scale = "column",
        RowSideColors = rc, margins = c(5,10),
        xlab = "specification variables", ylab =  "Car Models",
        main = "heatmap(<Mtcars data>, ..., scale = \"column\")")

## "no nothing"
heatmap(x, Rowv = NA, Colv = NA, scale = "column",
        main = "heatmap(*, NA, NA) ~= image(t(x))")



  heatmap(a )
a<-abundances(ASV_physeq_core_genus_rel_deseq2)
a<-as.matrix(a)
disease<-sample_data(ASV_physeq_core_genus_rel_deseq2)$disease
disease
corrplot(corr = a, )
heatmap(Ca,               symm = TRUE, margins = c(6,6)) # with reorder()

heatmap(Ca, Rowv = FALSE, symm = TRUE, margins = c(6,6)) # _NO_ reorder()
## slightly artificial with color bar, without and with ordering:
cc <- rainbow(nrow(Ca))
cc
heatmap(Ca, Rowv = FALSE, symm = TRUE, RowSideColors = cc, ColSideColors = cc,
        margins = c(6,6))
heatmap(Ca,		symm = TRUE, RowSideColors = cc, ColSideColors = cc,
        margins = c(6,6))

## For variable clustering, rather use distance based on cor():
symnum( cU <- cor(USJudgeRatings) )

hU <- heatmap(cU, Rowv = FALSE, symm = TRUE, col = topo.colors(16),
              distfun = function(c) as.dist(1 - c), keep.dendro = TRUE)
## The Correlation matrix with same reordering:
round(100 * cU[hU[[1]], hU[[2]]])
## The column dendrogram:
utils::str(hU$Colv)
[Package stats version 
  