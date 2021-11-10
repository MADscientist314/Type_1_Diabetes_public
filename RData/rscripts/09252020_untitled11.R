library(microbiome)
plot_composition(x = ASV_physeq_core, 
                 sample.sort = "disease", 
                 otu.sort = "Abundance", 
                 plot.type = 'barplot', 
                 group_by = "SampleType",
                 verbose = T, 
                 average_by = "disease")
library(dplyr)
library(phyloseq)
pseq <- ASV_physeq_core
#pseq<-subset_samples(disease=="Control", physeq = pseq)
pseq<-tax_glom(physeq = pseq,taxrank = "Phylum",NArm = T)
taxa_names(pseq)<-get_taxa_unique(physeq = pseq, taxonomic.rank = "Phylum")
pseq<-transform(x = pseq, transform = "compositional")
plot_composition(pseq, sample.sort = "neatmap", 
                 otu.sort = "abundance", 
                 group_by = "Phylum",x.label = "patient",
                 verbose = T) +scale_fill_brewer(palette = "Paired")
library(shiny)
sample_data(pseq)
#
data(esophagus)
tree <- phy_tree(esophagus)
otu  <- otu_table(esophagus)
otutree0 <- phyloseq(otu, tree)
 plot_tree(otutree0)
otutree1 <- merge_taxa(otutree0, 1:8, 2)
 plot_tree(esophagus, ladderize="left")
 microbiome::
saveRDS(pseq,"pseq.rds")
save.image(file = "pseq.rda") 
#install.packages("shiny")
#shiny::runGitHub("shiny-phyloseq","joey711")
otu<-data.frame(otu_table(pseq))
otu<-as.tibble(otu)
rownames(otu)
a<-psmelt(pseq)
a<-summarise(.data = a, .groups = c(Sampletype,disease))
a$disease_SampleType
plot_bar(physeq = pseq, fill = "Phylum", title = "Type 1 Diabetes vs Control", facet_grid = disease~SampleType)
b<-a %>%
  group_by(disease_SampleType, Abundance, Phylum) %>%
  summarise(OTU, n = n())
b$n
ggplot(data = b,mapping = aes(x = disease_SampleType, y=Abundance))+
  geom_dotplot()

library(ggpubr)
b<-a%>%group_by(disease,SampleType,Phylum)%>%summarise(mean = mean(Abundance))
b
mean(Abundance))
pseq.ord<-ordinate(ASV_physeq_core, method = "NMDS", distance = "bray")
a<-plot_ordination(physeq = ASV_physeq_core, ordination = pseq.ord, type = "ssampled", color = "disease")
a+facet_wrap(facets = ~SampleType)+scale_color_brewer(palette = "Set1")
ggplot(a,mapping = aes(x = NMDS1,y = NMDS2,color=Phylum))+
         geom_point()+facet_grid(facets = disease~Phylum)+scale_color_brewer(palette = "Paired")

a
ggballoonplot(data = b, 
              x = "disease", 
              y = "Phylum", 
              size = "mean", 
              size.range = c(1,20),
              fill = "mean", ggtheme = theme_bw())+
  scale_fill_viridis_c()+
  # scale_fill_brewer(direction = -1,palette = "Paired")+
  facet_grid(facets = ~SampleType)
