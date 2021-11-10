library(microbiomeutilities)
library(phyloseq)
sample_variables(ASV_physeq_core)

ASV_physeq_core_ord<-ordinate(physeq = Infant_physeq_core, method = 'NMDS', distance = "bray")
Infant_physeq_core<-subset_samples(ASV_physeq_core, Host=="Infant")
Mother_physeq_core<-subset_samples(ASV_physeq_core, Host=="Mother")
Mother_physeq_core<-subset_samples(Mother_physeq_core, SampleType!="Anus")

microbiomeutilities::plot_taxa_composition(x = Infant_physeq_core, 
                                          sample.sort = 'neatmap', 
                                           taxonomic.level = "Phylum",
                                           otu.sort = "abundance",
                                           transform = "compositional",
                                           plot.type = 'barplot',
                                           average_by = "Antibiotics"
                                            )
microbiomeutilities::plot_taxa_composition(x = Infant_physeq_core, 
                                           sample.sort = 'neatmap', 
                                           taxonomic.level = "Phylum",
                                           otu.sort = "abundance",
                                           transform = "compositional",
                                           plot.type = 'barplot',
                                           average_by = "disease_SampleType"
)
ASV_physeq_core_ord<-ordinate(physeq = Infant_physeq_core, method = 'NMDS', distance = "bray")
Infant_physeq_core<-subset_samples(ASV_physeq_core, Host=="Infant")
Mother_physeq_core<-subset_samples(ASV_physeq_core, Host=="Mother")
Mother_physeq_core<-subset_samples(Mother_physeq_core, SampleType!="Anus")


sample_names(Infant_physeq_core)
meta<-as.data.frame(sample_data(infant_physeq_core))
meta$RODZAJ_PORODU
shapes=c(8,9,10,11,2,13,14,15,16,17,18)
Mother_physeq_core_ord<-ordinate(physeq = Mother_physeq_core, method = 'NMDS', distance = "bray")
shapes=c(8,9,10,11,2,13,14,15,16,17,18)
plot_ordination(Mother_physeq_core, Mother_physeq_core_ord, type = "samples", axes = 1:2,
                color = "disease_SampleType", shape = "SGA_AGA_LGA", label = NULL, title = NULL,
                justDF = FALSE)+scale_shape_manual(values=shapes)+theme_classic()
)
microbiomeutilities::plot_taxa_composition(x = Mother_physeq_core, 
                                           sample.sort = 'neatmap', 
                                           taxonomic.level = "Phylum",
                                           otu.sort = "abundance",
                                           transform = "compositional",
                                           plot.type = 'barplot',
                                           average_by = "SGA_AGA_LGA")

library(vegan)
dist<-vegdist()
rodzaj.aov<-aov