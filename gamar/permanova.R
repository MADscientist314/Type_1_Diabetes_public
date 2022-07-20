# February 13 STAT

ASV_physeq_trim<-prune_samples(sample_sums(ASV_physeq)>0, ASV_physeq)
ASV_physeq_trim
ASV_physeq_trim<- format_to_besthit(ASV_physeq_trim)
ASV_physeq_rel <-  microbiome::transform(ASV_physeq_trim, "compositional")
otu <- abundances(ASV_physeq_rel)
meta <- meta(ASV_physeq_rel)



library(vegan)
#> Loading required package: permute
#> Loading required package: lattice
#> This is vegan 2.5-3
#> 
#> Attaching package: 'vegan'
#> The following object is masked from 'package:microbiome':
#> 
#>     diversity
permanova <- adonis(t(otu) ~ disease_SampleType,
                    data = meta, permutations = 999, method = "bray")# P-value
permanova

#ear
ear<-subset_samples(ASV_physeq_trim,SampleType=="Ear")
ear
ear_rel<-transform(ear,transform = "compositional")
otu <- abundances(ear_rel)
meta <- meta(ear_rel)
ear_permanova <- adonis(t(otu) ~ disease_SampleType,
                    data = meta, permutations = 999, method = "bray")# P-value
ear_permanova

#cervix
Cervix<-subset_samples(ASV_physeq_trim,SampleType=="Cervix")
Cervix
Cervix_rel<-transform(Cervix,transform = "compositional")
otu <- abundances(Cervix_rel)
meta <- meta(Cervix_rel)
Cervix_permanova <- adonis(t(otu) ~ disease_SampleType,
                        data = meta, permutations = 999, method = "bray")# P-value
Cervix_permanova

#stool
Stool<-subset_samples(ASV_physeq_trim,SampleType=="Stool")
Stool
Stool_rel<-transform(Stool,transform = "compositional")
otu <- abundances(Stool_rel)
meta <- meta(Stool_rel)
Stool_permanova <- adonis(t(otu) ~ disease_SampleType,
                           data = meta, permutations = 999, method = "bray")# P-value
Stool_permanova

#vagina
Vagina<-subset_samples(ASV_physeq_trim,SampleType=="Vagina")
Vagina
Vagina_rel<-transform(Vagina,transform = "compositional")
otu <- abundances(Vagina_rel)
meta <- meta(Vagina_rel)
Vagina_permanova <- adonis(t(otu) ~ disease_SampleType,
                          data = meta, permutations = 999, method = "bray")# P-value
Vagina_permanova

#Anus
Anus<-subset_samples(ASV_physeq_trim,SampleType=="Anus")
Anus
Anus_rel<-transform(Anus,transform = "compositional")
otu <- abundances(Anus_rel)
meta <- meta(Anus_rel)
Anus_permanova <- adonis(t(otu) ~ disease_SampleType,
                           data = meta, permutations = 999, method = "bray")# P-value
Anus_permanova

#introitus = vestibule
Vestibule<-subset_samples(ASV_physeq_trim,SampleType=="Anus")
Vestibule
Vestibule_rel<-transform(Vestibule,transform = "compositional")
otu <- abundances(Vestibule_rel)
meta <- meta(Vestibule_rel)
Vestibule_permanova <- adonis(t(otu) ~ disease_SampleType,
                         data = meta, permutations = 999, method = "bray")# P-value
Vestibule_permanova



