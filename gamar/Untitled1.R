# To reproduce.
set.seed(1024)
# generate a tree contained 100 tip labels.
formatto

ASV_physeq_core<-subset_taxa(physeq = ASV_physeq_core,Kingdom=="Bacteria")
ASV_physeq_core

ASV_physeq_phylum<-tax_glom(physeq = ASV_physeq_core,taxrank = "Genus")


taxa_names(ASV_physeq_phylum)<-get_taxa_unique(physeq = ASV_physeq_phylum,taxonomic.rank = "Genus")

ASV_physeq_phylum_comp<-transform(ASV_physeq_phylum,transform = "compositional")

sam<-meta(ASV_physeq_phylum_comp)
sam
sample_data(ASV_physeq_phylum_comp)$rep<-paste0(sample_data(ASV_physeq_phylum_comp)$disease,"_",
                                           sample_data(ASV_physeq_phylum_comp)$Host,"_",
                                           sample_data(ASV_physeq_phylum_comp)$SampleType,"_",
                                           sample_data(ASV_physeq_phylum_comp)$Delivery)
ASV_physeq_phylum_comp_merged<-merge_samples(x =  ASV_physeq_phylum_comp,
                                       group = paste0(sample_data(ASV_physeq_phylum_comp)$disease,"_",
                                        sample_data(ASV_physeq_phylum_comp)$Host,"_",
                                        sample_data(ASV_physeq_phylum_comp)$SampleType,"_",
                                        sample_data(ASV_physeq_phylum_comp)$Delivery),
                                        fun = mean)

microbiome::bes
library(microbiome)
data(atlas1006)

p0.f <- add_besthit(atlas1006, sep=":")

dietswap2<-aggregate_taxa(x = ASV_physeq,level = "Genus",verbose = T)
dietswap2
dietswap2<-aggregate_taxa(x = dietswap,level = "Genus",verbose = T)
test<-ASV_physeq

d2<-data.frame(tax_table(test))%>%mutate(Domain="Bacteria",
                                         Species=rownames(d2))%>%
  select(Domain,Kingdom,Phylum,Class,Order,Family,Genus,Species)
d2<-tax_table(object = as.matrix(d2))
tax_table(test)<-d2
remove_taxa()
p0.f <- add_besthit(test, sep=":")
taxa_names(p0.f)




ASV_physeq_phylum_comp_merged_small<-core(x = ASV_physeq_phylum_comp_merged,detection = 1/100,prevalence = 40/100)
ASV_physeq_phylum_comp_merged_small
ASV_physeq_phylum_comp_merged_small<-  prune_taxa(x = ASV_physeq_phylum_comp_merged,
                                                  taxa =  top_taxa(x = ASV_physeq_phylum_comp_merged,
                                                                   n = 20))
                                           
ASV_physeq_phylum_comp_merged_small
sam<-meta(ASV_physeq_phylum_comp_merged_small)%>%select(-rep)
sam<-sam%>%mutate(SampleID=rownames(sam))#%>%select(-rep)
sam
sam<-sam%>%separate(SampleID,into = c("disease","Host","SampleType","Delivery"),sep = "_")
sam
sample_data(ASV_physeq_phylum_comp_merged_small)<-sample_data(sam)
sample_data(ASV_physeq_phylum_comp_merged_small)

tr <- phy_tree(ASV_physeq_phylum_comp_merged_small)
tr$tip.label

library(treeio)
# generate three datasets, which are the same except the third column name.
# To reproduce.
set.seed(1024)
# generate three datasets, which are the same except the third column name.
library(ggtree)
tree <- read.nhx(textConnection(treetext))
df <-psmelt(ASV_physeq_phylum_comp_merged_small)%>%
  mutate(id=OTU,value=Abundance,group=disease)#%>%filter(Abundance>0.5)
df2<-df%>%group_by(OTU,Genus,disease)%>%summarise(Abundance=mean(Abundance))%>%ungroup()
df2
ASV_physeq_phylum_merged_comp
layouts<-c('rectangular', 'dendrogram', 'slanted', 'ellipse', 'roundrect', 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle', 'daylight' ,'ape')
# plot tree
# #lapply(layouts, function(X) ggtree(tr, layout=X, open.angle=0))

tr
p<-ggtree(tr) + geom_tiplab() #+ 
#  geom_text(aes(label=df2$Genus), hjust=-.5)
p
p <- ggtree(tr, layout="fan", open.angle=0,)+ geom_tiplab(offset = 0.1) 
p
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.

my_pal<-pal_aaas()(10)

# the first ring.

p1 <- p +
  geom_fruit(
    data = df,
    geom = geom_bar,
    mapping = aes(y=OTU, x=Abundance, fill=SampleType),
    orientation = "y",offset = 0.15,
    stat = "identity")+
  scale_fill_aaas()+
  ggnewscale::new_scale_fill()
#
p1
# # the second ring
# # geom_fruit_list is a list, which first element must be layer of geom_fruit.
p2 <- p1 +  geom_fruit(
  data=df,
  geom=geom_bar,
  mapping=aes(y=OTU,
              x=Abundance,
              color=disease,
              fill=disease),
  orientation="y",
  axis.params = list(axis = "y", 
                     text.angle = 0, 
                     text.size = 8, 
                     text = "das;kfj", 
                     title= "dissadfease", 
                     title.size = 30, 
                     title.height = 0.1, 
                     title.angle = 0, 
                     title.color = "black",
                     nbreak = 4, 
                     line.size = 0.2, 
                     line.color = "grey", 
                     line.alpha = 1,),
  # offset = 0.03,
  # position = position_stack(),
  stat="identity") +
  scale_color_manual(values=c(my_pal[10],my_pal[6]))+
  scale_fill_manual(values=c(my_pal[10],my_pal[6]))+
  new_scale_fill()+new_scale_color()




#
p2
df$value
p3 <- p2 + 
  geom_fruit(
    data=df,
    geom=geom_boxplot,
    mapping = aes(
      y=OTU,
      x=value,
      group=label,
      fill=Genus,label=Genus,
    ),
    # size=.2,
    # outlier.size=0.5,
    # outlier.stroke=0.08,
    # outlier.shape=21,
    # axis.params=list(
    #   axis       = "x",
    #   text.size  = 1.8,
    #   hjust      = 1,
    #   vjust      = 0.5,
    #   nbreak     = 3,
    # ),
  #  grid.params=list()
  ) +  scale_fill_d3("category20")+
#  facet_grid(SampleType~disease)+
  ggnewscale::new_scale_fill()

p3


data.frame(melt_simple)
df
p4 <- p3 +
  geom_fruit(
    data=df,
    geom=geom_tboxplot,
    mapping = aes(
      y=OTU,
      x=value,
      group=label,
      fill=Genus,label=Genus,
    ),
    # size=.2,
    # outlier.size=0.5,
    # outlier.stroke=0.08,
    # outlier.shape=21,
    # axis.params=list(
    #   axis       = "x",
    #   text.size  = 1.8,
    #   hjust      = 1,
    #   vjust      = 0.5,
    #   nbreak     = 3,
    # ),
    #  grid.params=list()
  ) +  scale_fill_d3("category20")+
  #  facet_grid(SampleType~disease)+
  ggnewscale::new_scale_fill()


geom_taxalink() 

p4

df
p5 <- p4 +
  scale_fill_discrete(
    name="Genus",
    guide=guide_legend(keywidth=0.8, keyheight=0.8, ncol=1)
  ) +
  theme(
    legend.title=element_text(size=9), 
    legend.text=element_text(size=7) 
  )





# library(ggtreeExtra)
# library(ggtree)
# library(phyloseq)
# library(dplyr)
#
# data("GlobalPatterns")
# GP <- GlobalPatterns
# GP <- prune_taxa(taxa_sums(GP) > 600, GP)
# sample_data(GP)$human <- get_variable(GP, "SampleType") %in%
#   c("Feces", "Skin")
# mergedGP <- merge_samples(GP, "SampleType")
# mergedGP <- rarefy_even_depth(mergedGP,rngseed=394582)
# mergedGP <- tax_glom(mergedGP,"Order")
#
mergedCP
melt_simple <- psmelt(ASV_physeq_phylum_comp_merged_small) %>%
   filter(Abundance > 0.01)%>%group_by(disease,Host,SampleType,Genus)%>%summarise(val=mean(Abundance*100))%>%ungroup()
  # select(OTU, val=Abundance)

melt_simple
p <- ggtree(ASV_physeq_phylum_comp_merged_small, layout="fan", open.angle=10) + 
  geom_tippoint(mapping=aes(color=Genus), 
                size=1.5,
                show.legend=FALSE)+scale_color_d3(palette = "category20")
p
p1 <- rotate_tree(p, -90)
p1
melt_simple$disease
melt_simple$Class
melt_simple


p2 <- p1 +
  geom_fruit(
    data=melt_simple,
    geom=geom_boxplot,
    mapping = aes(
      y=OTU,
      x=val,group=SampleType,
      fill=SampleType)) 

p2

p <- p +
  scale_fill_discrete(
    name="Genus",
    guide=guide_legend(keywidth=0.8, keyheight=0.8, ncol=1)
  ) +
  theme(
    legend.title=element_text(size=9), 
    legend.text=element_text(size=7) 
  )
p

