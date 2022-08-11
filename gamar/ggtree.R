library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(dplyr)
library(microbiome)

mergedGP <- rarefy_even_depth(ASV_physeq,rngseed=394582)
mergedGP <- tax_glom(mergedGP,"Phylum")
mergedGP
mergedGP <- microbiome::transform(mergedGP,"compositional")

melt_simple <- psmelt(mergedGP) %>%
  filter(Abundance < 120) %>%
  select(OTU, val=Abundance)
melt_simple


m2<-melt_simple%>%
  group_by(patient,disease,Host,Delivery,Phylum)%>%
  summarise(Abundance=sum(Abundance))%>%ungroup()%>%
  rename("Abundance"="val")%>%filter(val>1)

melt_simple
m2
p <- ggtree(mergedGP, layout="fan", open.angle=10) +
  geom_tippoint(mapping=aes(color=Phylum), 
                size=1.5,
                show.legend=FALSE)
p
p <- rotate_tree(p, -90)
melt_simple
p <- p +
  geom_fruit(
    data=melt_simple,
    geom=geom_boxplot,
    mapping = aes(
      y=Phylum,
      x=val,
      #group=disease,
      # fill=disease,
      # label=Phylum,
    ),
    size=.2,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
    ),
    grid.params=list()
  ) 

p





p <- p +
  scale_fill_discrete(
    name="Phyla",
    guide=guide_legend(keywidth=0.8, keyheight=0.8, ncol=1)
  ) +
  theme(
    legend.title=element_text(size=9), 
    legend.text=element_text(size=7) 
  )

p



p3 <- p + 
  new_scale_fill() + 
  geom_fruit(
    data=melt_simple,
    geom=geom_tile,
    mapping=aes(y=ID, x=Pos, fill=Type),
    offset=0.08,   # The distance between external layers, default is 0.03 times of x range of tree.
    pwidth=0.25 # width of the external layer, default is 0.2 times of x range of tree.
  ) + 
  scale_fill_manual(
    values=c("#339933", "#dfac03"),
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=3)
  ) 
p3

