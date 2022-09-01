library(parallel)
library(doParallel)
library(microbiome)
library(tidyverse)
library(ggpubr)
library(cowplot)
library(ggsci)
library(scales)
#tmp<-subset_samples(ASV_physeq_core,Host=="Infant")
sam<-data.frame(sample_data(ASV_physeq))%>%mutate(controlled=as_factor(HbA1C_1trym_2<7.2))
sample_data(tmp)<-sample_data(sam)

sam

sam$controlled
# ###### VST stuff ########
# Run DESEQ2 VST normalaiztion
deseq_counts<- phyloseq_to_deseq2(physeq = tmp,design = ~1)
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)
# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(tmp)
vst_tmp <- phyloseq(vst_count_phy, tax_table(tmp),sample_info_tab_phy,phy_tree(tmp))
#tmp<-prune_taxa(taxa = taxa_sums(tmp) > 0,x = tmp)
#vst_tmp<-prune_taxa(taxa = taxa_sums(vst_tmp) > 0,x = tmp)
cl<-makePSOCKcluster(detectCores())
registerDoParallel(cl)
wunifrac<-UniFrac(physeq = tmp,
                  normalized = T,
                  parallel = T,
                  weighted = T)
vst_wunifrac<-UniFrac(physeq = vst_tmp,
                      normalized = F,
                      parallel = T,
                      weighted = T)

unifrac<-UniFrac(physeq = tmp,
                 normalized = T,
                 parallel = T,
                 weighted = F)
vst_unifrac<-UniFrac(physeq = vst_tmp,
                     normalized = F,
                     parallel = T,
                     weighted = F)



tmp_wunifrac_ord<-ordinate(physeq = tmp,method = "PCoA",distance = wunifrac)
vst_tmp_wunifrac_ord<-ordinate(physeq = vst_tmp,method = "PCoA",distance = vst_wunifrac)
tmp_unifrac_ord<-ordinate(physeq = tmp,method = "PCoA",distance = unifrac)
vst_tmp_unifrac_ord<-ordinate(physeq = vst_tmp,method = "PCoA",distance = vst_unifrac)


p_wunifrac<-plot_ordination(tmp, tmp_wunifrac_ord,justDF = T)
vp_wunifrac<-plot_ordination(vst_tmp, vst_tmp_wunifrac_ord,justDF = T)
p_unifrac<-plot_ordination(tmp, tmp_unifrac_ord,justDF = T)
vp_unifrac<-plot_ordination(vst_tmp, vst_tmp_unifrac_ord,justDF = T)

my_pal2<-pal_aaas(palette = "default",alpha = 1)(10)
# show_col(my_pal2)
################################################################################
disease_wunifrac_p<-ggplot(p_wunifrac,
                           mapping = aes(x = Axis.1,y=Axis.2,color=disease,shape =disease,group=SampleType))+
  geom_point(size=2) +
  scale_color_manual(values=rev(my_pal))+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_wunifrac_ord$values$Eigenvalues[2]/tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("disease"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
disease_wunifrac_p
ggsave(plot = disease_wunifrac_p,filename = "infant_disease_wunifrac_pcoa.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")

disease_wunifrac_vp<-ggplot(vp_wunifrac,
                            mapping = aes(x = Axis.1,y=Axis.2,color=disease,shape =disease,group=SampleType))+
  geom_point(size=2) +
  scale_color_manual(values=rev(my_pal))+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_wunifrac_ord$values$Eigenvalues[2]/vst_tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("disease"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
disease_wunifrac_vp
ggsave(plot = disease_wunifrac_vp,filename = "vst_infant_disease_wunifrac_pcoa.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
#infant untransformed disease unweighted
disease_unifrac_p<-ggplot(p_unifrac,
                          mapping = aes(x = Axis.1,y=Axis.2,color=disease,shape =disease,group=SampleType))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=rev(my_pal2))+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_unifrac_ord$values$Eigenvalues[2]/tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("disease"," Unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
disease_unifrac_p
ggsave(plot = disease_unifrac_p,filename = "infant_disease_unifrac_pcoa.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
#infant untransformed disease unweighted

disease_unifrac_vp<-ggplot(vp_unifrac,
                           mapping = aes(x = Axis.1,y=Axis.2,color=disease,shape =disease,group=SampleType))+
  geom_point(size=2) +
  scale_color_manual(values=rev(my_pal))+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_unifrac_ord$values$Eigenvalues[2]/vst_tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("disease"," unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
disease_unifrac_vp
ggsave(plot = disease_unifrac_vp,filename = "vst_infant_disease_unifrac_pcoa.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")

################################################################################
Delivery_wunifrac_p<-ggplot(p_wunifrac,
                            mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =Delivery,group=SampleType))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_wunifrac_ord$values$Eigenvalues[2]/tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("Delivery"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
Delivery_wunifrac_p
ggsave(plot = Delivery_wunifrac_p,filename = "infant_Delivery_wunifrac_pcoa.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")

Delivery_wunifrac_vp<-ggplot(vp_wunifrac,
                             mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =Delivery,group=SampleType))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_wunifrac_ord$values$Eigenvalues[2]/vst_tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("Delivery"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
Delivery_wunifrac_vp
ggsave(plot = Delivery_wunifrac_vp,filename = "vst_infant_Delivery_wunifrac_pcoa.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
#infant untransformed Delivery unweighted
Delivery_unifrac_p<-ggplot(p_unifrac,
                           mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =Delivery,group=SampleType))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_unifrac_ord$values$Eigenvalues[2]/tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("Delivery"," Unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
Delivery_unifrac_p
ggsave(plot = Delivery_unifrac_p,filename = "infant_Delivery_unifrac_pcoa.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
#infant untransformed Delivery unweighted

Delivery_unifrac_vp<-ggplot(vp_unifrac,
                            mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =Delivery,group=SampleType))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_unifrac_ord$values$Eigenvalues[2]/vst_tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("Delivery"," unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
Delivery_unifrac_vp
ggsave(plot = Delivery_unifrac_vp,filename = "vst_infant_Delivery_unifrac_pcoa.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
################################################################################
################################################################################
controlled_wunifrac_p<-ggplot(p_wunifrac,
                              mapping = aes(x = Axis.1,y=Axis.2,color=controlled,shape =controlled,group=SampleType))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_wunifrac_ord$values$Eigenvalues[2]/tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("controlled"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
controlled_wunifrac_p
ggsave(plot = controlled_wunifrac_p,filename = "infant_controlled_wunifrac_pcoa.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")

controlled_wunifrac_vp<-ggplot(vp_wunifrac,
                               mapping = aes(x = Axis.1,y=Axis.2,color=controlled,shape =controlled,group=SampleType))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_wunifrac_ord$values$Eigenvalues[2]/vst_tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("controlled"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
controlled_wunifrac_vp
ggsave(plot = controlled_wunifrac_vp,filename = "vst_infant_controlled_wunifrac_pcoa.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
#infant untransformed controlled unweighted
controlled_unifrac_p<-ggplot(p_unifrac,
                             mapping = aes(x = Axis.1,y=Axis.2,color=controlled,shape =controlled,group=SampleType))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_unifrac_ord$values$Eigenvalues[2]/tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("controlled"," Unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
controlled_unifrac_p
ggsave(plot = controlled_unifrac_p,filename = "infant_controlled_unifrac_pcoa.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
#infant untransformed controlled unweighted

controlled_unifrac_vp<-ggplot(vp_unifrac,
                              mapping = aes(x = Axis.1,y=Axis.2,color=controlled,shape =controlled,group=SampleType))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_unifrac_ord$values$Eigenvalues[2]/vst_tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("controlled"," unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
controlled_unifrac_vp
ggsave(plot = controlled_unifrac_vp,filename = "vst_infant_controlled_unifrac_pcoa.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
################################################################################
################################################################################
#changing the color scheme
################################################################################
disease_wunifrac_p<-ggplot(p_wunifrac,
                           mapping = aes(x = Axis.1,y=Axis.2,color=SampleType,shape =disease,group=disease))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal2)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=disease))+
  scale_color_manual(values=rev(my_pal))+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_wunifrac_ord$values$Eigenvalues[2]/tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("disease"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
disease_wunifrac_p
ggsave(plot = disease_wunifrac_p,filename = "infant_disease_wunifrac_pcoa_varcolor.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")

disease_wunifrac_vp<-ggplot(vp_wunifrac,
                            mapping = aes(x = Axis.1,y=Axis.2,color=disease,shape =disease,group=disease))+
  geom_point(size=2) +
  scale_color_manual(values=rev(my_pal))+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_wunifrac_ord$values$Eigenvalues[2]/vst_tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("disease"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
disease_wunifrac_vp
ggsave(plot = disease_wunifrac_vp,filename = "vst_infant_disease_wunifrac_pcoa_varcolor.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
#infant untransformed disease unweighted
disease_unifrac_p<-ggplot(p_unifrac,
                          mapping = aes(x = Axis.1,y=Axis.2,color=disease,shape =disease,group=disease))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=rev(my_pal2))+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_unifrac_ord$values$Eigenvalues[2]/tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("disease"," Unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
disease_unifrac_p
ggsave(plot = disease_unifrac_p,filename = "infant_disease_unifrac_pcoa_varcolor.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
#infant untransformed disease unweighted

disease_unifrac_vp<-ggplot(vp_unifrac,
                           mapping = aes(x = Axis.1,y=Axis.2,color=disease,shape =disease,group=disease))+
  geom_point(size=2) +
  scale_color_manual(values=rev(my_pal))+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_unifrac_ord$values$Eigenvalues[2]/vst_tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("disease"," unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
disease_unifrac_vp
ggsave(plot = disease_unifrac_vp,filename = "vst_infant_disease_unifrac_pcoa_varcolor.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")

################################################################################
Delivery_wunifrac_p<-ggplot(p_wunifrac,
                            mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =disease))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=disease,linetype=Delivery))+
  scale_color_manual(values=my_pal)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_wunifrac_ord$values$Eigenvalues[2]/tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("Delivery"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")+facet_wrap(~SampleType)
Delivery_wunifrac_p
Delivery_wunifrac_p<-ggplot(p_unifrac,
                            mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =disease))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=Delivery,linetype=disease))+
  scale_color_manual(values=my_pal)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_wunifrac_ord$values$Eigenvalues[2]/vst_tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("Delivery"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")+facet_wrap(~SampleType)
Delivery_wunifrac_p




ggsave(plot = Delivery_wunifrac_p,filename = "infant_Delivery_wunifrac_pcoa_varcolor.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")

Delivery_wunifrac_vp<-ggplot(vp_wunifrac,
                             mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =Delivery,group=disease))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_wunifrac_ord$values$Eigenvalues[2]/vst_tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("Delivery"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
Delivery_wunifrac_vp
ggsave(plot = Delivery_wunifrac_vp,filename = "vst_infant_Delivery_wunifrac_pcoa_varcolor.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
#infant untransformed Delivery unweighted
Delivery_unifrac_p<-ggplot(p_unifrac,
                           mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =Delivery,group=disease))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_unifrac_ord$values$Eigenvalues[2]/tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("Delivery"," Unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
Delivery_unifrac_p
ggsave(plot = Delivery_unifrac_p,filename = "infant_Delivery_unifrac_pcoa_varcolor.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
#infant untransformed Delivery unweighted

Delivery_unifrac_vp<-ggplot(vp_unifrac,
                            mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =Delivery,group=disease))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_unifrac_ord$values$Eigenvalues[2]/vst_tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("Delivery"," unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
Delivery_unifrac_vp
ggsave(plot = Delivery_unifrac_vp,filename = "vst_infant_Delivery_unifrac_pcoa_varcolor.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
################################################################################
################################################################################
controlled_wunifrac_p<-ggplot(p_wunifrac,
                              mapping = aes(x = Axis.1,y=Axis.2,color=stool_collection_day,shape =controlled,group=disease))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=controlled))+
  scale_color_manual(values=my_pal)+facet_grid(~SampleType)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_wunifrac_ord$values$Eigenvalues[2]/tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("controlled"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
controlled_wunifrac_p
p_wunifrac
ggsave(plot = controlled_wunifrac_p,filename = "infant_controlled_wunifrac_pcoa_varcolor.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")

controlled_wunifrac_vp<-ggplot(vp_wunifrac,
                               mapping = aes(x = Axis.1,y=Axis.2,color=controlled,shape =controlled,group=disease))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_wunifrac_ord$values$Eigenvalues[2]/vst_tmp_wunifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("controlled"," Weighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
controlled_wunifrac_vp
ggsave(plot = controlled_wunifrac_vp,filename = "vst_infant_controlled_wunifrac_pcoa_varcolor.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
#infant untransformed controlled unweighted
controlled_unifrac_p<-ggplot(p_unifrac,
                             mapping = aes(x = Axis.1,y=Axis.2,color=controlled,shape =controlled,group=disease))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_unifrac_ord$values$Eigenvalues[2]/tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("controlled"," Unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
controlled_unifrac_p
ggsave(plot = controlled_unifrac_p,filename = "infant_controlled_unifrac_pcoa_varcolor.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
#infant untransformed controlled unweighted

controlled_unifrac_vp<-ggplot(vp_unifrac,
                              mapping = aes(x = Axis.1,y=Axis.2,color=controlled,shape =controlled,group=disease))+
  geom_point(size=2) +
  scale_color_manual(values=my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=SampleType))+
  scale_color_manual(values=my_pal2)+
  #geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(vst_tmp_unifrac_ord$values$Eigenvalues[2]/vst_tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("controlled"," unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
controlled_unifrac_vp
ggsave(plot = controlled_unifrac_vp,filename = "vst_infant_controlled_unifrac_pcoa_varcolor.png",device = "png",width = 6,height = 4,units = "in",dpi = 600,path = "../gamar/cowplots/")
################################################################################
################################################################################



# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# Delivery_p<-ggplot(p,mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =Delivery,group=Delivery))+
#   geom_point(size=1) +
#   #stat_ellipse(type = "norm")+
#   geom_smooth(se = F,method = "glm")+
#   coord_fixed(sqrt(tmp_ord$values$Eigenvalues[2]/tmp_ord$values$Eigenvalues[1])) +
#   ggtitle(paste0("Delivery Mode"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
#   scale_color_manual(values=my_pal)+
#   theme_cowplot() + theme(legend.position="top")
# controlled_p<-ggplot(p,mapping = aes(x = Axis.1,y=Axis.2,color=controlled,shape =controlled,group=controlled))+
#   geom_point(size=1) +
#   geom_smooth(se = F,method = "glm")+
#   coord_fixed(sqrt(tmp_ord$values$Eigenvalues[2]/tmp_ord$values$Eigenvalues[1])) +
#   ggtitle(paste0("1st Trimester HbA1c > 7.2"))+
#   scale_color_manual(values=my_pal)+
#   theme_cowplot() + theme(legend.position="top")
# disease_vp<-ggplot(vp,mapping = aes(x = Axis.1,y=Axis.2,color=disease,shape =disease,group=disease))+
#   geom_point(size=1) +
#   geom_smooth(se = F,method = "glm")+
#   coord_fixed(sqrt(vst_tmp_ord$values$Eigenvalues[2]/vst_tmp_ord$values$Eigenvalues[1])) +
#   ggtitle(paste0("Diabetes"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
#   scale_color_manual(values=rev(my_pal))+
#   theme_cowplot() + theme(legend.position="top")
# Delivery_vp<-ggplot(vp,mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =Delivery,group=Delivery))+
#   geom_point(size=1) +
#   geom_smooth(se = F,method = "glm")+
#   coord_fixed(sqrt(vst_tmp_ord$values$Eigenvalues[2]/vst_tmp_ord$values$Eigenvalues[1])) +
#   ggtitle(paste0("Delivery Mode"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
#   scale_color_manual(values=my_pal)+
#   theme_cowplot() + theme(legend.position="top")
# controlled_vp<-ggplot(vp,mapping = aes(x = Axis.1,y=Axis.2,color=controlled,shape =controlled,group=controlled))+
#   geom_point(size=1) +
#   coord_fixed(sqrt(vst_tmp_ord$values$Eigenvalues[2]/vst_tmp_ord$values$Eigenvalues[1])) +
#   ggtitle(paste0("1st Trimester HbA1c > 7.2"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
#   geom_smooth(se = F,method = "glm")+
#   scale_color_manual(values=my_pal)+
#   theme_cowplot() + theme(legend.position="top")
# p2<-cowplot::plot_grid(disease_vp,Delivery_vp,controlled_vp,
#                        nrow = 1,
#                        labels=c(paste0(SAMTYPE," ",ORD," Weighted UniFrac = ",UNI)))
# # p2<-cowplot::plot_grid(disease_p,Delivery_p,controlled_p,disease_vp,Delivery_vp,controlled_vp,
# #                        nrow = 2,
# #                        label_y = c(0,1),
# #                        label_x = c(1,0),
# #                        labels=c(paste0("VST ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI),paste0("Normalized ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI)))#,
# print(p2)
# ggsave(plot = p2,filename = paste0(SAMTYPE,"_",ORD,"_Weighted_",UNI,".png"),path = "cowplots",device = "png",height = 4,width=12,units = "in",dpi = 600)
# }
# 
# # meths<-c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
# # meths<-c("DCA","NMDS", "MDS", "PCoA")
# 
# lapply(X=levels(sam$SampleType),FUN =  function(X){
# lapply(X = c(TRUE,FALSE),FUN =  function(Y){
#   lapply(X = c("PCoA"),FUN = function(Z){
#     #lapply(X = c("disease","Delivery","controlled"),FUN = function(A){
#     ordme(SAMTYPE = X,UNI = Y,ORD = Z)})})})