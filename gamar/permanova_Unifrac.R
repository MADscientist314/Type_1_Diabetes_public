tmp<-ASV_physeq_core
stopImplicitCluster()
library(snow)
library(microbiome)
library(doParallel)
library(tidyverse)
library(vegan)
library(broom)
library(cowplot)

cl <- makeCluster(rep("localhost",20), type = "SOCK")
registerDoParallel(cl)
run<-function(S)
  {
  tmp<-ASV_physeq_core
  vst_tmp<-vst_physeq
  
  manifest<-sam%>%mutate(keep=SampleType==S)
  sample_data(tmp)<-sample_data(manifest)
  sample_data(vst_tmp)<-sample_data(manifest)
  
  tmp<-subset_samples(physeq = tmp, keep=="TRUE")
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
  ##########################################################################
  wunifrac<-UniFrac(physeq = tmp,
                    normalized = T,
                    parallel = T,
                    weighted = T)
  vst_wunifrac<-UniFrac(physeq = vst_tmp,
                        normalized = T,
                        parallel = T,
                        weighted = T)
  unifrac<-UniFrac(physeq = tmp,
                   parallel = T,
                   weighted = F)
  vst_unifrac<-UniFrac(physeq = vst_tmp,
                       parallel = T,
                       weighted = F)
  sam<-meta(tmp)
 ##############################################################################
  a<-tidy(adonis2(wunifrac~disease,data = sam,
               permutations = 9999,
               parallel = cl))%>%mutate(SampleType=unique(sam$SampleType),variable="disease",analysis="Weighted Unifrac")
  b<-tidy(adonis2(vst_wunifrac~disease,data = sam,
               permutations = 9999,
               parallel = cl))%>%mutate(SampleType=unique(sam$SampleType),variable="disease",analysis="VST Weighted Unifrac")
  c<-tidy(adonis2(unifrac~disease,data = sam,
               permutations = 9999,
               parallel = cl))%>%mutate(SampleType=unique(sam$SampleType),variable="disease",analysis="Unweighted Unifrac")
  d<-tidy(adonis2(vst_unifrac~disease,data = sam,
               permutations = 9999,
               parallel = cl))%>%mutate(SampleType=unique(sam$SampleType),variable="disease",analysis="VST Unweighted Unifrac")
  #######################################
  e<-tidy(adonis2(wunifrac~Delivery,data = sam,
               permutations = 9999,
               parallel = cl))%>%mutate(SampleType=unique(sam$SampleType),variable="Delivery",analysis="Weighted Unifrac")
  f<-tidy(adonis2(vst_wunifrac~Delivery,data = sam,
               permutations = 9999,
               parallel = cl))%>%mutate(SampleType=unique(sam$SampleType),variable="Delivery",analysis="VST Weighted Unifrac")
  g<-tidy(adonis2(unifrac~Delivery,data = sam,
               permutations = 9999,
               parallel = cl))%>%mutate(SampleType=unique(sam$SampleType),variable="Delivery",analysis="Unweighted Unifrac")
  h<-tidy(adonis2(vst_unifrac~Delivery,data = sam,
               permutations = 9999,
               parallel = cl))%>%mutate(SampleType=unique(sam$SampleType),variable="Delivery",analysis="VST Unweighted Unifrac")
  #######################################
  i<-tidy(adonis2(wunifrac~HbA1C_1trym_2,data = sam,
               permutations = 9999,
               parallel = cl))%>%mutate(SampleType=unique(sam$SampleType),variable="HbA1C_1trym_2",analysis="Weighted Unifrac")
  j<-tidy(adonis2(vst_wunifrac~HbA1C_1trym_2,data = sam,
               permutations = 9999,
               parallel = cl))%>%
    mutate(SampleType=unique(sam$SampleType),variable="HbA1C_1trym_2",analysis="VST Weighted Unifrac")
  k<-tidy(adonis2(unifrac~HbA1C_1trym_2,data = sam,
               permutations = 9999,
               parallel = cl))%>%
    mutate(SampleType=unique(sam$SampleType),variable="HbA1C_1trym_2",analysis="Unweighted Unifrac")
  l<-tidy(adonis2(vst_unifrac~HbA1C_1trym_2,data = sam,
               permutations = 9999,
               parallel = cl))%>%
    mutate(SampleType=unique(sam$SampleType),variable="HbA1C_1trym_2",analysis="VST unweighted Unifrac")
  #######################################
  res<-rbind(a,b,c,d,e,f,g,h,i,j,k,l)
  print(data.frame(res))
  ########################################
  ##########################################################################
  tmp_wunifrac_pcoa_ord<-ordinate(physeq = tmp,method = "PCoA",distance = wunifrac)
  vst_tmp_wunifrac_pcoa_ord<-ordinate(physeq = vst_tmp,method = "PCoA",distance = vst_wunifrac)
  tmp_unifrac_pcoa_ord<-ordinate(physeq = tmp,method = "PCoA",distance = unifrac)
  vst_tmp_unifrac_pcoa_ord<-ordinate(physeq = vst_tmp,method = "PCoA",distance = vst_unifrac)
  eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples
  
  
  tmp_wunifrac_pcoa_ord<-ordinate(physeq = tmp,method = "PCoA",distance = wunifrac)
  vst_tmp_wunifrac_pcoa_ord<-ordinate(physeq = vst_tmp,method = "PCoA",distance = vst_wunifrac)
  tmp_unifrac_pcoa_ord<-ordinate(physeq = tmp,method = "PCoA",distance = unifrac)
  vst_tmp_unifrac_pcoa_ord<-ordinate(physeq = vst_tmp,method = "PCoA",distance = vst_unifrac)
  
  tmp_wunifrac_mds_ord<-ordinate(physeq = tmp,method = "mds",distance = wunifrac)
  vst_tmp_wunifrac_mds_ord<-ordinate(physeq = vst_tmp,method = "mds",distance = vst_wunifrac)
  tmp_unifrac_mds_ord<-ordinate(physeq = tmp,method = "mds",distance = unifrac)
  vst_tmp_unifrac_mds_ord<-ordinate(physeq = vst_tmp,method = "mds",distance = vst_unifrac)
  
  tmp_wunifrac_NMDS_ord<-ordinate(physeq = tmp,method = "NMDS",distance = wunifrac)
  vst_tmp_wunifrac_NMDS_ord<-ordinate(physeq = vst_tmp,method = "NMDS",distance = vst_wunifrac)
  tmp_unifrac_NMDS_ord<-ordinate(physeq = tmp,method = "NMDS",distance = unifrac)
  vst_tmp_unifrac_NMDS_ord<-ordinate(physeq = vst_tmp,method = "NMDS",distance = vst_unifrac)
  #######################################
  # and making our new ordination of just basalts with our adonis statistic
  plot_ordination(vst_tmp, tmp_wunifrac_pcoa_ord, color="disease") + geom_point(size=1) + 
    annotate("text", x=25, y=62, label="Permutational ANOVA = 0.003") + 
    coord_fixed(sqrt(tmp_wunifrac_pcoa_ord$values$Eigenvalues[2]/tmp_wunifrac_pcoa_ord$values$Eigenvalues[1])) + ggtitle("PCoA - basalts only") + 
    scale_color_manual(values=my_pal)+
    theme_bw() + theme(legend.position="right")
  }
glm()
#######################################
#######################################
unique(sam$SampleType)
levels(sam$SampleType)
lapply(levels(sam$SampleType), function(X)run(S =  X))
lapply(levels(sam$SampleType), function(X)run2(S =  X))

#######################################
#######################################

meths<-c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
meths<-c("DCA","NMDS", "MDS", "PCoA")

sam<-sam%>%mutate(controlled=as_factor(HbA1C_1trym_2<7.2))
ordme<-function(SAMTYPE,UNI,ORD){
  print(paste0(SAMTYPE,", Weighted = ",UNI,", Ordination= ",ORD))
  tmp<-ASV_physeq_core
  manifest<-sam%>%mutate(keep=SampleType==SAMTYPE)
  sample_data(tmp)<-sample_data(manifest)
  tmp<-subset_samples(physeq = tmp, keep=="TRUE")
  
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
  unifrac<-UniFrac(physeq = tmp,
                   normalized = T,
                   parallel = T,
                   weighted = UNI)
  vst_unifrac<-UniFrac(physeq = vst_tmp,
                       normalized = F,
                       parallel = T,
                       weighted = UNI)
  
  tmp_ord<-ordinate(physeq = tmp,method = ORD,distance = unifrac)
  vst_tmp_ord<-ordinate(physeq = vst_tmp,method = ORD,distance = vst_unifrac)
  p<-plot_ordination(tmp, tmp_ord,justDF = T)
  #colnames(p)<-c("Axis.1","Axis.2","SampleID","patient","disease","Host","SampleType","HbA1C_1trym_2","Delivery","controlled","keep")
  vp<-plot_ordination(vst_tmp, vst_tmp_ord,justDF = T)
  #colnames(vp)<-c("Axis.1","Axis.2","SampleID","patient","disease","Host","SampleType","HbA1C_1trym_2","Delivery","controlled","keep")
  
  disease_p<-ggplot(p,mapping = aes(x = Axis.1,y=Axis.2,color=disease,shape =disease,group=disease))+
    geom_point(size=1) +
    #stat_ellipse(type = "norm")+
    geom_smooth(se = F,method = "glm")+
    coord_fixed(sqrt(tmp_ord$values$Eigenvalues[2]/tmp_ord$values$Eigenvalues[1])) +
    ggtitle(paste0("Diabetes"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
    scale_color_manual(values=rev(my_pal))+
    theme_cowplot() + theme(legend.position="top")
  Delivery_p<-ggplot(p,mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =Delivery,group=Delivery))+
    geom_point(size=1) +
    #stat_ellipse(type = "norm")+
    geom_smooth(se = F,method = "glm")+
    coord_fixed(sqrt(tmp_ord$values$Eigenvalues[2]/tmp_ord$values$Eigenvalues[1])) +
    ggtitle(paste0("Delivery Mode"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
    scale_color_manual(values=my_pal)+
    theme_cowplot() + theme(legend.position="top")
  controlled_p<-ggplot(p,mapping = aes(x = Axis.1,y=Axis.2,color=controlled,shape =controlled,group=controlled))+
    geom_point(size=1) +
    geom_smooth(se = F,method = "glm")+
    coord_fixed(sqrt(tmp_ord$values$Eigenvalues[2]/tmp_ord$values$Eigenvalues[1])) +
    ggtitle(paste0("1st Trimester HbA1c > 7.2"))+
    scale_color_manual(values=my_pal)+
    theme_cowplot() + theme(legend.position="top")
  disease_vp<-ggplot(vp,mapping = aes(x = Axis.1,y=Axis.2,color=disease,shape =disease,group=disease))+
    geom_point(size=1) +
    geom_smooth(se = F,method = "glm")+
    coord_fixed(sqrt(vst_tmp_ord$values$Eigenvalues[2]/vst_tmp_ord$values$Eigenvalues[1])) +
    ggtitle(paste0("Diabetes"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
    scale_color_manual(values=rev(my_pal))+
    theme_cowplot() + theme(legend.position="top")
  Delivery_vp<-ggplot(vp,mapping = aes(x = Axis.1,y=Axis.2,color=Delivery,shape =Delivery,group=Delivery))+
    geom_point(size=1) +
    geom_smooth(se = F,method = "glm")+
    coord_fixed(sqrt(vst_tmp_ord$values$Eigenvalues[2]/vst_tmp_ord$values$Eigenvalues[1])) +
    ggtitle(paste0("Delivery Mode"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
    scale_color_manual(values=my_pal)+
    theme_cowplot() + theme(legend.position="top")
  controlled_vp<-ggplot(vp,mapping = aes(x = Axis.1,y=Axis.2,color=controlled,shape =controlled,group=controlled))+
    geom_point(size=1) +
    coord_fixed(sqrt(vst_tmp_ord$values$Eigenvalues[2]/vst_tmp_ord$values$Eigenvalues[1])) +
    ggtitle(paste0("1st Trimester HbA1c > 7.2"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
    geom_smooth(se = F,method = "glm")+
    scale_color_manual(values=my_pal)+
    theme_cowplot() + theme(legend.position="top")
  p2<-cowplot::plot_grid(disease_vp,Delivery_vp,controlled_vp,
                         nrow = 1,
                         labels=c(paste0(SAMTYPE," ",ORD," Weighted UniFrac = ",UNI)))
  # p2<-cowplot::plot_grid(disease_p,Delivery_p,controlled_p,disease_vp,Delivery_vp,controlled_vp,
  #                        nrow = 2,
  #                        label_y = c(0,1),
  #                        label_x = c(1,0),
  #                        labels=c(paste0("VST ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI),paste0("Normalized ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI)))#,
  print(p2)
  ggsave(plot = p2,filename = paste0(SAMTYPE,"_",ORD,"_Weighted_",UNI,".png"),path = "cowplots",device = "png",height = 4,width=12,units = "in",dpi = 600)
}

# meths<-c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
# meths<-c("DCA","NMDS", "MDS", "PCoA")

lapply(X=levels(sam$SampleType),FUN =  function(X){
lapply(X = c(TRUE,FALSE),FUN =  function(Y){
lapply(X = c("PCoA"),FUN = function(Z){
#lapply(X = c("disease","Delivery","controlled"),FUN = function(A){
ordme(SAMTYPE = X,UNI = Y,ORD = Z)})})})

