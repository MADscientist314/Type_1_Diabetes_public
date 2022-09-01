tmp2<-ASV_physeq#subset_samples(physeq = ASV_physeq,Host=="Infant")


tt<-data.frame(sample_data(t))%>%
  select(patient,disease,SampleType,Host,Glu1,newborn_weight,Antibiotics_2)%>%
  filter(Host=="Infant")%>%
  filter(SampleType=="Stool")%>%
  filter(!is.na(Glu1))
tt<-tt%>%mutate(SampleID=rownames(tt))
tt
dim(tt)
dim(t2)
tt
t2<-data.frame(sample_data(tmp2))
t2
t3<-full_join(tt,t2)
t3
t4<-t3%>%filter(SampleID%in%vp_unifrac_outlier$SampleID)%>%mutate(outlier="outlier")
t5<-t3%>%filter(!SampleID%in%vp_unifrac_outlier$SampleID)%>%mutate(outlier="norm")

t6<-full_join(t4,t5)
t6
rownames(t6)<-t6$SampleID
sample_data(tmp2)<-sample_data(t6)

tmp3<-prune_taxa(taxa_sums(tmp2) > 1, tmp2) 

tmp3<-aggregate_taxa(x =tmp3,level = "Genus",verbose = T)
tmp3<-transform(tmp3,"compositional")

tmp3







ggboxplot(data = out,x = "outlier",y = "Abundance",color = "outlier",add = "jitter")+
  facet_wrap(~Genus,scales = "free")+
  stat_compare_means(comparisons = my_comparisons)
unique(out$Genus)
keep<-c("Bacteroides","Bifidobacterium","Clostridium_sensu_stricto_1","Enterococcus","Escherichia/Shigella","Parabacteroides","Staphylococcus","Unknown")

out2<-out%>%filter(Genus%in%keep)
ggviolin(trim = T,
         data = out2,
         x = "outlier",
         y = "Abundance",color="outlier",
         shape = "disease",
         group="outlier",add="jitter")+
  facet_wrap(~Genus,scales = "free")+scale_color_manual(values = rev(my_pal))+
  stat_compare_means(comparisons = my_comparisons)


ggbarplot(data = out2,x = "outlier",y ="Abundance",color =  "Genus")+scale_color_aaas()

plot_composition(x = tmp4,sample.sort = "outlier",verbose = T,average_by = "outlier")+scale_fill_d3()
tmp<-tmp2
tmp2
taxa_sums(tmp)
sample_data(tmp)$Glu1

out2

out
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
tmp<-subset_taxa(physeq = tmp,Genus!="Escherichia/Shigella")
tmp
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
vp_unifrac$Antibiotics_2

controlled_wunifrac_p<-ggplot(p_unifrac,
                              mapping = aes(x = Axis.1,y=Axis.2,color=Antibiotics_2,shape =Delivery,group=disease))+
  geom_point(size=2) +
#  scale_color_viridis_c()+
  scale_color_manual(values = my_pal)+
  ggnewscale::new_scale_color()+
  stat_ellipse(type = "t",mapping = aes(color=disease))+
  scale_color_manual(values = my_pal)+
  
  facet_grid(disease~Delivery)+
  geom_smooth(se = F,method = "glm")+
  coord_fixed(sqrt(tmp_unifrac_ord$values$Eigenvalues[2]/tmp_unifrac_ord$values$Eigenvalues[1])) +
  ggtitle(paste0("controlled"," Unweighted Unifrac"))+#ANN," ",ORD," ",SAMTYPE," Weighted UniFrac = ",UNI," ")) +
  theme_pubr() + theme(legend.position="right")
controlled_wunifrac_p
vp_unifrac_outlier<-vp_unifrac%>%filter(Axis.1<0)
vp_unifrac_outlier
u<-data.frame(sample_data(t))%>%mutate(SampleID=rownames(u))

u$SampleID
vp_unifrac_outlier<-u%>%filter(SampleID%in%vp_unifrac_outlier$SampleID)

vp_unifrac_outlier

tmp2
w

