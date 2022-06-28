library(tidyverse)
library(phyloseq)
library(mosaic)
library(ggsci)
library(ggpubr)
library(microbiome)
library(DESEq2)


#set some visualization themes and functions
theme_set(theme_bw())


# Running the DESeq2 analysis
#make a list of the phyloseq objects subsampled by SampleType
list<-c(ASV_physeq_core_Anus,ASV_physeq_core_Cervix,ASV_physeq_core_Ear,ASV_physeq_core_Introitus,ASV_physeq_core_Stool,ASV_physeq_core_Vagina)

#Agglomerate to a genus level
tax_glom<-lapply(list, function(x) tax_glom(physeq = x,taxrank = "Genus"))
tax_glom

# rename the taxa for the genus names
rename<-function(X)
{
  pseq2<-X
  tax<-data.frame(tax_table(pseq2))
  taxa_names(pseq2)<-tax$Genus
  return(pseq2)
}

renamed<-lapply(tax_glom, function(X) rename(X = X))
renamed
#run deseq2 on the renamed stuff
deseq_me<-function(X){
  print(paste0("LOG: running DESEQ2 on ", unique(meta(X)$SampleType)))
  if(is.na(meta(X)$Delivery))
    {
    ds <- phyloseq_to_deseq2(X, ~ disease)
    }
  else
    {
    print(paste0("controlling for delivery"))
    ds <- phyloseq_to_deseq2(X, ~ Delivery+disease)
    }
  ds
  #convert the disease column to factor
  ds$disease <- factor(ds$disease, levels = c("Control","T1D"))
  ds$disease <- relevel(ds$disease, ref = "Control")
  #run deseq2
  dds<-DESeq(ds)
  dds
  # #get the results
  res <- results(dds, contrast=c("disease","T1D","Control"))
  #res = results(dds, cooksCutoff = FALSE)
  # resultsNames(dds)
  print(res)
  # #Log fold change shrinkage for visualization and ranking
  #resLFC <- lfcShrink(dds, coef="disease_T1d_vs_Control")# type="apeglm")
  # resLFC
  
  #Investigate test results table
  alpha = 0.01 #significance level cutoff
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(X)[rownames(sigtab), ], "matrix"))
  print(sigtab)
  return(sigtab)
  #print a MA plot
  print(plotMA(res,
               main = paste0(unique(meta(X)$SampleType), "MA-plot"),
               xlab = "mean of normalized counts", ylim=c(-2,2)))
  #print(plotMA(resLFC, ylim=c(-2,2)))
  #####################
  # #visualiztion#########
  # 
  # #Let's look at the OTUs that were significantly different between the two tissues. 
  # #The following makes a nice ggplot2 summary of the results.
  # 
  # # Phylum order
  # x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  # # Genus order
  # x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
  # ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  #   theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

  return(sigtab)
}

sam

res_list<-lapply(renamed,deseq_me)
#index teh sampletypes
res_list[[1]]$SampleType<-"Anus"
res_list[[2]]$SampleType<-"Cervix"
res_list[[3]]$SampleType<-"Ear"
res_list[[4]]$SampleType<-"Vaginal Introitus"
res_list[[5]]$SampleType<-"Stool"
res_list[[6]]$SampleType<-"Midvaginal"

res2<-data.frame(unlist(res_list[[1]]))
res2<-full_join(res_list[[1]],res_list[[2]])
res2<-full_join(res2,res_list[[3]])
res2<-full_join(res2,res_list[[4]])
res2<-full_join(res2,res_list[[5]])
res2<-full_join(res2,res_list[[6]])

res2<-as_tibble(res2)

res3<-left_join(res2,sam)%>%
  select(baseMean, log2FoldChange,lfcSE,stat,pvalue,padj,Kingdom,Phylum,Class,Order,Family,Genus,SampleType,disease,Delivery)
res3
my_pal<-get_palette(palette = "aaas",k = 35)
res2$SampleType
res3
res
ggscatter(res2,
              x= "SampleType",
              y = "Genus",
              size = "log2FoldChange",
              size.range = c(1,7),
              color = "SampleType",
          title = paste0("Log2 fold Changes amongst Differential Taxa associated with Type 1 Diabetes vs Control"),
              fill="SampleType")+
  scale_color_npg()+
  scale_fill_npg()+
  # scale_color_manual(values = my_pal)+
  # scale_fill_manual(values = my_pal)+
  # scale_color_viridis_c()+
  # scale_fill_viridis_c()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, face = "italic", hjust =1, vjust=1))+
  theme(axis.text.y = element_text(angle = 0, face = "italic", hjust =1, vjust=1))
range(res3$log2FoldChange)
my_pal<-get_palette(palette = "aaas",k = 6)

plots<-lapply(rep(1:6), function(X) ggscatter(data = res_list[[X]],
                                      x= "log2FoldChange",
                                      y = "Genus",
                                      size = "log2FoldChange",
                                      #size.range = c(-20,20),
                                      color = "SampleType",
                                      fill="SampleType")+
                scale_size(limits = c(-20,20))+
                xscale(.scale = "log2")+#xlim(c(-20,20))+
                scale_color_manual(values = my_pal[X])+
                scale_fill_manual(values = my_pal[X])+
                theme_bw()+
                theme(axis.text.x = element_text(angle = 45, face = "italic", hjust =1, vjust=1))+
                theme(axis.text.y = element_text(angle = 0, face = "italic", hjust =1, vjust=1)))
plots
###############################################
###############################################
###############################################
Sys.time()
#"2022-01-21 17:21:04 CST"
# I stopped her trying to figure out cowplots with same scale sizes and colors
###############################################
###############################################
###############################################



library(ggplot2)
library(cowplot)
anus<-ggscatter(data = res_list[[1]],
                x= "log2FoldChange",
                y = "Genus",
                size = "log2FoldChange",
                xlab = "",
                ylab = "",
                #size.range = c(-25,25),
                color = "SampleType",
                fill="SampleType")+
  scale_size(limits = c(-25,25))+
  scale_color_manual(values = my_pal[1])+
  scale_fill_manual(values = my_pal[1])+
  theme_bw(base_size = 25)+
  xlim(c(-25,25))+
  theme(axis.text.x = element_text(angle = 45, face = "italic", hjust =0, vjust=0))+
  theme(axis.text.y = element_text(angle = 0, face = "italic", hjust =0, vjust=0))
anus<-ggpar(p = anus,legend = "none")

cervix<-ggscatter(data = res_list[[2]],
                  x= "log2FoldChange",
                  y = "Genus",
                  size = "log2FoldChange",
                  #size.range = c(-25,25),
                  xlab = "",
                  ylab = "",
                  color = "SampleType",
                  fill="SampleType")+
  scale_size(limits = c(-25,25))+
  scale_color_manual(values = my_pal[2])+
  scale_fill_manual(values = my_pal[2])+
  xlim(c(-25,25))+
  theme_bw(base_size = 25)+
  theme(axis.text.x = element_text(angle = 45, face = "italic", hjust =0, vjust=0))+
  theme(axis.text.y = element_text(angle = 0, face = "italic", hjust =0, vjust=0))
cervix<-ggpar(p = cervix,legend = "none")

ear<-ggscatter(data = res_list[[3]],
                      x= "log2FoldChange",
                      y = "Genus",
                      size = "log2FoldChange",
                      #size.range = c(-25,25),
                      xlab = "",
                      ylab = "",
                      color = "SampleType",
                      fill="SampleType")+
  scale_size(limits = c(-25,25))+
  xlim(c(-25,25))+
  scale_color_manual(values = my_pal[3])+
  scale_fill_manual(values = my_pal[3])+
  theme_bw(base_size = 25)+
  theme(axis.text.x = element_text(angle = 45, face = "italic", hjust =0, vjust=0))+
  theme(axis.text.y = element_text(angle = 0, face = "italic", hjust =0, vjust=0))
ear<-ggpar(p = ear,legend = "none")

intro<-ggscatter(data = res_list[[4]],
                      x= "log2FoldChange",
                      y = "Genus",
                      size = "log2FoldChange",
                      #size.range = c(-25,25),
                      xlab = "",
                      ylab = "",
                      color = "SampleType",
                      fill="SampleType")+
  scale_size(limits = c(-25,25))+
  xlim(c(-25,25))+
  scale_color_manual(values = my_pal[4])+
  scale_fill_manual(values = my_pal[4])+
  theme_bw(base_size = 25)+
  theme(axis.text.x = element_text(angle = 45, face = "italic", hjust =0, vjust=0))+
  theme(axis.text.y = element_text(angle = 0, face = "italic", hjust =0, vjust=0))
intro<-ggpar(p = intro,legend = "none")
stool<-ggscatter(data = res_list[[5]],
                      x= "log2FoldChange",
                      y = "Genus",
                      size = "log2FoldChange",
                      #size.range = c(-25,25),
                      xlab = "",
                      ylab = "",
                      color = "SampleType",
                      fill="SampleType")+
  scale_size(limits = c(-25,25))+
  xlim(c(-25,25))+
  scale_color_manual(values = my_pal[5])+
  scale_fill_manual(values = my_pal[5])+
  theme_bw(base_size = 25)+
  theme(axis.text.x = element_text(angle = 45, face = "italic", hjust =0, vjust=0))+
  theme(axis.text.y = element_text(angle = 0, face = "italic", hjust =0, vjust=0))
stool<-ggpar(p = stool,legend = "none")


midvaginal<-ggscatter(data = res_list[[6]],
                  x= "log2FoldChange",
                  y = "Genus",
                  size = "log2FoldChange",
                  #size.range = c(-25,25),
                  xlab = "",
                  ylab = "",
                  color = "SampleType",
                  fill="SampleType")+
  scale_size(limits = c(-25,25))+
  xlim(c(-25,25))+
  scale_color_manual(values = my_pal[6])+
  scale_fill_manual(values = my_pal[6])+
  theme_bw(base_size = 25)+
  theme(axis.text.x = element_text(angle = 45, face = "italic", hjust =0, vjust=0))+
  theme(axis.text.y = element_text(angle = 0, face = "italic", hjust =0, vjust=0))
midvaginal<-ggpar(p = midvaginal,legend = "none")


plot_grid(intro,midvaginal,cervix,anus,ear,stool, labels = c('A)', 'B)','C)','D)','E)','F)'))


plot_grid(intro,midvaginal,cervix,anus, labels = c('A)', 'B)','C)','D)'),label_size = 24)

plot_grid(ear,stool, labels = c('E)','F)'),label_size = 24)


res_list[[6]]









## S4 method for signature 'DESeqDataSet'
plotMA(object, 
       alpha = 0.1, 
       main = paste0(unique(meta(X)$SampleType), "MA-plot")m
       xlab = "mean of normalized counts", 
       ylim, 
       MLE = FALSE, )

## S4 method for signature 'DESeqResults'
plotMA(object, alpha, main = "",
       xlab = "mean of normalized counts", ylim, MLE = FALSE, ...)
