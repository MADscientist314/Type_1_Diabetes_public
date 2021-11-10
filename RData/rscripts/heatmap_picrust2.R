#install.packages("matrixTests")
#############HEATMAPP##############
library(mosaic)
test<-summary2%>%
  pivot_wider(names_from = OTU,values_from="perc_samples")
test<-data.frame(test)
rownames(test)<-paste(test$Delivery,test$disease,sep="_")
test$Delivery<-NULL
test$disease<-NULL
test_mat<-as.matrix(t(test))
test_mat

library(matrixTests)
#write.table(test2,"picrust2_percentsamples.tsv",sep="\t")

df_test<-df5%>%
  select(Abundance,OTU,Delivery,disease,subject,Host)%>%
  pivot_wider(names_from = c(subject,Delivery,Host),values_from=Abundance)%>%select(-disease)

df_test<-df5%>%
  select(Abundance,OTU,Delivery,disease,subject,Host)%>%
  pivot_wider(names_from = OTU,values_from=Abundance)

dim(df_test)


#x = numeric matrix.
#g = vector specifying group membership for each observation of x.
otu_stats<-col_bartlett(x = df_test[,5:176], g=c(paste(df_test$disease,df_test$Delivery,sep="_")))
head(otu_stats)

otu_stats$pbonf <- p.adjust(otu_stats$pvalue, "bonferroni")
sig_otu_stats<-otu_stats%>%filter(pbonf<0.05)
sig_otu_stats<-as_tibble(sig_otu_stats,rownames = "OTU")
sig_otu_stats
summary2<-summary2%>%filter(OTU%in%sig_otu_stats$OTU)
summary2
test<-summary2%>%
  pivot_wider(names_from = OTU,values_from="perc_samples")
test
test<-data.frame(test,check.names = F,stringsAsFactors = T)
test
rownames(test)<-paste(test$Delivery,test$disease,sep="_")
test$Delivery<-NULL
test$disease<-NULL
test_mat<-as.matrix(t(test))
dim(test_mat)




write.table(sig_otu_stats,"sig_path_stats.tsv",sep="\t",row.names = T,quote = F)
write.table(test,"sig_path_perc.tsv",sep="\t",row.names = T,quote = F)
write.table(summary2,"summary2_sig_path.tsv",sep="\t",row.names = T,quote = F)

#install.packages('pheatmap')
library(pheatmap)
colnames(test)
rownames(test)
pheatmap(mat = t(test_mat),fontsize_number = 5,)
pheatmap(mat = test_mat,fontsize_number = 5,)







library(matrixTests)

## Draw heatmaps
pheatmap(test)
pheatmap(test, kmeans_k = 2)
pheatmap(test, scale = "row", clustering_distance_rows = "correlation")
pheatmap(t(test), color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
pheatmap(test, cluster_row = FALSE)
pheatmap(test, legend = FALSE)

# Show text within cells
pheatmap(test, display_numbers = TRUE)
pheatmap(test, display_numbers = TRUE, number_format = "\%.1e")
pheatmap(test, display_numbers = matrix(ifelse(test > 5, "*", ""), nrow(test)))
pheatmap(test, cluster_row = FALSE, legend_breaks = -1:4, legend_labels = c("0",
                                                                            "1e-4", "1e-3", "1e-2", "1e-1", "1"))

# Fix cell sizes and save to file with correct size
pheatmap(test, cellwidth = 15, cellheight = 12, main = "Example heatmap")
pheatmap(test, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "test.pdf")

# Generate annotations for rows and columns
annotation_col = data.frame(
  CellType = factor(rep(c("CT1", "CT2"), 5)), 
  Time = 1:5
)
rownames(annotation_col) = paste("Test", 1:10, sep = "")

annotation_row = data.frame(
  GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
)
rownames(annotation_row) = paste("Gene", 1:20, sep = "")

# Display row and color annotations
pheatmap(test, annotation_col = annotation_col)
pheatmap(test, annotation_col = annotation_col, annotation_legend = FALSE)
pheatmap(test, annotation_col = annotation_col, annotation_row = annotation_row)

# Change angle of text in the columns
pheatmap(test, annotation_col = annotation_col, annotation_row = annotation_row, angle_col = "45")
pheatmap(test, annotation_col = annotation_col, angle_col = "0")

# Specify colors
ann_colors = list(
  Time = c("white", "firebrick"),
  CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
  GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
)

pheatmap(test, annotation_col = annotation_col, annotation_colors = ann_colors, main = "Title")
pheatmap(test, annotation_col = annotation_col, annotation_row = annotation_row, 
         annotation_colors = ann_colors)
pheatmap(test, annotation_col = annotation_col, annotation_colors = ann_colors[2]) 

# Gaps in heatmaps
pheatmap(test, annotation_col = annotation_col, cluster_rows = FALSE, gaps_row = c(10, 14))
pheatmap(test, annotation_col = annotation_col, cluster_rows = FALSE, gaps_row = c(10, 14), 
         cutree_col = 2)

# Show custom strings as row/col names
labels_row = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
               "", "", "Il10", "Il15", "Il1b")

pheatmap(test, annotation_col = annotation_col, labels_row = labels_row)

# Specifying clustering from distance matrix
drows = dist(test, method = "minkowski")
dcols = dist(t(test), method = "minkowski")
pheatmap(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)

# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap(test, clustering_callback = callback)

## Not run: 
# Same using dendsort package
library(dendsort)

callback = function(hc, ...){dendsort(hc)}
pheatmap(test, clustering_callback = callback)

## End(Not run)



