df4_nodes <- load("df4_nodes.Rdata")
df4_edges <- load("df4_edges.Rdata")


control_edges<-df4_edges[1:20,]
t1d_edges<-df4_edges[21:40,]
control_nodes<-df4_nodes[c(1:10,21:22),]
t1d_nodes<-df4_nodes[c(11:22),]
t1d_nodes
control_nodes
#Create graph for Louvain
graph <- graph_from_data_frame(control_edges, directed = FALSE)

#Louvain Comunity Detection
cluster <- cluster_louvain(graph)
cluster_df <- data.frame(as.list(membership(cluster)))
cluster_df <- as.data.frame(t(cluster_df))
cluster_df$label <- rownames(cluster_df)
cluster_df
graph

#Create group column
control_nodes <- left_join(control_nodes, cluster_df, by = "label")


df4_nodes<-df4_nodes%>%
  mutate(shape="image",
         shadow=T,
         label=gsub("Anus","Rectum",label),
         label=gsub("Vagina","Midvaginal",label),
         image=gsub("Rectum","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",label),
         image=gsub("Cervix","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",image),
         image=gsub("Midvaginal","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",image),
         image=gsub("Introitus","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",image),
         image=gsub("Ear","https://raw.githubusercontent.com/twbs/icons/3e96eef5cf71207d37608532b5f81290d953dd53/icons/ear.svg",image),
         image=gsub("Stool","https://raw.githubusercontent.com/FortAwesome/Font-Awesome/2360bd54ca4abe8e013d424e6679a397e9b717c8/svgs/solid/poop.svg",image),
         image=gsub("Unknown","https://raw.githubusercontent.com/FortAwesome/Font-Awesome/2360bd54ca4abe8e013d424e6679a397e9b717c8/svgs/solid/circle-question.svg",image),
         color=gsub("T1D","#AC2121",disease),
         color=gsub("Control","#1E90FF",color))

df4_edges<-df4_edges%>%mutate(label=round(width,1))
visNetwork(df4_nodes, df4_edges,
           main = list(text = "Maternal-dyad source sink contributions",
           style = "font-family:Arial;color:#003087;font-size:40px;font-weight:bold;text-align:center;"),
           width = 2100,
           height = 1050,
           submain = list(text = c("Nodes are Louvain clustered and colored by disease.<br> Edges are colored by disease and size weighted by source-sink contribution"),
           style = "font-family:Arial;color:#2b2d42;font-size:24px;text-align:justify;"))%>%
  visNodes(physics = F,shapeProperties = list(useBorderWithImage = F))
