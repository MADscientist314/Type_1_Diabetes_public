#sourcetracker2 input stuff
rm(list = ls())


df2<-as_tibble(read.table("sourcetracker_results2.txt",header = T,sep = "\t"))%>%
  rename(patient=Description)%>%
  rename(disease=Study)%>%
  rename(Sink=SampleType)%>%
  select(SampleID,patient,disease,Sink,Delivery,Antibiotics,
         Anus,Cervix,Introitus,Vagina,Unknown)%>%
  pivot_longer(cols =-c(SampleID,patient,disease,Sink,Delivery,Antibiotics),
               names_to = "Source",
               values_to = "contribution")
df2
df3<-df2%>%
  pivot_longer(cols = c(Source,Sink),names_to="env",values_to="type.label")%>%
  select(-c(SampleID,patient))%>%group_by(disease,Delivery,Antibiotics,env,type.label)%>%
  summarise(contribution=mean(contribution))%>%ungroup()
df3


##########################################################################
#Viznet tutorial



# Libraries ---------------------------------------------------------------
library(visNetwork)
library(geomnet)
library(igraph)

# Data Preparation --------------------------------------------------------
df2
df4<-df2%>%
  select(disease,Delivery,Source,Sink,contribution)%>%
  group_by(disease,Delivery,Source,Sink)%>%
  summarise(contribution=mean(contribution))%>%
  ungroup()
df4

df5<-df4%>%
  mutate(index=rep(1:40),
         Source=gsub("Anus","Rectum",Source),
         Source=gsub("Vagina","Midvaginal",Source))%>%
  pivot_longer(cols = c(Sink,Source),names_to="env",values_to="SampleType")%>%
  mutate(color=gsub("T1D","#AC2121",disease),
    color=gsub("Control","#1E90FF",color),
    image=gsub("Rectum","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",SampleType),
    image=gsub("Cervix","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",image),
    image=gsub("Midvaginal","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",image),
    image=gsub("Introitus","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",image),
    image=gsub("Ear","https://raw.githubusercontent.com/twbs/icons/3e96eef5cf71207d37608532b5f81290d953dd53/icons/ear.svg",image),
    image=gsub("Stool","https://raw.githubusercontent.com/FortAwesome/Font-Awesome/2360bd54ca4abe8e013d424e6679a397e9b717c8/svgs/solid/poop.svg",image),
    image=gsub("Unknown","https://raw.githubusercontent.com/FortAwesome/Font-Awesome/2360bd54ca4abe8e013d424e6679a397e9b717c8/svgs/solid/circle-question.svg",image),
    shape="image",
    shadow=T,
    width=contribution*100)
    
df5
# from=paste0(disease,"_",Delivery,"_",SampleType),
# to=paste0(Sink),

df5
df4_edges<-df5%>%group_by(index,contribution)%>%select(index,SampleType,env,contribution)%>%pivot_wider(names_from =env,values_from=SampleType)%>%ungroup()
df4_edges<-left_join(df4_edges,df5)%>%select(-c(env,SampleType))%>%distinct_all()
df4_edges
data.frame(df4_edges)
df4
df4_nodes<-df4%>%
  #  group_by(disease)%>%
  pivot_longer(cols = c(Sink,Source),names_to="env",values_to="SampleType")%>%
  mutate(id=paste0(disease,"_",Delivery,"_",SampleType))%>%
  select(id,disease,Delivery,SampleType,env)%>%
  distinct_all()
head(df4_nodes)
head(df4_edges)

unique(df4_nodes$id)
unique(df4_edges$to)
#Create graph for Louvain
graph <- graph_from_data_frame(df4_edges, directed = FALSE)

#Louvain Comunity Detection
cluster <- cluster_louvain(graph)
cluster_df <- data.frame(as.list(membership(cluster)))
cluster_df <- as.data.frame(t(cluster_df))

cluster_df<-cluster_df%>%mutate(id=rownames(cluster_df))%>%rename(group=V1)
cluster_df

df4_nodes
#Create group column
df4_nodes <- left_join(df4_nodes, cluster_df, by = "id")
df4_nodes

df4_nodes<-df4_nodes%>%
  mutate(shape="image",
         shadow=T,
         label=gsub("Anus","Rectum",SampleType),
         label=gsub("Vagina","Midvaginal",label),
         image=gsub("Rectum","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",label),
         image=gsub("Cervix","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",image),
         image=gsub("Midvaginal","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",image),
         image=gsub("Introitus","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",image),
         image=gsub("Ear","https://raw.githubusercontent.com/twbs/icons/3e96eef5cf71207d37608532b5f81290d953dd53/icons/ear.svg",image),
         image=gsub("Stool","https://raw.githubusercontent.com/FortAwesome/Font-Awesome/2360bd54ca4abe8e013d424e6679a397e9b717c8/svgs/solid/poop.svg",image),
         image=gsub("Unknown","https://raw.githubusercontent.com/FortAwesome/Font-Awesome/2360bd54ca4abe8e013d424e6679a397e9b717c8/svgs/solid/circle-question.svg",image),
         color=gsub("T1D","#AC2121",disease),
         color=gsub("Control","#1E90FF",color))%>%
  mutate(label=paste0(disease,"_",SampleType),
         label=gsub("Control_Ear","Ear",label),
         label=gsub("T1D_Ear","Ear",label),
         label=gsub("Control_Stool","Stool",label),
         label=gsub("T1D_Stool","Stool",label))%>%
  select(-SampleType)%>%
  distinct_all()

df4_edges
df4_nodes
df4_nodes2<-df4_nodes%>%select(id,label,group,shape,image,shadow,color)
df4_edges<-df4_edges%>%mutate(label=round(width,1))
df4_nodes$label
df4_edges

visNetwork(df4_nodes, df4_edges,
           main = list(text = "Maternal-dyad source sink contributions",
           style = "font-family:Arial;color:#003087;font-size:40px;font-weight:bold;text-align:center;"),
           width = 2100,
           height = 1050,
           submain = list(text = c("Nodes are Louvain clustered and colored by disease.<br> Edges are colored by disease and size weighted by source-sink contribution"),
           style = "font-family:Arial;color:#2b2d42;font-size:24px;text-align:justify;"))%>%
  visNodes(physics = F,shapeProperties = list(useBorderWithImage = F))



# sapply(X =c("pearson","spearman","kendall"),
#        FUN = function(X) 
#          paste0(X,"=",round(cor(x = iris$Sepal.Width,iris$Sepal.Length,use = "everything",method = X),digits = 2)))

       