
# Libraries ---------------------------------------------------------------
library(visNetwork)
library(geomnet)
library(igraph)
library(tidyverse)

setwd("K:/github/Type_1_Diabetes_public/RData/sourcetracker/")
#sourcetracker2 input stuff
rm(list = ls())


df2<-as_tibble(read.table("sourcetracker_results.txt",header = T,sep = "\t"))
df2<-df2%>%mutate(patient=Description,
                  disease=Study,
                  Sink=SampleType)%>%
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



# Data Preparation --------------------------------------------------------
#index the samples and fix the nomenclature
#take the mean
df4<-df2%>%
  select(disease,Delivery,Source,Sink,contribution)%>%
  group_by(disease,Delivery,Source,Sink)%>%
  summarise(contribution=mean(contribution))%>%
  mutate(Source=gsub("Anus","Rectum",Source),
         Source=gsub("Vagina","Midvaginal",Source),
         id=paste0(disease,"_",Delivery,"_",Source))%>%
  ungroup()

library(ggsci)
library(scales)
show_col(pal_aaas()(10))
color_pal<-pal_aaas()(7)
df4_edges<-df4%>%
  group_by(id)%>%
  mutate(width=contribution*10,
         label=as.character(round(contribution*100,1)),
        # length=width,
         title=label,
         # dashes=gsub("Ear",TRUE,Sink),
         # dashes=gsub("Stool",FALSE,dashes),
        smooth=T,
        shadow=T)%>%
  ungroup()%>%
  mutate(to=Sink,
         from=id)%>%
  select(from,
         to,
         label,
         width,
         #dashes,
         title,
         smooth,
         shadow
         )%>%distinct_all()
as_tibble(df4_edges)

df4
df4_nodes<-df4%>%
  select(-id)%>%#remove the id column so there are duplicate id
  #pivot the sources and sinks so you can match the images
  pivot_longer(cols = c(Source,Sink),names_to = "name",values_to = "id2")%>%
  #change the Sampletype names to images
  mutate(
    image=gsub("Rectum","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",id2),
    image=gsub("Cervix","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",image),
    image=gsub("Midvaginal","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",image),
    image=gsub("Introitus","https://github.com/MADscientist314/Type_1_Diabetes_public/blob/main/references/Picture2.png?raw=true",image),
    image=gsub("Ear","https://raw.githubusercontent.com/twbs/icons/3e96eef5cf71207d37608532b5f81290d953dd53/icons/ear.svg",image),
    image=gsub("Stool","https://raw.githubusercontent.com/FortAwesome/Font-Awesome/2360bd54ca4abe8e013d424e6679a397e9b717c8/svgs/solid/poop.svg",image),
    image=gsub("Unknown","https://raw.githubusercontent.com/FortAwesome/Font-Awesome/2360bd54ca4abe8e013d424e6679a397e9b717c8/svgs/solid/circle-question.svg",image))%>%
  # fix the nomenclature to match the ids in the df4_edges 
  # gsub the Ear and Stool to get rid of the disease and delivery affiliation
  # adda a label column
  mutate(id=paste0(disease,"_",Delivery,"_",id2),
         id=gsub("Control_Csection_Ear","Neonate Ear",id),
         id=gsub("T1D_Csection_Ear","Neonate Ear",id),
         id=gsub("Control_Vaginal_Ear","Neonate Ear",id),
         id=gsub("T1D_Vaginal_Ear","Neonate Ear",id),
         id=gsub("Control_Csection_Stool","Neonate Stool",id),
         id=gsub("T1D_Csection_Stool","Neonate Stool",id),
         id=gsub("Control_Vaginal_Stool","Neonate Stool",id),
         id=gsub("T1D_Vaginal_Stool","Neonate Stool",id),
         label=id)%>%
  
  #add colors from the disease columns, shadow, and shape column
  mutate(color=case_when(str_detect(label, 'Cervix') ~ '#3B4992FF',
                         str_detect(label, 'Introitus') ~ '#EE0000FF',
                         str_detect(label, 'Midvaginal') ~ '#008B45FF',
                         str_detect(label, 'Rectum') ~ '#631879FF',
                         str_detect(label, 'Ear') ~ '#008280FF',
                         str_detect(label, 'Stool') ~ '#5F559BFF',
                         str_detect(label, 'Unknown') ~ '#1B1919FF'),
         Host=case_when(str_detect(label, 'Ear|Stool') ~ 'Infant',
                        !str_detect(label, 'Ear|Stool') ~ 'Mother'),
         shape="image",
         shadow=TRUE,
         color.background = color,
         color.highlight.background = "black",
         font.color =color)%>%
  #arrange so neonate samples are on the bottom and format to select only the things you need in the order you want
#  mutate(id2=factor(id2,levels = c( "Introitus","Midvaginal","Cervix","Unknown","Rectum", "Ear", "Stool")))%>%
  arrange(desc(Host),label)%>%
  select(id,label,shape,image,shadow,color,color.background,font.color)%>%#,color.highlight.background,font.color)%>%
  separate(col = label,into = "group",sep = "_",remove = F,extra = "drop")%>%
  distinct_all()

tail(df4_nodes$label)
df4_nodes<-df4_nodes%>%separate(id,into=c("disease","Delivery","Source"),remove = F)
tail(df4_nodes)
df4_nodes<-df4_nodes%>%mutate(label=gsub("Control_","",label),
                              label=gsub("T1D_","",label),
                              label=gsub("Csection_","",label),
                              label=gsub("Vaginal_","",label))

#get rid of the control and t1d in the label
df4_nodes<-df4_nodes%>%mutate(label=gsub("_","\n",label),
                              disease=factor(disease,levels = c("Control","T1D","Neonate Ear","Neonate Stool")),
                              Source=factor(Source,levels = c( "Introitus","Midvaginal","Cervix","Rectum","Unknown")))
df4_nodes
# 
# legend <- df4_nodes%>%select() = c("lightblue", "red"),
#                      label = c("reverse", "depends"), arrows =c("to", "from"))



#make x and y coornidates
p<-data.frame(index=rep(1:10))
p
p$x<-as.numeric(sapply(X = p$index, function(Y) 600*cos((Y*(2*pi))/7+(Y))))
p$y<-as.numeric(lapply(X = p$index, function(Y) 600*sin((Y*(2*pi))/7+(Y))))
p<-p%>%arrange(y,x)
p[11,]<-c(11,-300,-25)
p[12,]<-c(12,300,25)
p

library(ggpubr)

tail(df4_nodes)
df4_nodes<-df4_nodes%>%group_by(disease,Source)%>%arrange(disease,desc(Source),Delivery)%>%ungroup()
data.frame(df4_nodes)
df4_nodes$x<-c(p$x[1:10],p$x[1:12])#,p$x[11],p$x[11],p$x[12],p$x[12])
df4_nodes$y<-c(p$y[1:10],p$y[1:12])#,p$y[11],p$y[11],p$y[12],p$y[12])
#df4_nodes<-df4_nodes%>%select(-c(disease,Delivery,Source))
# 
# colnames(df4_nodes)
# data.frame(df4_nodes)
# 
# head(df4_nodes)
# head(df4_edges)
# 
# ####### supplement code to give the edges font color
# 
# 
# # #Create graph for Louvain
# # graph <- graph_from_data_frame(df4_edges, directed = FALSE)
# # #Louvain Comunity Detection
# # cluster <- cluster_louvain(graph)
# # cluster_df <- data.frame(as.list(membership(cluster)))
# # cluster_df <- as.data.frame(t(cluster_df))
# # cluster_df<-cluster_df%>%mutate(id=rownames(cluster_df))%>%rename(group=V1)
# # cluster_df
# # df4_nodes
# # #Create group column
# # df4_nodes <- left_join(df4_nodes, cluster_df, by = "id")
# # df4_nodes
tail(df4_nodes)
df4_edges<-df4_edges%>%mutate(to=paste0("Neonate ",to))

control_edges<-df4_edges[1:20,]
t1d_edges<-df4_edges[21:40,]
control_nodes<-df4_nodes[c(1:10,21:22),]
t1d_nodes<-df4_nodes[c(11:22),]
ggscatter(t1d_nodes,x = "x",y = "y",label = "label")



control_nodes
visNetwork(control_nodes, control_edges,
           main = list(text = "Maternal-dyad source sink contributions",
           style = "font-family:Arial;color:#003087;font-size:40px;font-weight:bold;text-align:center;"),
           footer = "Fig.1 minimal example",
           height = "1200px",
           #width = "100%",
           width = "1200px",
           #height = 1050,
           submain = list(text = c("Nodes are Louvain clustered and colored by disease.<br> Edges are colored by disease and size weighted by source-sink contribution"),
           style = "font-family:Arial;color:#2b2d42;font-size:24px;text-align:justify;"))%>%
 # visGroups(groupname = "Control", color = "#1E90FF")%>%
  # red triangle for group "B"
  visNodes(physics = F,
           #shape = df4_nodes$shape,
           #label = paste0("<b>",df4_nodes$label,
           # color = list(background =df4_nodes$color.background, 
           #              border = df4_nodes$color,
           #              highlight=list(background =df4_nodes$color.background,
           #                             border = df4_nodes$color),
           #              hover=list(background =df4_nodes$color.background,
           #                         border = df4_nodes$color)),
            font=list(
              color=control_nodes$color,
              size=32,face="bold"),#,strokeWidth=1,strokeColor="grey50"),#face=c("arial","bold")),
           # group = df4_nodes$group,
           # x = df4_nodes$x,
           # y=df4_nodes$y,
           # id = df4_nodes$id,
           labelHighlightBold=TRUE,
           shapeProperties = list(useBorderWithImage = F))%>%
  visEdges(dashes = T,label = control_edges$label,
           font=list(
                color=control_nodes$color,
             size=24,face="bold",
             face="arial",weight="bold"))%>%
  visSave(file = "control_sourcetracker_network.html")





visNetwork(t1d_nodes, t1d_edges,
           main = list(text = "Maternal-dyad source sink contributions",
                       style = "font-family:Arial;color:#003087;font-size:40px;font-weight:bold;text-align:center;"),
           footer = "Fig.1 minimal example",
           height = "1200px",
           #width = "100%",
           width = "1200px",
           #height = 1050,
           submain = list(text = c("Nodes are Louvain clustered and colored by disease.<br> Edges are colored by disease and size weighted by source-sink contribution"),
                          style = "font-family:Arial;color:#2b2d42;font-size:24px;text-align:justify;"))%>%
  # visGroups(groupname = "t1d", color = "#1E90FF")%>%
  # red triangle for group "B"
  visNodes(physics = F,
           #shape = df4_nodes$shape,
           #label = paste0("<b>",df4_nodes$label,
           # color = list(background =df4_nodes$color.background, 
           #              border = df4_nodes$color,
           #              highlight=list(background =df4_nodes$color.background,
           #                             border = df4_nodes$color),
           #              hover=list(background =df4_nodes$color.background,
           #                         border = df4_nodes$color)),
           font=list(
             color=t1d_nodes$color,
             size=32,face="bold"),#,strokeWidth=1,strokeColor="grey50"),#face=c("arial","bold")),
           # group = df4_nodes$group,
           # x = df4_nodes$x,
           # y=df4_nodes$y,
           # id = df4_nodes$id,
           labelHighlightBold=TRUE,
           shapeProperties = list(useBorderWithImage = F))%>%
  visEdges(dashes = T,label = t1d_edges$label,
           font=list(
             color=t1d_nodes$color,
             size=24,face="bold",
             face="arial",weight="bold"))%>%
  visSave(file = "t1d_sourcetracker_network.html")

df4_edges
