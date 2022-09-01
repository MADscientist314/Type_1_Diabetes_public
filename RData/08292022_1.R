library(mosaic)
#investifgete 10 and 8.7 vaginal
# 2nd and 3rd trimester a1c and newborn weights
t<-readRDS("ASV_physeq_core.rds")
ASV_physeq_core
sam<-data.frame(sample_data(ASV_physeq_core))%>%mutate(controlled=as_factor(HbA1C_1trym_2>=7.2))
t1dsam<-sam%>%filter(disease=="T1D")%>%group_by(patient,disease)%>%select(patient,disease,HbA1C_1trym_2,Delivery,controlled)%>%distinct_all()
sam
tally(Delivery~controlled,sam)





p<-ggscatter(data = t1dsam,title = "Type 1 Diabetic 1st Trimester HbA1c values x delivery mode colored by controlled <7.2",
             x = "HbA1C_1trym_2",
             y="Delivery",
             color = "controlled",
             palette = my_pal,label = t1dsam$HbA1C_1trym_2,repel = T,
             shape = "Delivery")+rotate_x_text()

ggpar(p = p,xticks.by = 0.2)
tally(~disease,data=t1dsam,"data.frame")%>%filter(Freq>0)



library(mosaic)
tab<-tally(~disease+Delivery+controlled,data=sam,"data.frame")%>%filter(Freq>0)
ggbarplot(data = tab,
          x = "disease",
          y = "Freq",
          fill="controlled",
          label = tab$Freq,
          palette = my_pal,
          position = position_dodge(),
          ylab = "Patients (n)",xlab = "Diabetes",
          facet.by = "Delivery",title = "Diabetes by delivery mode colored by 1st Trimester HbA1c >7.2",
          ggtheme = theme_pubr(legend = "right"))

################################################################################



vec<-as_tibble(data.frame(sample_data(ASV_physeq_core)))%>%group_by(patient,disease)%>%
  filter(disease=="T1D")%>%
  select(patient,disease,HbA1C_1trym_2,Delivery,controlled)%>%distinct_all()
vec%>%filter(controlled=="FALSE")
sam<-sam%>%select(patient,disease,HbA1C_1trym_2,Delivery,controlled)%>%distinct_all()
sam
vec
sam%>%filter(patient==10)
vec%>%filter(patient==10)

tally(disease~controlled+Delivery,vec)
tally(disease~controlled+Delivery,sam)


tmp2<-as_tibble(data.frame(sample_data(t)))%>%select(patient,disease,tydzien_porodu_2,newborn_weight,starts_with("HbA1C_"))%>%distinct_all()%>%select(-c(HbA1C_1trym_2))
tmp2

vec


tmp3<-full_join(vec,tmp2)
tmp3


tmp4<-tmp3%>%
  #filter(Delivery=="Vaginal")%>%
  #filter(HbA1C_1trym_2>7.2)%>%
  pivot_longer(starts_with("HbA1C_"),names_to = "trimester",values_to = "HbA1c")%>%
  mutate(trimester=gsub("HbA1C_1trym_2","1st",trimester),
  trimester=gsub("HbA1C_2trym_2","2nd",trimester),
  trimester=gsub("HBA1C_POROD_2","Delivery",trimester),
  trimester=ordered(trimester,levels = c("1st","2nd","Delivery")))

tmp4<-tmp4%>%mutate(patient=as_factor(patient))
tmp4_fix<-tmp4%>%filter(disease=="Control")%>%select(-HbA1c)
tmp4_fix2<-tmp4%>%filter(disease!="Control")

tmp4<-full_join(tmp4_fix,tmp4_fix2)%>%distinct_all()
tmp4
tmp5<-tmp4%>%filter(disease=="T1D")%>%distinct_all()
tmp5

tmp3
gghistogram(data = tmp3,x = "newborn_weight",color = "controlled",fill = "Delivery")

#favstats on T1D only
# a1c values
# chisq test on delivery
library(gghighlight)

p<-ggplot(tmp5,mapping = aes(x = trimester,y =HbA1c,
                          group=patient,
                          color=controlled)) +
  geom_line() +
  geom_point(mapping=aes(shape=Delivery))+
  scale_color_manual(values = my_pal)+
  scale_fill_manual(values = my_pal)+
  theme_pubr(legend = "right")+facet_grid(~Delivery)+
  ggtitle(label = "Outlier T1D patients with vaginal delivery 1st trimester HbA1c values > 7.2")
ggpar(p = p,yticks.by = 0.2)


vec
sam%>%filter(disease=="T1D")
tmp3.5<-tmp3%>%filter(disease=="T1D")%>%ungroup()%>%select(patient,disease,HbA1C_1trym_2,Delivery,controlled,newborn_weight)%>%distinct_all()

tmp3.5

p<-ggplot(tmp3.5,mapping = aes(y = newborn_weight,x =HbA1C_1trym_2,
                             group=patient,
                             color=controlled)) +
  geom_line() +
  geom_point(mapping=aes(shape=Delivery,))+
  geom_label(mapping = aes(label=newborn_weight))+
  scale_color_manual(values = my_pal)+
  scale_fill_manual(values = my_pal)
  
p
ggpar(p = p,
      legend = "right",
      ggtheme = theme_minimal(),
      xticks.by = 0.2,
      xlab = "1st trimester HbA1c",
      ylab = "newborn weight",
      title = "T1D 1st trimester HbA1c x newborn weight",
      subtitle = "
  color = 1st trimester HbA1c values > 7.2 
  shape = Delivery mode")

tmp3.5
11+22+2+14
