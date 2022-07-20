tmp<-read.csv("h1b.txt", header = T, sep = "\t", row.names = NULL)
a<-as.tibble(tmp)
#a$patient<-as.character(a$patient)
mdata <- melt(a, id=c("patient"))
mdata<-filter(mdata,!is.na(value))

mdata$value
ggplot(data = mdata, mapping = aes(x = variable, y = value), fill=as.character(patient), color=as.character(patient)) + #, group=patient))+
  geom_violin(mapping = aes(x = variable, y = value), alpha=1/4)+
  geom_jitter(width = 0.1, height = 0.1,mapping = aes(color=as.character(patient), fill=as.character(patient), show.legend = F))+
  #  geom_smooth(method = lm, se = FALSE, show.legend = F)+
  #  facet_wrap(nrow = 5, ncol = 10, facets = ~ patient)+
  theme(axis.text = element_text(angle = 90))+
  ggtitle(label = "HbA1C by trimester")

#mapping = aes(x = variable, y = value, color=patient, fill=patient)) +
#geom_violin(mapping = aes(),alpha=1/4)+

#geom_point(mapping = aes(x = variable, y = value, color=patient, fill=patient, group=patient))
#  geom_dotplot(binaxis = "y", stackdir = "center")

a<-lmer(formula =  ~ variable * value, data = mdata)
class(mdata$value)
mean(~ mdata$variable * mdata$value)
library(mosaic)
favstats(variable ~ value,
         data = mdata)

gf_boxplot(variable ~ value,
           data = mdata)
mod<-lm(value ~ variable,data = mdata)
anova(mod)
capture.output(msummary(mod), file = "H1bAC_lm_Anova.tsv")
capture.output(TukeyHSD(mod), file = "H1bAC_Tukey.tsv")
TukeyHSD(mod)