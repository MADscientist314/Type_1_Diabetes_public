otu_core<-data.frame(otu_table(ASV_physeq_core))
taxa_core<-data.frame(tax_table(ASV_physeq_core))
sample_core<-data.frame(sample_data(ASV_physeq_core))
write.csv(otu_core,'otu_core.csv')
write.csv(taxa_core,'taxa_core.csv')
write.csv(sample_core,'sample_core.csv')


#shiny

install.packages("shiny")
library(shiny)
shiny::runGitHub("shiny-phyloseq","joey711")
