#################### PACKAGE phyloseq script for Microbial 16 amplicon sequence analysis ###################

### For questions on the scripts below, please contact the author at pietervanveelen2@gmail.com

### DATA ANALYSIS

#load package after successful installation, and other useful packages
library("phyloseq")
library("ggplot2")
library("scales")
library("grid")
library("multcomp")
library("nlme")
library("gridExtra")
library("vegan")
library("biom")
library("PMCRM")
library("MuMIn")


############## CREATE DATAFRAMES, A PHYLOSEQ OBJECT, PRE-PROCESS DATA ######################## 
## REMOVE CHLOROPLAST AND MITOCHONDRIAL OTU'S using phyloseq if not done using QIIME [Current data: DONE IN QIIME WITH filter_taxa_from_otu_table.py]

ASV_physeq_core <- subset_taxa(ASV_physeq_core, Class != "c__Chloroplast")
ASV_physeq_core <- subset_taxa(ASV_physeq_core, Family != "f__mitochondria")
sam<-read.csv(file = "sample_coreJuly15.txt", header = T, row.names = 1, sep = "\t")
sam<-sample_data(sam)
sample_data(ASV_physeq_core)<-sam
sam

# first we need to make a DESeq2 object
deseq_counts<-phyloseq_to_deseq2(ASV_physeq_core, design = ~ disease)
deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
vst_trans_count_tab <- assay(deseq_counts_vst)
# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))
euc_clust <- hclust(euc_dist, method="ward.D2")
euc_dend <- as.dendrogram(euc_clust, hang=0.1)
# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(ASV_physeq_core)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)
sample_variables(vst_physeq)
################### ALPHA DIVERISTY - OTU RICHNESS & SHANNON DIVERSITY ###################################################
# Estimate richness on all sample types
richness_samples5000 <- estimate_richness(ASV_physeq_core)
dim(richness_samples5000)
# sort rownames: cut the "X" and leave the second part of the split
rownames(richness_samples5000) <- sapply(strsplit(rownames(richness_samples5000),"X"),"[",2)

# add sample_id as number from rownames(object)
richness_samples5000[,"sample_id"] <- rownames(richness_samples5000)
richness_samples5000
# order by sample ID
richness_samples5000 <-richness_samples5000[order(richness_samples5000$sample_id),]
meta<-data.frame(sample_data(ASV_physeq_core))

# add metadata to test patterns
library(gridExtra)
meta$sample_id<-rownames(meta)
richness_samples5000<-merge(meta,richness_samples5000, by="sample_id")

Observed_plot <-ggplot(data=richness_samples5000,
                       aes( x=SampleType,y=Observed, fill=disease)) +
  labs(list(x = "SampleType", y = "Observed",fill="disease")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
  geom_boxplot(show.legend=TRUE)

Observed_plot
Shannon_plot <-ggplot(data=richness_samples5000,
                      aes(x=SampleType, y=Shannon, fill=disease)) +
  labs(list(x = "SampleType", y = "Shannon diversity")) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
  geom_boxplot(show.legend=TRUE)

Shannon_plot
grid.arrange(Observed_plot, Shannon_plot, ncol=2)
#dev.off()


vec<-as.vector(sample_variables(ASV_physeq_core))
vec
### ANOVA stats with Observed and Shannon plots per sample type
## Dobson (1990) Page 93: Randomized Controlled Trial :

glm.D93 <- glm(data = richness_samples5000,Observed ~ disease + BMI_before, family = poisson())
anova(glm.D93)
summary(glm.D93)
## Computing AIC [in many ways]:
(A0 <- AIC(glm.D93))
(ll <- logLik(glm.D93))
A1 <- -2*c(ll) + 2*attr(ll, "df")
A2 <- glm.D93$family$aic(counts, mu=fitted(glm.D93), wt=1) +
  2 * length(coef(glm.D93))
stopifnot(exprs = {
  all.equal(A0, A1)
  all.equal(A1, A2)
  all.equal(A1, glm.D93$aic)
})




# A Gamma example, from McCullagh & Nelder (1989, pp. 300-2)
clotting <- data.frame(
  u = c(5,10,15,20,30,40,60,80,100),
  lot1 = c(118,58,42,35,27,25,21,19,18),
  lot2 = c(69,35,26,21,18,16,13,12,12))
summary(glm(lot1 ~ log(u), data = clotting, family = Gamma))
summary(glm(lot2 ~ log(u), data = clotting, family = Gamma))
## Aliased ("S"ingular) -> 1 NA coefficient
(fS <- glm(lot2 ~ log(u) + log(u^2), data = clotting, family = Gamma))
tools::assertError(update(fS, singular.ok=FALSE), verbose=interactive())
## -> .. "singular fit encountered"

## Not run: 
## for an example of the use of a terms object as a formula
demo(glm.vr)


 ## End(Not run)

## phylotype richness:
observed3 <-  lm(Observed~disease*SampleType, richness_samples5000)
anova(observed3)

# drop interaction
observed4 <- lm(Observed~disease+SampleType, richness_samples5000)
anova(observed4)

# check residuals
hist(residuals(observed4))
plot(residuals(observed4) ~ fitted(observed4))

'''
I have no idea what this is doing
'''
# Pairwise contrasts 
library(multcomp)
# Testing group differences per species and per sample type (for Fig. 1a)
richness_samples5000$IntFac = interaction(richness_samples5000$disease, richness_samples5000$Delivery, drop=T)
richness.species.type <- lm(Observed ~ IntFac, richness_samples5000)
summary(richness.species.type)
summary(glht(richness.species.type, mcp(IntFac="Tukey")))


##################################################
## Shannon diversity:
shannon3 <- lm(Shannon~disease*SampleType, richness_samples5000)
anova(shannon3)

# drop interaction
shannon4 <- lm(Shannon~disease+SampleType, richness_samples5000)
anova(shannon4)

# check residuals
hist(residuals(shannon4))

# Pairwise contrasts
# Testing group differences per species and per sample type (for Fig. 1a)
shannon.species.type <- lm(Shannon~IntFac, richness_samples5000)
summary(shannon.species.type)
summary(glht(shannon.species.type, mcp(IntFac="Tukey")))

# Summary stats per species per type
# Observed phylotype richness
library(doBy)
summaryBy(Observed~disease+SampleType,data=richness_samples5000, FUN = function(x) { c(m = mean(x), s = sd(x), n=length(x), sem= sd(x)/sqrt(length(x)))} )

# Shannon diversity
summaryBy(Shannon~Species+Type,data=richness_samples5000, FUN = function(x) { c(m = mean(x), s = sd(x), n=length(x), sem= sd(x)/sqrt(length(x)))} )


## testing individual differences at the Nest level using LMM

library(nlme)
richness_samples5000$
# LMM model with a random intercept for each Nest 
fm1 <- lme(distance ~ age, data = Orthodont) # random is ~ age
fm2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1)
summary(fm1)
summary(fm2)

observed.mod1 <- lme(Observed ~ SampleType,random = ~ 1, data=richness_samples5000, na.action = (na.exclude))

# Linear model with fixed effects Species and Type 
#(note the order, in case of Type I SS = sequential SS: so Type effect already corrected for Species effect)
observed.mod0 <- gls(Observed~SampleType+disease, data=richness_samples5000)

# Test the significance of adding ~1|Nest to the goodness-of-fit using LR-test
anova.lme(observed.mod0, observed.mod1)
Model df      AIC      BIC    logLik   Test  L.Ratio p-value
observed.mod0     1  7 1349.342 1367.852 -667.6708                        
observed.mod1     2  8 1344.836 1365.991 -664.4179 1 vs 2 6.505791  0.0108

# estimate the proportion of variance explained by the fixed (R2m) and fixed+random (R2c) effects
#Nakagawa, S., and H. Schielzeth. 2013. A general and simple method for obtaining R2 from generalized linear mixed-effects models. 
#Methods in Ecology and Evolution 4(2): 133-142. DOI: 10.1111/j.2041-210x.2012.00261.x
library(MuMIn)

r.squaredGLMM(observed.mod0)
R2m       R2c 
0.3996826 0.5372349 
# ~40% of variance within nests explained by Sample Type (also 40% for model without 'Species')


(0.5372349-0.3996826)*100
[1] 13.75523 # % of total variance explained by Nest



## Only bird associated samples ______________________

# filter only bird-associated sample types
richness_samples5000_bird <- subset(richness_samples5000, Delivery== "Csection")# Type==  "broodpatch" | Type==  "feathers"))

# only random intercepts
observed_bird_mod0 <- lme(Observed~1, random = ~1|disease, data=richness_samples5000_bird)

# Linear model with fixed effects Species and Type (note the order, in case of Type I SS = sequential SS: so Type effect already corrected for Species effect)
observed_bird_mod1 <- gls(Observed~SampleType+disease, data = richness_samples5000_bird)

# LMM model with a random intercept for each Nest 
observed_bird_mod2 <- lme(Observed~SampleType+disease, random=~1|disease, data = richness_samples5000_bird)

# Test the significance of adding ~1|Nest to the goodness-of-fit using LR-test
anova(observed_bird_mod1,observed_bird_mod2)

# all sample types
r.squaredGLMM(observed.bird.mod1)
R2m       R2c 
0.3996826 0.5372349 

# variance explained by random effect (% of total variance)
(0.5372349-0.3996826)*100
13.75523
# = ~14%()

# bird only
r.squaredGLMM(observed_bird_mod2)
R2m       R2c 
0.4253293 0.6116700 

# variance explained by random effect (% of total variance)
(0.6116700-0.4253293)*100
[1] 18.63407
# = ~19%


library(plyr)
library(dplyr)
library(tidyverse)
## Only environmental samples
richness_samples_no_eggs_env <- subset(richness_samples5000, SampleType== "Vaginal"| Delivery==  "Vaginal")

summaryBy(Observed~disease,data=richness_samples_no_eggs_env, FUN = function(x) { c(m = mean(x), s = sd(x), n=length(x), sem= sd(x)/sqrt(length(x)), CIlow = mean(x)-(2*(sd(x)/sqrt(length(x)))), CIhigh = mean(x)+(2*(sd(x)/sqrt(length(x)))))} )

# LM 
observed_env_mod0 <- gls(Observed~Species+Type, data=richness_samples_no_eggs_env)

# LMM with nest as random effect
observed_env_mod1 <- lme(Observed~Species+Type, random=~1|Nest, data=richness_samples_no_eggs_env)

anova(observed_env_mod0, observed_env_mod1)


## calculate coefficients of variation within sample type for both species of lark (b=woodlark, v=skylark)

# subset for both species
b <- subset(richness_samples5000, Species=="woodlark")
v <- subset(richness_samples5000, Species=="skylark")

# calculate means
b.mean.observed <- aggregate(b$Observed, by=list(Type=b$Type), mean)
v.mean.observed <- aggregate(v$Observed, by=list(Type=v$Type), mean)

# calculate sd
b.sd.observed <- aggregate(b$Observed, by=list(Type=b$Type), sd)
v.sd.observed <- aggregate(v$Observed, by=list(Type=v$Type), sd)

# merge means and sd for each species
b.observed.summary <- cbind(b.mean.observed, b.sd.observed[,2])
v.observed.summary <- cbind(v.mean.observed, v.sd.observed[,2])

# assign colnames
names(b.observed.summary) <- c("Type", "mean.obs","sd.obs")
names(v.observed.summary) <- c("Type", "mean.obs","sd.obs")

# calculate biased estimate of Coeff Var
b.observed.summary[["coeff.var"]] <- b.observed.summary$sd.obs/b.observed.summary$mean.obs
v.observed.summary[["coeff.var"]] <- v.observed.summary$sd.obs/v.observed.summary$mean.obs

# add N per group
b.observed.summary[["N"]] <- c(11,12,13,13,12)
v.observed.summary[["N"]] <- c(11,13,11,7,7)

# calculate unbiased estimator Coeff Var
b.observed.summary[["unbias.coeff.var"]] <- b.observed.summary$coeff.var*(1+(1/(4*b.observed.summary$N)))
v.observed.summary[["unbias.coeff.var"]] <- v.observed.summary$coeff.var*(1+(1/(4*v.observed.summary$N)))

# add larks as factor
b.observed.summary[["Species"]] <- "B"
v.observed.summary[["Species"]] <- "V"

# merge data frames
observed.summary <- rbind(b.observed.summary, v.observed.summary)
observed.summary$Type <- factor(observed.summary$Type, levels = c("cloacal", "broodpatch","feathers","nest material", "soil", "eggshell"))
# Plot Coeff Var per species for each type

y1 <- observed.summary$unbias.coeff.var
y2 <- shannon.summary$unbias.coeff.var

pdf("Coefficient of variation_per type_per species_shannon1.pdf", width = 6, height = 6)
ggplot(shannon.summary, aes(x=Type, y=unbias.coeff.var, color=Species, group=Species)) +
labs(list(x = "Type", y = "Coefficient of variation")) + 
theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
scale_x_discrete(labels=c("Cloacal gut","Brood patch skin", "Feather", "Nest lining", "Surface soil")) +
scale_colour_manual(values=c("#E69F00", "#56B4E9"), 
name="",
breaks=c("B", "V"),
labels=c("Woodlark", "Skylark")) +
geom_line(show.legend=TRUE, size = 2) +
geom_point(size=2)
dev.off()




################### OTU CO-OCCURRENCE USING VENN DIAGRAMS ###################################################

# create venn diagram for all lark-associated samples #

# input: vst_physeq

## both larks combined
# subset all types individually
library(microbiome)
sample_data(vst_physeq)
meta(vst_physeq)$SampleType
v <- subset_samples(vst_physeq, SampleType=="Vagina")
i <- subset_samples(vst_physeq, SampleType==" Introitus")
c <- subset_samples(vst_physeq, SampleType=="Cervix")
a <- subset_samples(vst_physeq, SampleType=="Anus")
e <- subset_samples(vst_physeq, SampleType=="Ear")
s <- subset_samples(vst_physeq, SampleType=="Stool")

# prune 0 abundance taxa
v.pruned <- prune_taxa(taxa_sums(v) > 0, v)
i.pruned <- prune_taxa(taxa_sums(i) > 0, i)
c.pruned <- prune_taxa(taxa_sums(c) > 0, c)
a.pruned <- prune_taxa(taxa_sums(a) > 0, a)
e.pruned <- prune_taxa(taxa_sums(e) > 0, e)
s.pruned <- prune_taxa(taxa_sums(s) > 0, s)

# extract OTU names
OTUs_v.pruned <- rownames(otu_table(v.pruned))
OTUs_i.pruned <- rownames(otu_table(i.pruned))
OTUs_c.pruned <- rownames(otu_table(c.pruned))
OTUs_a.pruned <- rownames(otu_table(a.pruned))
OTUs_e.pruned <- rownames(otu_table(e.pruned))
OTUs_s.pruned <- rownames(otu_table(s.pruned))

## plot venn diagrams

# cloaca, broodpatch and feathers
pdf("Venn_larks_cl.bp.f_venn.pdf", useDingbats = F)
cl.bp.f_venn <- venn_diagram3(OTUs_v.pruned,OTUs_i.pruned,OTUs_c.pruned, "Vagina", "Introitus", "Cervix", colors= c("green","white","red"), euler=FALSE)
dev.off()
library(VennDiagram)
# nest material, soil
pdf("Venn_larks_nm.s_venn.pdf", useDingbats = F)
nm.s_venn <- venn_diagram2(OTUs_nm.pruned,OTUs_s.pruned, "Nest material", "Soil", colors= c("grey","blue"), euler=FALSE)
dev.off()

# cloacal, nest, soil
pdf("Venn_larks_cl.nm.s_venn.pdf", useDingbats = F)
cl.nm.s_venn <- venn_diagram3(OTUs_cl.pruned,OTUs_nm.pruned,OTUs_s.pruned, "Cloaca", "Nest material", "Soil", colors= c("green","grey","blue"), euler=FALSE)
dev.off()

# brood patch, feather, nest, soil
pdf("Venn_larks_bp.f.nm.s_venn.pdf", useDingbats = F)
bp.f.nm.s_venn <- venn_diagram4(OTUs_bp.pruned, OTUs_f.pruned, OTUs_nm.pruned, OTUs_s.pruned, "Brood patch skin", "Feather", "Nest material", "Soil", colors= c("white","red","grey","blue"), euler=FALSE)
dev.off()



w <- subset_samples(vst_physeq, Larks=="B")
s <- subset_samples(vst_physeq, Larks=="V")

# woodlark only subsets
w.cl <- subset_samples(w, Type=="cloacal")
w.bp <- subset_samples(w, Type=="broodpatch")
w.f <- subset_samples(w, Type=="feathers")
w.nm <- subset_samples(w, Type=="nest material")
w.s <- subset_samples(w, Type=="soil")

# skylark only subsets
s.cl <- subset_samples(s, Type=="cloacal")
s.bp <- subset_samples(s, Type=="broodpatch")
s.f <- subset_samples(s, Type=="feathers")
s.nm <- subset_samples(s, Type=="nest material")
s.s <- subset_samples(s, Type=="soil")

# prune taxa with 0 abundance
w.cl.pruned <- prune_taxa(taxa_sums(w.cl) > 0, w.cl)
w.bp.pruned <- prune_taxa(taxa_sums(w.bp) > 0, w.bp)
w.f.pruned <- prune_taxa(taxa_sums(w.f) > 0, w.f)
w.nm.pruned <- prune_taxa(taxa_sums(w.nm) > 0, w.nm)
w.s.pruned <- prune_taxa(taxa_sums(w.s) > 0, w.s)

s.cl.pruned <- prune_taxa(taxa_sums(s.cl) > 0, s.cl)
s.bp.pruned <- prune_taxa(taxa_sums(s.bp) > 0, s.bp)
s.f.pruned <- prune_taxa(taxa_sums(s.f) > 0, s.f)
s.nm.pruned <- prune_taxa(taxa_sums(s.nm) > 0, s.nm)
s.s.pruned <- prune_taxa(taxa_sums(s.s) > 0, s.s)


# create lists of rownames otu table (OTU names)
OTUs_w.cl.pruned <- rownames(otu_table(w.cl.pruned))
OTUs_w.bp.pruned <- rownames(otu_table(w.bp.pruned))
OTUs_w.f.pruned <- rownames(otu_table(w.f.pruned))
OTUs_w.nm.pruned <- rownames(otu_table(w.nm.pruned))
OTUs_w.s.pruned <- rownames(otu_table(w.s.pruned))

OTUs_s.cl.pruned <- rownames(otu_table(s.cl.pruned))
OTUs_s.bp.pruned <- rownames(otu_table(s.bp.pruned))
OTUs_s.f.pruned <- rownames(otu_table(s.f.pruned))
OTUs_s.nm.pruned <- rownames(otu_table(s.nm.pruned))
OTUs_s.s.pruned <- rownames(otu_table(s.s.pruned))



## venn diagrams

# woodlark only

# cloaca, broodpatch and feathers
pdf("Venn_woodlarks_cl.bp.f_venn.pdf", useDingbats = F)
w.cl.bp.f_venn <- venn_diagram3(OTUs_w.cl.pruned,OTUs_w.bp.pruned,OTUs_w.f.pruned, "Cloaca", "Brood patch skin", "Feather", colors= c("green","white","red"), euler=FALSE)
dev.off()

# nest material, soil
pdf("Venn_woodlarks_nm.s_venn.pdf", useDingbats = F)
w.nm.s_venn <- venn_diagram2(OTUs_w.nm.pruned,OTUs_w.s.pruned, "Nest material", "Soil", colors= c("grey","blue"), euler=FALSE)
dev.off()

# cloacal, nest, soil
pdf("Venn_woodlarks_cl.nm.s_venn.pdf", useDingbats = F)
w.cl.nm.s_venn <- venn_diagram3(OTUs_w.cl.pruned,OTUs_w.nm.pruned,OTUs_w.s.pruned, "Cloaca", "Nest material", "Soil", colors= c("green","grey","blue"), euler=FALSE)
dev.off()

# brood patch, feather, nest, soil
pdf("Venn_woodlarks_bp.f.nm.s_venn.pdf", useDingbats = F)
w.bp.f.nm.s_venn <- venn_diagram4(OTUs_w.bp.pruned, OTUs_w.f.pruned, OTUs_w.nm.pruned, OTUs_w.s.pruned, "Brood patch skin", "Feather", "Nest material", "Soil", colors= c("white","red","grey","blue"), euler=FALSE)
dev.off()

# Skylark only

# cloaca, broodpatch and feathers
pdf("Venn_skylarks_cl.bp.f_venn.pdf", useDingbats = F)
s.cl.bp.f_venn <- venn_diagram3(OTUs_s.cl.pruned,OTUs_s.bp.pruned,OTUs_s.f.pruned, "Cloaca", "Brood patch skin", "Feather", colors= c("green","white","red"), euler=FALSE)
dev.off()

# nest material, soil
pdf("Venn_skylarks_nm.s_venn.pdf", useDingbats = F)
s.nm.s_venn <- venn_diagram2(OTUs_s.nm.pruned,OTUs_s.s.pruned, "Nest material", "Soil", colors= c("grey","blue"), euler=FALSE)
dev.off()

# cloacal, nest, soil
pdf("Venn_skylarks_cl.nm.s_venn.pdf", useDingbats = F)
s.cl.nm.s_venn <- venn_diagram3(OTUs_s.cl.pruned,OTUs_s.nm.pruned,OTUs_s.s.pruned, "Cloaca", "Nest material", "Soil", colors= c("green","grey","blue"), euler=FALSE)
dev.off()

# brood patch, feather, nest, soil
pdf("Venn_skylarks_bp.f.nm.s_venn.pdf", useDingbats = F)
s.bp.f.nm.s_venn <- venn_diagram4(OTUs_s.bp.pruned, OTUs_s.f.pruned, OTUs_s.nm.pruned, OTUs_s.s.pruned, "Brood patch skin", "Feather", "Nest material", "Soil", colors= c("white","red","grey","blue"), euler=FALSE)
dev.off()





################### STACKED BARPLOTS OF RELATIVE ABUNDANCES ###################################################

# agglomerate phylotypes into at phylum level groups
samples5000.phylum <- tax_glom(samples5000, taxrank="Phylum")
samples5000.class <- tax_glom(samples5000, taxrank="Class")

# split by lark species
samples5000.phylum.V <- subset_samples(samples5000.phylum, Larks =="V")
samples5000.phylum.B <- subset_samples(samples5000.phylum, Larks =="B")
samples5000.class.V <- subset_samples(samples5000.class, Larks =="V")
samples5000.class.B <- subset_samples(samples5000.class, Larks =="B")

# calculate mean abundance per sample type
samples_merged.phylum <- merge_samples(samples5000.phylum, "Type", fun=mean)
samples_merged.phylum.V <- merge_samples(samples5000.phylum.V, "Type", fun=mean)
samples_merged.phylum.B <- merge_samples(samples5000.phylum.B, "Type", fun=mean)
samples_merged.class.V <- merge_samples(samples5000.class.V, "Type", fun=mean)
samples_merged.class.B <- merge_samples(samples5000.class.B, "Type", fun=mean)

# calculate relative abundances
samples_merged.phylum.rel <- transform_sample_counts(samples_merged.phylum, function(x) 100* x / sum(x))
samples_merged.phylum.V.rel <- transform_sample_counts(samples_merged.phylum.V, function(x) 100* x / sum(x))
samples_merged.phylum.B.rel <- transform_sample_counts(samples_merged.phylum.B, function(x) 100* x / sum(x))
samples_merged.class.V.rel <- transform_sample_counts(samples_merged.class.V, function(x) 100* x / sum(x))
samples_merged.class.B.rel <- transform_sample_counts(samples_merged.class.B, function(x) 100* x / sum(x))

# reorder facor levels of sample type
sample_data(samples_merged.phylum.rel)$Type <-factor(sample_data(samples_merged.class.V.rel)$Type, levels = c("broodpatch", "cloacal", "feathers", "nest material", "soil"))
sample_data(samples_merged.phylum.V.rel)$Type <-factor(sample_data(samples_merged.class.V.rel)$Type, levels = c("broodpatch", "cloacal", "feathers", "nest material", "soil"))
sample_data(samples_merged.phylum.B.rel)$Type <-factor(sample_data(samples_merged.class.V.rel)$Type, levels = c("broodpatch", "cloacal", "feathers", "nest material", "soil"))
sample_data(samples_merged.class.V.rel)$Type <-factor(sample_data(samples_merged.class.V.rel)$Type, levels = c("broodpatch", "cloacal", "feathers", "nest material", "soil"))
sample_data(samples_merged.class.B.rel)$Type <-factor(sample_data(samples_merged.class.V.rel)$Type, levels = c("broodpatch", "cloacal", "feathers", "nest material", "soil"))

# create a vector with sample type names in the desired order
type.order <- c("cloacal", "broodpatch", "feathers", "nest material", "soil")


sample_data(samples_merged.phylum.V.rel) %>%
slice(match(type.order, Type))
sample_data(samples_merged.phylum.B.rel) %>%
slice(match(type.order, Type))
sample_data(samples_merged.class.V.rel) %>%
slice(match(type.order, Type))
sample_data(samples_merged.class.B.rel) %>%
slice(match(type.order, Type))

# split Phylum Proteobacteria in its Classes
phyla_noProteob.V <- subset_taxa(samples_merged.phylum.V.rel, Phylum !="p__Proteobacteria", na.rm = FALSE) 
phyla_noProteob.B <- subset_taxa(samples_merged.phylum.B.rel, Phylum !="p__Proteobacteria", na.rm = FALSE) 
classes_Proteo.V <- subset_taxa(samples_merged.class.V.rel, Phylum =="p__Proteobacteria", na.rm = FALSE) 
classes_Proteo.B <- subset_taxa(samples_merged.class.B.rel, Phylum =="p__Proteobacteria", na.rm = FALSE) 

otus_V <- merge_phyloseq_pair(otu_table(phyla_noProteob.V),otu_table(classes_Proteo.V))
otus_B <- merge_phyloseq_pair(otu_table(phyla_noProteob.B),otu_table(classes_Proteo.B))

samples_V <- merge_phyloseq_pair(sample_data(phyla_noProteob.V),sample_data(classes_Proteo.V))
samples_B <- merge_phyloseq_pair(sample_data(phyla_noProteob.B),sample_data(classes_Proteo.B))

tax_V <- merge_phyloseq_pair(tax_table(phyla_noProteob.V),tax_table(classes_Proteo.V))
tax_B <- merge_phyloseq_pair(tax_table(phyla_noProteob.V),tax_table(classes_Proteo.V))

# replace the phylum Proteobacteria name tags with their class names
tax_V[21:25,2] <- tax_V[21:25,3]
tax_B[21:25,2] <- tax_B[21:25,3]

#create new phyloseq objects with proteob classes merged with other phyla
physeq.V <- phyloseq(otus_V, samples_V, tax_V)
physeq.B <- phyloseq(otus_B, samples_B, tax_B)

# select only the abundant taxa
selected.V = filter_taxa(physeq.V, function(x) mean(x) > 0.1, prune = TRUE)
selected.B = filter_taxa(physeq.B, function(x) mean(x) > 0.1, prune = TRUE) 

pruned.V <- prune_taxa(names(sort(taxa_sums(selected.V), TRUE))[1:11], selected.V)
pruned.B <- prune_taxa(names(sort(taxa_sums(selected.V), TRUE))[1:11], selected.B)

# Give Type names
sample_data(pruned.V)$Type <- rownames(sample_data(pruned.V))
sample_data(pruned.B)$Type <- rownames(sample_data(pruned.B))

# prune
samples_merged.phylum.rel1 <- prune_taxa(names(sort(taxa_sums(samples_merged.phylum.rel), TRUE)[1:11], samples_merged.phylum.rel)
                 
                 # OPTION
                 phylum.rel.df.V <- psmelt(samples_merged.phylum.V.rel)  
                 phylum.rel.df.B <- psmelt(samples_merged.phylum.B.rel)  
                 class.rel.df.V <- psmelt(samples_merged.class.V.rel) 
                 class.rel.df.B <- psmelt(samples_merged.class.B.rel)     
                 
                 # OPTION write tables to manually replace the Proteobacteria phylum with its classes
                 write.table(phylum.rel.df.V, file = "Relative_abundances_Skylarks_Phyla_allSamples.txt", sep = "\t", dec=".")
                 write.table(phylum.rel.df.B, file = "Relative_abundances_Woodlarks_Phyla_allSamples.txt", sep = "\t", dec=".")
                 write.table(class.rel.df.V, file = "Relative_abundances_Skylarks_Classes_allSamples.txt", sep = "\t", dec=".")
                 write.table(class.rel.df.B, file = "Relative_abundances_Woodlarks_Classes_allSamples.txt", sep = "\t", dec=".")
                 
                 ## RUN "stacked_barplot_phyloseq_joey711.R" first to build the plots.
                 # Google for the R script on Github of joey711 (https://github.com/joey711/phyloseq/issues/442)
                 
                 # barplot of phyla (and proteobacteria classes) for both species across all sample types
                 pdf("Barplot_OTU_rel_abund_perLark_per_Type_ProteoClasses.no.eggs_neworder2.pdf", width = 10, height = 4)
                 V <- plot_ordered_bar(pruned.V, x="Type",fill = "Phylum")  + 
                   scale_y_continuous(name = "Relative abundance") + 
                   labs(title = "Skylark") +
                   scale_x_discrete(labels=c("Cloacal gut","Brood patch skin", "Feather", "Nest lining", "Surface soil")) +
                   scale_fill_manual(values = palette11)
                 B <- plot_ordered_bar(pruned.B, "Type",fill = "Phylum")  + 
                   scale_y_continuous(name = "Relative abundance") + 
                   labs(title = "Woodlark") +
                   scale_x_discrete(labels=c("Cloacal gut","Brood patch skin","Feather", "Nest lining", "Surface soil"))  +
                   scale_fill_manual(values = palette11)
                 grid.arrange(V, B, ncol=2)
                 dev.off()
                 
                 
                 
                 ############## Test for differential abundance with ANCOM ###################
                 
                 # Load all dependecies:
                 
                 install.packages("doParallel")
                 install.packages("DT")
                 install.packages("exactRankTests")
                 install.packages("foreach")
                 install.packages("ggplot2")
                 install.packages("Rcpp")
                 install.packages("shiny")
                 
                 # Install Ancom.R
                 install.packages("~/AncomR/ancom.R_1.1-2.tar.gz")
                 
                 # Load Ancom.R package
                 library("ancom.R")
                 
                 # community table from phyloseq object
                 otu.table <- t(otu_table(samples5000))
                 otu.table <- otu.table[order(as.numeric(rownames(otu.table))),]
                 
                 # Transform by sample ordered transposed matrix into data.frame
                 df_otu_table <- as.data.frame(otu.table)
                 
                 # attach metadata
                 mapping5000 <- as.data.frame(sample_data(samples5000))
                 mapping5000 <- mapping5000[order(mapping5000$X.SampleID),]
                 
                 # Add Sample type column to ordered OTU dataframe
                 df_otu_table$Group <- mapping5000$Type
                 
                 # create separate data frames for each sample type
                 otu_cloacal <- subset(df_otu_table, Group == "cloacal")
                 otu_broodp <- subset(df_otu_table, Group == "broodpatch")
                 otu_feather <- subset(df_otu_table, Group == "feathers")
                 otu_nest <- subset(df_otu_table, Group == "nest material")
                 otu_soil <- subset(df_otu_table, Group == "soil")
                 
                 # delete all grouping columns to remove 'type' specification in 'Group' and replace them by Lark specification
                 otu_cloacal <- otu_cloacal[,-ncol(otu_cloacal)]
                 otu_broodp <- otu_broodp[,-ncol(otu_broodp)]
                 otu_feather <- otu_feather[,-ncol(otu_feather)]
                 otu_nest <- otu_nest[,-ncol(otu_nest)]
                 otu_soil <- otu_soil[,-ncol(otu_soil)]
                 
                 # add Lark labels in factor 'Group'
                 otu_cloacal$Group <- subset(mapping5000_no_eggs, Type == "cloacal")$Larks
                 otu_broodp$Group <- subset(mapping5000_no_eggs, Type == "broodpatch")$Larks
                 otu_feather$Group <- subset(mapping5000_no_eggs, Type == "feathers")$Larks
                 otu_nest$Group <- subset(mapping5000_no_eggs, Type == "nest material")$Larks
                 otu_soil$Group <- subset(mapping5000_no_eggs, Type == "soil")$Larks
                 
                 # Run ANCOM on each OTU table with FDR p-value correction
                 cloacal_ancom.out <- ANCOM(real.data = otu_cloacal, sig = 0.05, multcorr = 2, tau = 0.02, theta = 0.1)
                 broodpatch_ancom.out <- ANCOM(real.data = otu_broodp, sig = 0.05, multcorr = 2, tau = 0.02, theta = 0.1)
                 feather_ancom.out <- ANCOM(real.data = otu_feather, sig = 0.05, multcorr = 2, tau = 0.02, theta = 0.1)
                 nest_ancom.out <- ANCOM(real.data = otu_nest, sig = 0.05, multcorr = 2, tau = 0.02, theta = 0.1)
                 soil_ancom.out <- ANCOM(real.data = otu_soil, sig = 0.05, multcorr = 2, tau = 0.02, theta = 0.1)
                 
                 
                 
                 
                 
                 
                 #################### Unifrac and Bray Curtis distances for within and between sample types ####################
                 
                 library(vegan)
                 
                 # change level order of sample 'Type' factor
                 sample_data(vst_physeq)$Type <- factor(sample_data(vst_physeq)$Type, levels = c("cloacal","broodpatch","feathers" ,"nest material","soil")) 
                 
                 # change colname "Species" to "Larks" to prevent confusion with bacterial Species later
                 colnames(sample_data(vst_physeq))[7] <- "Larks"
                 
                 # calculate weighted UniFrac distance matrix
                 # run adonis PERMANOVA
                 library(vegan)
                 colnames(df) <- "Larks"
                 
                 df = as(sample_data(vst_physeq), "data.frame")
                 wUniFrac  <- UniFrac(vst_physeq, weighted=TRUE, normalized = TRUE)
                 lark.type.adonis = adonis(wUniFrac ~  Larks + Type + RingnumberFemale,strata = df$Larks, df)
                 lark.type.adonis
                 
                 
                 
                 # calculate Bray-Curtis dissimilarity matrix
                 # run adonis PERMANOVA
                 d.bray <- vegdist(t(otu_table(samples5000)), method="bray", binary = FALSE, na.rm = FALSE)
                 lark.bray.adonis = adonis(d.bray ~ Larks, df)
                 lark.bray.adonis
                 
                 # create dataframes with PCoA 1,2,and 3 of Bray-Curtis indices
                 bray_df.12 <- plot_ordination(vst_physeq, bray_within_pcoa.df, axes = c(1,2), type = "samples", color = "Type", justDF=TRUE)
                 bray_df.13 <- plot_ordination(vst_physeq, bray_within_pcoa.df, axes = c(1,3), type = "samples", color = "Type", justDF=TRUE)
                 
                 # subset PCoA data by Lark species    
                 bray_df.B12 <- subset(bray_df.12, Larks=="woodlark")
                 bray_df.B13 <- subset(bray_df.13, Larks=="woodlark")
                 bray_df.V12 <- subset(bray_df.12, Larks=="skylark")
                 bray_df.V13 <- subset(bray_df.13, Larks=="skylark")
                 
                 pdf("FigS3_CH2_Bray_allSamples.pdf", useDingbats = F, width = 12, height = 9)
                 bray.B12 <- ggplot() +
                   geom_point(data=bray_df.B12, aes(x=Axis.1, y=Axis.2, color=Type, size = 3)) +
                   ylim(-0.6, 0.4) +
                   ggtitle("WL12")
                 bray.B13 <- ggplot() +
                   geom_point(data=bray_df.B13, aes(x=Axis.1, y=Axis.3, color=Type, size = 3)) +
                   ylim(-0.6, 0.4)+
                   ggtitle("WL13")
                 bray.V12 <- ggplot() +
                   geom_point(data=bray_df.V12, aes(x=Axis.1, y=Axis.2, color=Type, size = 3)) +
                   ylim(-0.6, 0.4)+
                   ggtitle("SL12")
                 bray.V13 <- ggplot() +
                   geom_point(data=bray_df.V13, aes(x=Axis.1, y=Axis.3, color=Type, size = 3)) +
                   ylim(-0.6, 0.4)+
                   ggtitle("SL13")
                 grid.arrange(bray.B12,bray.V12,bray.B13,bray.V13, nrow=2, ncol=2)
                 dev.off()
                 
                 #wUniFrac Larks
                 wUniFr_within_pcoa <- ordinate(vst_physeq, method = "PCoA", distance ="unifrac", weighted = TRUE)
                 wUniFr_df <- plot_ordination(vst_physeq, wUniFr_within_pcoa, axes = c(1,2), type = "samples", color = "Type", justDF=TRUE)
                 
                 # subset PCoA data by Lark species   
                 wUniFr_df.B <- subset(wUniFr_df, Larks=="woodlark")
                 wUniFr_df.V <- subset(wUniFr_df, Larks=="skylark")
                 
                 pdf("Fig2_CH2_wUniFrac_allSamples.pdf", useDingbats = F, width = 6, height = 9)
                 wUnifrac.B <- ggplot() +
                   geom_point(data=wUniFr_df.B, aes(x=Axis.1, y=Axis.2, color=Type, size = 3)) +
                   ylim(-0.7, 0.3)
                 wUnifrac.V <- ggplot() +
                   geom_point(data=wUniFr_df.V, aes(x=Axis.1, y=Axis.2, color=Type, size = 3)) +
                   ylim(-0.7, 0.3) 
                 grid.arrange(wUnifrac.B,wUnifrac.V, nrow=2)
                 dev.off()
                 
                 
                 
                 
                 ################ Phylogenetic dispersion of samples in PCoA based on weighted UniFrac ##########################
                 
                 # create separate phyloseq objects for both lark species
                 samples5000.W <- subset_samples(ASV_physeq_core, disease=="Control")
                 samples5000.S <- subset_samples(ASV_physeq_core, Larks=="T1D")
                 
                 # Make pairwise distance matrices using Weighted UniFrac method
                 samples_.wUniFr.W <- UniFrac(samples5000.W, weighted=TRUE, normalized = TRUE)
                 samples_.wUniFr.S <- UniFrac(samples5000.S, weighted=TRUE, normalized =TRUE)
                 
                 # create a long format df with weighted UniFrac distances for each combination of samples
                 m.W <- as.matrix(samples_.wUniFr.W) 
                 m2.W <- melt(m.W)[melt(upper.tri(m.W))$value,]
                 names(m2.W) <- c("sample1", "sample2", "distance")
                 m2.W[["Larks"]] <- "woodlark"
                 
                 m.S <- as.matrix(samples_.wUniFr.S) 
                 m2.S <- melt(m.S)[melt(upper.tri(m.S))$value,]
                 names(m2.S) <- c("sample1", "sample2", "distance")
                 m2.S[["Larks"]] <- "skylark"
                 
                 # collate dataframes of both larks
                 wUniFr_distances_larks <- rbind(m2.W,m2.S)
                 
                 # Make dataframe for joining SampleIDs in 'sample1' with SampleType from sample_data(samples5000) 
                 sampletype1_join <- sample_data(vst_physeq)[,c(1,6)]
                 sampletype1_join <- sampletype1_join@.Data
                 sampletype1_join <- as.data.frame(sampletype1_join)
                 colnames(sampletype1_join) <- c("sample1" ,"Sample1Type")
                 
                 require(plyr)
                 
                 # join sample types for sample1 column
                 wUniFr_distances_larks2 <- join(wUniFr_distances_larks, sample1_join, by="sample1",type="left")
                 
                 # repeat for 'sample2' column
                 sample2_join <- sample1_join
                 colnames(sample2_join)[2] <- "Sample2Type"
                 colnames(sample2_join) <- c("sample2", "Sample2Type")
                 wUniFr_distances_larks3 <- join(wUniFr_distances_larks2, sample2_join, by="sample2",type="left")
                 
                 # save resulting dataframe
                 write.csv(wUniFr_distances_larks3, file="DATA_weighted UniFrac distances for sample Types_no_eggs.csv")
                 
                 # IN EXCEL: 
                 # Sort dataframe in Excel by Sample1Type, Sample2Type, sample1, sample2
                 # Manually added factor "within" type comaprison 'yes/no', "between" 'yes/no', for comparisons within or between sample types, resp.
                 # add factor "concatenated varnames" e.g. 'cloacal - feathers' for sample types of sample1 and sample2 for each comparison, resp.
                 
                 # READ new csv file:
                 wUniFr_betw_types <- as.data.frame(read.csv("Data_weighted_UniFrac_dist_sorted_concatVarname.csv", header=T, sep=";"))
                 
                 # within versus between type variation
                 
                 # linear model to test differential within type distances versus between type distances
                 within.vs.between.lm <- lm(distance ~ within*Larks, wUniFr_betw_types)
                 
                 # ANOVA
                 anova(within.vs.between.lm)
                 summary(glht(within.vs.between.lm, within="Tukey"))
                 
                 # pairwise summary for interaction of 'within' and 'larks' 
                 wUniFr_betw_types$IntFac = interaction(wUniFr_betw_types$within, wUniFr_betw_types$Larks, drop=T)
                 within.vs.between.lm.final <- lm(distance~IntFac, wUniFr_betw_types)
                 summary(within.vs.between.lm.final)
                 summary(glht(within.vs.between.lm.final, mcp(IntFac="Tukey")))
                 
                 # pairwise summary for interaction of 'each sample type by sample type comparison' and 'larks' 
                 wUniFr_betw_types$IntFac2 = interaction(wUniFr_betw_types$concat_varname, wUniFr_betw_types$Larks, drop=T)
                 within.vs.between.lm.final <- lm(distance~IntFac2, wUniFr_betw_types)
                 summary(within.vs.between.lm.final)
                 summary(glht(within.vs.between.lm.final, mcp(IntFac2="Tukey")))
                 
                 
                 # plot boxplot for all concatenated varnames 
                 
                 # subset rows containing sample comparisons between different sample types
                 between_wUniFr_betw_types <- subset(wUniFr_betw_types, between == "yes")
                 
                 # rename double paired level names in concat_varname
                 wUniFr_betw_types$concat_varname <- mapvalues(wUniFr_betw_types$concat_varname, from=c(
                   "broodpatch - cloacal",
                   "feathers - broodpatch", 
                   "nest material - broodpatch", 
                   "soil - broodpatch", 
                   "feathers - cloacal", 
                   "nest material - cloacal",
                   "soil - cloacal",
                   "nest material - feathers",
                   "soil - feathers", 
                   "soil - nest material"
                 ),
                 to=c(
                   "cloacal - broodpatch", 
                   "broodpatch - feathers", 
                   "broodpatch - nest material", 
                   "broodpatch - soil", 
                   "cloacal - feathers",
                   "cloacal - nest material",
                   "cloacal - soil",
                   "feathers - nest material",
                   "feathers - soil",
                   "nest material - soil"
                 ))
                 
                 # test if no. of levels is correct. should be 15.
                 levels(wUniFr_betw_types$concat_varname)
                 
                 # reorder pairwise comparisons
                 between_wUniFr_betw_types$concat_varname <- factor(between_wUniFr_betw_types$concat_varname, 
                                                                    levels = c("cloacal - broodpatch",
                                                                               "cloacal - feathers",
                                                                               "cloacal - nest material",
                                                                               "cloacal - soil",
                                                                               "broodpatch - feathers",
                                                                               "broodpatch - nest material",
                                                                               "broodpatch - soil",
                                                                               "feathers - nest material",
                                                                               "feathers - soil",
                                                                               "nest material - soil"))
                 
                 
                 # Plot within and between sample type dispersion for both larks (weighted UniFrac distances)
                 pdf("Weighted Unifrac distances between types.3.pdf", width = 10, height = 6, useDingbats = F ) 
                 ggplot(data=between_wUniFr_betw_types,
                        aes(x=concat_varname, y=distance, fill=Larks)) +
                   labs(list(x = "Pairwise sample type comparison", y = "Normalized weighted UniFrac distance")) + 
                   theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
                   scale_fill_manual(values=c("#E69F00", "#56B4E9"), 
                                     name="",
                                     breaks=c("woodlark", "skylark"),
                                     labels=c("Woodlark", "Skylark")) +
                   geom_boxplot(show.legend=TRUE)
                 dev.off()
                 
                 # statistics of phylogenetic dispersion       
                 lm.pairw.wUniFr <- lm(distance~concat_varname*Larks,between_wUniFr_betw_types)
                 lm.pairw.wUniFr2 <- lm(distance~concat_varname+Larks,between_wUniFr_betw_types)
                 
                 anova(lm.pairw.wUniFr)
                 
                 hist(residuals(lm.pairw.wUniFr))
                 
                 between_wUniFr_betw_types$IntFac = interaction(between_wUniFr_betw_types$concat_varname, between_wUniFr_betw_types$Larks, drop=T)
                 lm.pairw.wUniFr3 <- lm(distance~IntFac, between_wUniFr_betw_types)
                 summary(glht(lm.pairw.wUniFr, mcp(IntFac="Tukey")))
                 
                 
                 ### Phylogenetic dispersion WITHIN sample types (based on weighted UniFrac) ####
                 
                 
                 # plot boxplot for all concatenated varnames: 
                 
                 # subset all within-sample-type comparisons
                 within_wUniFr_within_types <- subset(wUniFr_betw_types, within == "yes")
                 
                 # reorder pairwise comparisons
                 within_wUniFr_within_types$concat_varname <- factor(within_wUniFr_within_types$concat_varname, levels = c("cloacal - cloacal","broodpatch - broodpatch","feathers - feathers","nest material - nest material","soil - soil"))
                 
                 
                 # Plot within and between variation for both larks in Weighted Unifrac distances
                 pdf("Weighted Unifrac distances within types.2.pdf", width = 8, height = 6,useDingbats = F) 
                 ggplot(data=within_wUniFr_within_types,
                        aes(x=concat_varname, y=distance, fill=Larks)) +
                   labs(list(x = "Pairwise sample type comparison", y = "Normalized weighted UniFrac distance")) + 
                   theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
                   scale_fill_manual(values=c("#E69F00", "#56B4E9"), 
                                     name="",
                                     breaks=c("woodlark", "skylark"),
                                     labels=c("Woodlark", "Skylark")) +
                   geom_boxplot(show.legend=TRUE)
                 dev.off()
                 
                 
                 # linear model of within type variance
                 
                 lm.within.type <- lm(distance~concat_varname*Larks, within_wUniFr_within_types)
                 
                 anova(lm.within.type)
                 
                 hist(residuals(lm.within.type))
                 
                 within_wUniFr_within_types$IntFac = interaction(within_wUniFr_within_types$concat_varname, within_wUniFr_within_types$Larks, drop=T)
                 lm.within.type.final <- lm(distance~IntFac, within_wUniFr_within_types)
                 summary(lm.within.type.final)
                 summ.within <- summary(glht(lm.within.type.final, mcp(IntFac="Tukey")))
                 summ.within
                 
                 # test differences in phylogenetic dispersion between females
                 within_wUniFr_within_types$IntFac = interaction(within_wUniFr_within_types$concat_varname, within_wUniFr_within_types$Larks, drop=T)
                 lm.within.type.final <- lm(distance~IntFac, within_wUniFr_within_types)
                 lm.within.with.nest <- lmer(distance ~ IntFac + (1|Nest), REML=F, within_wUniFr_within_types)
                 
                 AIC(lm.within.type.final, lm.within.with.nest)
                 
                 
                 
                 
                 ########## Partitioning around medoids ###############################
                 
                 ## partitioning around medoids ##
                 #
                 #following :http://statweb.stanford.edu/~susan/papers/Pregnancy/PNAS_Vaginal_Analysis.html
                 
                 library(phyloseq)
                 library(vegan)
                 library("cluster"); packageVersion("cluster")
                 library("igraph"); packageVersion("igraph")
                 library("markovchain"); packageVersion("markovchain")
                 library(ggpubr)
                 wUniFrac <- phyloseq::distance(vst_physeq, method="bray")
                 ord = ordinate(vst_physeq, method = "MDS", distance = "bray")
                 
                 pdf("EXTRA_plot_scree_PAM.pdf", useDingbats = F)
                 plot_scree(ord) + ggtitle("MDS-bray ordination eigenvalues")
                 dev.off()
                 
                 evs <- ord$value$Eigenvalues
                 print(evs[1:20])
                 [1] 5.9620421 2.8610595 2.0894529 0.9535996 0.8715639 0.7434830 0.5914184
                 [8] 0.5050820 0.3821053 0.3703655 0.3268519 0.2772767 0.2518184 0.2342529
                 [15] 0.2143129 0.1961665 0.1777192 0.1681265 0.1563137 0.1489535
                 
                 print(tail(evs))
                 
                 h_sub5 <- hist(evs[6:length(evs)], 100)
                 pdf("EXTRA_h_sub5.pdf", useDingbats = F)
                 plot(h_sub5$mids, h_sub5$count, log="y", type='h', lwd=10, lend=2)
                 dev.off()
                 
                 
                 NDIM <- 8
                 x <- ord$vectors[,1:NDIM]  # rows=sample, cols=MDS axes, entries = value
                 pamPCoA = function(x, k) {
                   list(cluster = pam(x[,1:2], k, cluster.only = TRUE))
                 }
                 gs = clusGap(x, FUN = pamPCoA, K.max = 12, B = 50)
                 pdf("EXTRA_clusgap_PAM.pdf", useDingbats = F)
                 plot_clusgap(gs) + scale_x_continuous(breaks=c(seq(0, 12, 2)))
                 dev.off()
                 
                 # levelling after 7 clusters
                 
                 K <- 5
                 x <- ord$vectors[,1:NDIM]
                 clust <- as.factor(pam(x, k=K, cluster.only=T))
                 # SWAPPING THE ASSIGNMENT OF 2 AND 3 TO MATCH RAVEL CST ENUMERATION
                 #clust[clust==2] <- NA
                 #clust[clust==3] <- 2
                 #clust[is.na(clust)] <- 3
                 sample_data(ASV_physeq_core)$CST <- clust
                 CSTs <- as.character(seq(K))
                 library(RColorBrewer)
                 CSTColors <- brewer.pal(8,"Paired")[c(1,3,2,5,4)] # Length 7 for consistency with pre-revision CST+ coloration
                 names(CSTColors) <- CSTs
                 CSTColorScale <- scale_colour_manual(name = "CST", values = CSTColors[1:5])
                 CSTFillScale <- scale_fill_manual(name = "CST", values = CSTColors[1:5])
                 library(gridExtra)
                 pdf("EXTRA_PAM_Clustering.pdf", useDingbats = F, width = 10, height = 12)
                a <- plot_ordination(ASV_physeq_core, ord, color="CST", shape="disease") + facet_wrap( ~ SampleType) + CSTColorScale
                b <- plot_ordination(ASV_physeq_core, ord, axes=c(3,4), color="CST", shape="disease") +  facet_wrap( ~ SampleType) + CSTColorScale
                 grid.arrange(a,b, nrow=2)
                 
                 dev.off()
                 
                 pdf("EXTRA_plot_ordination_wUNiFrac_NMDS.pdf", useDingbats = F)
                 c<-plot_ordination(vst_physeq, ordinate(vst_physeq, method="NMDS", distance=wUniFrac), color="CST") + CSTColorScale + ggtitle("NMDS -- wUniFrac -- By Cluster")
                 d<-plot_ordination(vst_physeq, ordinate(vst_physeq, method="NMDS", distance=wUniFrac), color="SampleType", shape = "disease") + ggtitle("NMDS -- wUniFrac -- By Cluster")
                 grid.arrange(c,d,nrow=2)
                 dev.off()
                 
                 
                 
                 
                 ################ TEST OF VARIANCE-ADJUSTED UNIFRAC vs. WEIGHTED UNIFRAC ####################
                 
                 ## Variance-adjusted UniFrac ##
                 
                 # we tested the qualitative interpretation of VW-UniFrac on grouping of Host species (Larks) and sample type (Type)
                 # we followed the fillowing script: https://rdrr.io/cran/GUniFrac/man/GUniFrac.html
                 # input phyloseq object = vst_physeq (rarefied data to 5000 reads/sample)
                 
                 require(ade4)
                 
                 otu.tab <- t(otu_table(vst_physeq))
                 treefile <- phy_tree(vst_physeq)
                 
                 # Calculate the UniFracs
                 unifracs <- GUniFrac(otu.tab, treefile, alpha=c(0, 0.5, 1))$unifracs
                 
                 Larks <- as(sample_data(vst_physeq), "data.frame")$Larks
                 Type <- as(sample_data(vst_physeq), "data.frame")$Type
                 
                 # create 
                 dw <- unifracs[, , "d_1"]		  # Weighted UniFrac
                 du <- unifracs[, , "d_UW"]		# Unweighted UniFrac	
                 dv <- unifracs[, , "d_VAW"]		# Variance adjusted weighted UniFrac
                 d0 <- unifracs[, , "d_0"]     # GUniFrac with alpha 0  
                 d5 <- unifracs[, , "d_0.5"]   # GUniFrac with alpha 0.5 
                 
                 
                 # test using perMANOVA
                 library(vegan)
                 
                 adonis(as.dist(dw) ~ Larks + Type)
                 adonis(as.dist(dv) ~ Larks + Type)
                 adonis(as.dist(d5) ~ Larks + Type)
                 
                 # RESULTS
                 adonis(as.dist(dw) ~ Larks + Type)
                 
                 
                 # distance matrices
                 dist.dw <- as.dist(dw)
                 dist.dv <- as.dist(dv)
                 dist.d5 <- as.dist(d5)
                 
                 # ordinate dm's
                 ord.dw <- ordinate(vst_physeq, method="PCoA", dist.dw)
                 ord.dv <- ordinate(vst_physeq, method="PCoA", dist.dv)
                 ord.d5 <- ordinate(vst_physeq, method="PCoA", dist.d5)
                 
                 # PCoA plots
                 library(gridExtra)
                 pdf("EXTRA_PCoA_different_UniFrac_metrics.pdf", useDingbats = F, width = 16, height = 5)
                 plot.dw <- plot_ordination(vst_physeq, ord.dw, axes = c(1,2), type = "samples", color = "Type", shape = "Larks") + guides(colour = "none")
                 plot.dv <- plot_ordination(vst_physeq, ord.dv, axes = c(1,2), type = "samples", color = "Type", shape = "Larks") + guides(colour = "none")
                 plot.d5 <- plot_ordination(vst_physeq, ord.d5, axes = c(1,2), type = "samples", color = "Type", shape = "Larks")
                 grid.arrange(plot.dw, plot.dv, plot.d5, ncol=3)
                 dev.off()
                 
                 pcoa_dv <- s.class(cmdscale(d5, k=2), fac = Larks) 
                 pcoa_d5 <- s.class(cmdscale(d5, k=2), fac = Larks) 
                 grid.arrange(plot.dw, plot.dv, plot.d5, ncol=3)
                 
                 
                 
                 
                 ############## PHYLOGENETIC STRUCTURE - MNTD and NTI calculations ########################
                 
                 ### NULL model = independent swap method (null model 4 in Kembel, 2009) ###
                 
                 library(phyloseq)
                 library(picante)
                 
                 w <- subset_samples(vst_physeq, Larks=="B")
                 s <- subset_samples(vst_physeq, Larks=="V")
                 
                 cl <- subset_samples(vst_physeq, Type=="cloacal")
                 bp <- subset_samples(vst_physeq, Type=="broodpatch")
                 f <- subset_samples(vst_physeq, Type=="feathers")
                 nm <- subset_samples(vst_physeq, Type=="nest material")
                 s <- subset_samples(vst_physeq, Type=="soil")
                 
                 samp.w <- t(otu_table(w))@.Data
                 samp.s <- t(otu_table(s))@.Data
                 
                 samp.cl <- t(otu_table(cl))@.Data
                 samp.bp <- t(otu_table(bp))@.Data
                 samp.f <- t(otu_table(f))@.Data
                 samp.nm <- t(otu_table(nm))@.Data
                 samp.s <- t(otu_table(s))@.Data
                 
                 phy <- phy_tree(samples5000)
                 # remove al redundant tips and nodes in the phy tree given the OTU's absence in samp 
                 # (always perform this after subsetting 'samp' based on metadata category)
                 # filter all 0 OTU total counts out
                 w.pruned <- prune_taxa(taxa_sums(w) > 0, w)
                 s.pruned <- prune_taxa(taxa_sums(s) > 0, s)
                 
                 cl.pruned <- prune_taxa(taxa_sums(cl) > 0, cl)
                 bp.pruned <- prune_taxa(taxa_sums(bp) > 0, bp)
                 f.pruned <- prune_taxa(taxa_sums(f) > 0, f)
                 nm.pruned <- prune_taxa(taxa_sums(nm) > 0, nm)
                 s.pruned <- prune_taxa(taxa_sums(s) > 0, s)
                 
                 
                 phy.w <- phy_tree(w.pruned)
                 phy.s <- phy_tree(s.pruned)
                 
                 phy.cl <- phy_tree(cl.pruned)
                 phy.bp <- phy_tree(bp.pruned)
                 phy.f <- phy_tree(f.pruned)
                 phy.nm <- phy_tree(nm.pruned)
                 phy.s <- phy_tree(s.pruned)
                 
                 # We also need to make sure the species are arranged in the some order in the community data and the phylogeny. 
                 # This is an important step - several functions in picante assume that the community or trait data and phylogeny data have species arranged in the same order, so it's good to always make sure we've done so before running any analysis. 
                 # The following command sorts the columns of samp to be in the same order as the tip labels of the phylogeny:
                 samporder.cl <- samp.cl[, phy.cl$tip.label]
                 samporder.bp <- samp.bp[, phy.bp$tip.label]
                 samporder.f <- samp.f[, phy.f$tip.label]
                 samporder.nm <- samp.nm[, phy.nm$tip.label]
                 samporder.s <- samp.s[, phy.s$tip.label]
                 samporder.cl[1:10,1:10]
                 
                 # Let's visualize our data. Now let's see how taxa from the six communities in the Phylocom example data set are arranged on the tree. The following commands set up the layout of the plot to have 2 rows and 3 columns, and then plot a black dot for the species present in each of the six communities:
                 
                 #pdf("Visualize_taxa_on_tree_NTI_calc.pdf", useDingbats = F)
                 #par(mfrow = c(2, 3))
                 #for (i in row.names(samporder)) {
                 #  plot(prunedphy, show.tip.label = FALSE, main = i)
                 #  tiplabels(tip = which(samporder[i, ] > 0), pch = 19, cex = 2)
                 #}
                 #dev.off()
                 
                 #Now we will calculate NTI (Nearest Taxon Index) for our different communities. Rembember, NTI values indicate phylogenetic overdispersion, and positive NTI values indicate phylogenetic clustering.
                 #First we need to make a phylogenetic distance matrix.
                 
                 phydist.cl <- cophenetic(phy.cl)
                 phydist.bp <- cophenetic(phy.bp)
                 phydist.f <- cophenetic(phy.f)
                 phydist.nm <- cophenetic(phy.nm)
                 phydist.s <- cophenetic(phy.s)
                 
                 # Take a look at phydist.  This is a matrix where the rows and columns are the taxa and the elements of the matrix are the phylogenetic distance between those pairs of taxa.  
                 
                 
                 # The fifth column is the rank of the score against randomized communities; the sixth column is the negative NRI; and the seventh column is the one tailed p-value for significantly high NRI.  So that communities 1, 2 and 3 are significantly clustered and community 5 is significantly spread out.
                 
                 # To calculate the NTI:
                 MNTD_NTI_all_cl5000 <- ses.mntd(samporder.cl, phydist.cl,null.model="independentswap",runs = 999, iterations = 1000, abundance.weighted=FALSE)
                 MNTD_NTI_all_bp5000 <- ses.mntd(samporder.bp, phydist.bp,null.model="independentswap", runs = 999, iterations = 1000, abundance.weighted=FALSE)
                 MNTD_NTI_all_f5000 <- ses.mntd(samporder.f, phydist.f,null.model="independentswap", runs = 999, iterations = 1000, abundance.weighted=FALSE)
                 MNTD_NTI_all_nm5000 <- ses.mntd(samporder.nm, phydist.nm,null.model="independentswap", runs = 999, iterations = 1000, abundance.weighted=FALSE)
                 MNTD_NTI_all_s5000 <- ses.mntd(samporder.s, phydist.s,null.model="independentswap", runs = 999, iterations = 1000, abundance.weighted=FALSE)
                 
                 # merge results
                 MNTD_NTI_all_persample5000 <- rbind(MNTD_NTI_all_cl5000,MNTD_NTI_all_bp5000,MNTD_NTI_all_f5000,MNTD_NTI_all_nm5000,MNTD_NTI_all_s5000)
                 
                 # transform z-scores in NTI by multiplying by -1
                 MNTD_NTI_all_persample5000[["NTI"]] <- MNTD_NTI_all_persample5000$mntd.obs.z*-1 
                 
                 # Results stored in write.table(MNTD_NTI_all_samples5000, "NTI_MNTD_null_vs_obs_allSamples5000_20161219.xlsx", row.names=TRUE,col.names=T,sep="\t")
                 write.table(MNTD_NTI_all_persample5000, "NTI_MNTD_null_vs_obs_all_persample5000_20170915.xlsx", row.names=TRUE,col.names=T,sep="\t")
                 
                 # add NTI data to metadata samples5000
                 mapping_NTI_persampletype <- sample_data(samples5000)
                 colnames(mapping_NTI_persampletype)[1] <- "X.SampleID"
                 # colnames(mapping_NTI_persampletype)[1] # "sampleID"
                 MNTD_NTI_all_persample5000$sampleID <- rownames(MNTD_NTI_all_persample5000)
                 nti_persampletype <- as.data.frame(mapping_NTI_persampletype@.Data)
                 colnames(nti_persampletype) <- mapping_NTI_persampletype@names
                 colnames(nti_persampletype)[1] <- "sampleID"
                 nti_persampletype[nti_persampletype$sampleID=="192",7] <- "B"
                 NTIdata_persampletype[NTIdata_persampletype$sampleID=="54",7] <- "B"
                 
                 
                 nrow(nti_persampletype) # 175
                 nrow(MNTD_NTI_all_persample5000) # 110
                 
                 
                 # --> copy data from MNTD_NTI_all_persample5000 to nti_persampletype for only matiching sampleID values
                 NTIdata_persampletype <- merge(nti_persampletype, MNTD_NTI_all_persample5000, "sampleID")
                 NTIdata_persampletype[NTIdata_persampletype$sampleID=="54",7] <- "B"
                 NTIdata_persampletype[NTIdata_persampletype$sampleID=="118",10] <- "B14.11"
                 dim(NTIdata_persampletype) # 110 21
                 
                 NTIdata_persampletype$Type <- factor(NTIdata_persampletype$Type, levels = c("cloacal", "broodpatch", "feathers", "nest material", "soil"))
                 
                 # Analysis of variance of NTI 
                 # lm with Larks x Type interaction
                 
                 lm_NTI_persampletype <- lm(NTI~Larks*Type, data=NTIdata_persampletype)
                 
                 anova(lm_NTI_persampletype)
                 
                 
                 # create contrast matrix for each level of Larks for each level of Type
                 
                 NTIdata_persampletype$Type <- factor(NTIdata_persampletype$Type, levels = c("cloacal", "broodpatch", "feathers", "nest material", "soil"))
                 levels(NTIdata_persampletype$Type)
                 
                 library(multcomp)
                 
                 X <- model.matrix(~ Larks * Type, data = NTIdata_persampletype)
                 
                 Tukey <- contrMat(table(NTIdata_persampletype$Type), "Tukey")
                 K1 <- cbind(Tukey, matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)))
                 rownames(K1) <- paste(levels(NTIdata_persampletype$Larks)[1], rownames(K1), sep = ":")
                 K2 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), Tukey)
                 rownames(K2) <- paste(levels(NTIdata_persampletype$Larks)[2], rownames(K2), sep = ":")
                 
                 K <- rbind(K1,K2)
                 colnames(K) <- c(colnames(Tukey), colnames(Tukey))  
                 
                 glht_NTI <- summary(glht(lm_NTI_persampletype, linfct = K))
                 glht_NTI
                 
                 
                 
                 # Calculate summary statistics
                 (library(doBy))
                 summary_NTI_persampletype <- summaryBy(NTI~Larks+Type,data=NTIdata_persampletype, FUN = function(x) { c(m = mean(x), s = sd(x), n=length(x))})
                 summary_MNTD_persampletype <- summaryBy(c(mntd.obs) ~Larks+Type,data=NTIdata_persampletype, FUN = function(x) { c(m = mean(x), s = sd(x), n=length(x))})
                 
                 summary_NTI_persampletype[["NTI.df"]] <- summary_NTI_persampletype$NTI.n-1
                 
                 # manually find t-value based on n in inverse t table (online)
                 summary_NTI_persampletype[["t"]] <-c(2.228, 2.201, 2.179,2.179,2.201, 2.228, 2.179, 2.228, 2.447, 2.447)
                 
                 # calculate SE 
                 summary_NTI_persampletype[["NTI.se"]] <- summary_NTI_persampletype$NTI.s/sqrt(summary_NTI_persampletype$NTI.n)
                 
                 # calculate CI abs
                 summary_NTI_persampletype[["NTI.CI"]] <- summary_NTI_persampletype$NTI.se*summary_NTI_persampletype$t
                 
                 # calculate CI range based on mean (low)
                 summary_NTI_persampletype[["NTI.lowCI"]] <- summary_NTI_persampletype$NTI.m-summary_NTI_persampletype$NTI.CI
                 
                 # calculate CI range based on mean (high)
                 summary_NTI_persampletype[["NTI.highCI"]] <- summary_NTI_persampletype$NTI.m+summary_NTI_persampletype$NTI.CI
                 
                 
                 summary_addMNTD <- cbind(summary_NTI_persampletype, summary_MNTD_persampletype[,3:5])
                 
                 colnames(summary_addMNTD)[12] <- "mntd.obs.m"
                 colnames(summary_addMNTD)[13] <- "mntd.obs.s"
                 colnames(summary_addMNTD)[14] <- "mntd.obs.n"
                 
                 
                 #calculate mntd SE
                 summary_addMNTD[["mntd.se"]] <- summary_addMNTD$mntd.obs.s/sqrt(summary_addMNTD$mntd.obs.n)
                 
                 # calculate mntd CI abs
                 summary_addMNTD[["mntd.CI"]] <- summary_addMNTD$mntd.se*summary_addMNTD$t
                 
                 # view summary
                 summary_NTI_persampletype
                 
                 # add NTI column to main data
                 NTIdata_persampletype$NTI <- NTIdata_persampletype$mntd.obs.z*-1
                 
                 summary_NTI_persampletype$Type <- factor(summary_NTI_persampletype$Type, levels = c("cloacal", "broodpatch","feathers","nest material", "soil"))
                 
                 pdf("NTI_perType_perSpecies_persample5000.nullmodel4_nonAbundWeighted_new2.pdf", useDingbats = F)
                 ggplot(aes(x=Type, y=NTI.m, fill=Larks), data=summary_NTI_persampletype) +
                   geom_point(data=summary_NTI_persampletype, aes(x=Type, y=NTI.m), position=position_dodge(width=0.75), size=3, color="black") +
                   geom_errorbar(aes(ymin=NTI.lowCI, ymax=NTI.highCI), width=.3, position=position_dodge(width=0.75),color="black") +
                   geom_point(data=NTIdata_persampletype, aes(x=Type, y=NTI, color=Larks), position=position_jitterdodge(), size=3, alpha=0.5) +
                   #geom_jitter(data=NTIdata_no_eggs, aes(x=Type, y=NTI, color=Larks), position=position_jitter(width=.3), size=3, alpha=0.5) +
                   geom_hline(yintercept=2) +
                   geom_hline(yintercept=-2) +
                   scale_y_continuous(limits=c(-4,6))
                 dev.off()      
                 
                 
                 
                 ############## PHYLOGENETIC STRUCTURE - MPD and NRI calculations ########################
                 
                 ### NULL model = independent swap method (null model 4 in Kembel, 2009) ###
                 
                 
                 library(phyloseq)
                 library(picante)
                 
                 #To calculate the NRI:
                 
                 #Now we will calculate NRI (Net Relatedness Index) Rembember, negative NRI values indicate a phylogenetic overdispersion, and positive NRI values indicate phylogenetic clustering.
                 #First we need to make a phylogenetic distance matrix.
                 
                 phydist.cl <- cophenetic(phy.cl)
                 phydist.bp <- cophenetic(phy.bp)
                 phydist.f <- cophenetic(phy.f)
                 phydist.nm <- cophenetic(phy.nm)
                 phydist.s <- cophenetic(phy.s)
                 
                 
                 # Take a look at phydist.  This is a matrix where the rows and columns are the taxa and the elements of the matrix are the phylogenetic distance between those pairs of taxa.  
                 
                 #phydist
                 
                 # To calculate the NRI:
                 
                 #ses.mpd(samporder, phydist,null.model="taxa.labels", abundance.weighted=TRUE)
                 
                 # The rows are the communities.  The first four columns should be pretty straight forward given the definition of the NRI from lecture.  To review:
                 # ntaxa Number of taxa in community
                 # mpd.obs Observed mean pairwise distance (MPD) in community
                 # mpd.rand.mean Mean MPD in null communities
                 # mpd.rand.sd Standard deviation of MPD in null communities
                 # mpd.obs.rank Rank of observed MPD vs. null communities
                 # mpd.obs.z Standardized effect size of MPD vs. null communities (equivalent to -NRI)
                 
                 # The fifth column is the rank of the score against randomized communities; the sixth column is the negative NRI; and the seventh column is the one tailed p-value for significantly high NRI.  So that communities 1, 2 and 3 are significantly clustered and community 5 is significantly spread out.
                 
                 # To calculate the NRI:
                 MPD_NRI_all_cl5000 <- ses.mpd(samporder.cl, phydist.cl,null.model="independentswap",runs = 999, iterations = 1000, abundance.weighted=FALSE)
                 MPD_NRI_all_bp5000 <- ses.mpd(samporder.bp, phydist.bp,null.model="independentswap", runs = 999, iterations = 1000, abundance.weighted=FALSE)
                 MPD_NRI_all_f5000 <- ses.mpd(samporder.f, phydist.f,null.model="independentswap", runs = 999, iterations = 1000, abundance.weighted=FALSE)
                 MPD_NRI_all_nm5000 <- ses.mpd(samporder.nm, phydist.nm,null.model="independentswap", runs = 999, iterations = 1000, abundance.weighted=FALSE)
                 MPD_NRI_all_s5000 <- ses.mpd(samporder.s, phydist.s,null.model="independentswap", runs = 999, iterations = 1000, abundance.weighted=FALSE)
                 
                 # abundance-weighted NRI:
                 #MPD_NRI_all_cl5000 <- ses.mpd(samporder.cl, phydist.cl,null.model="independentswap",runs = 999, iterations = 1000, abundance.weighted=TRUE)
                 #MPD_NRI_all_bp5000 <- ses.mpd(samporder.bp, phydist.bp,null.model="independentswap", runs = 999, iterations = 1000, abundance.weighted=TRUE)
                 #MPD_NRI_all_f5000 <- ses.mpd(samporder.f, phydist.f,null.model="independentswap", runs = 999, iterations = 1000, abundance.weighted=TRUE)
                 #MPD_NRI_all_nm5000 <- ses.mpd(samporder.nm, phydist.nm,null.model="independentswap", runs = 999, iterations = 1000, abundance.weighted=TRUE)
                 #MPD_NRI_all_s5000 <- ses.mpd(samporder.s, phydist.s,null.model="independentswap", runs = 999, iterations = 1000, abundance.weighted=TRUE)
                 
                 # merge results
                 MPD_NRI_all_persample5000 <- rbind(MPD_NRI_all_cl5000,MPD_NRI_all_bp5000,MPD_NRI_all_f5000,MPD_NRI_all_nm5000,MPD_NRI_all_s5000)
                 
                 # transform z-scores in NRI by multiplying by -1
                 MPD_NRI_all_persample5000[["NRI"]] <- MPD_NRI_all_persample5000$mpd.obs.z*-1 
                 
                 # Results stored in write.table(MPD_NRI_all_samples5000, "NRI_MPD_null_vs_obs_allSamples5000_20161219.xlsx", row.names=TRUE,col.names=T,sep="\t")
                 write.table(MPD_NRI_all_persample5000, "NRI_MPD_null_vs_obs_all_persample5000_20170915.xlsx", row.names=TRUE,col.names=T,sep="\t")
                 
                 # add NRI data to metadata samples5000
                 mapping_NRI_persampletype <- sample_data(samples5000)
                 colnames(mapping_NRI_persampletype)[1] <- "X.SampleID"
                 # colnames(mapping_NRI_persampletype)[1] # "sampleID"
                 MPD_NRI_all_persample5000$sampleID <- rownames(MPD_NRI_all_persample5000)
                 nri_persampletype <- as.data.frame(mapping_NRI_persampletype@.Data)
                 colnames(nri_persampletype) <- mapping_NRI_persampletype@names
                 colnames(nri_persampletype)[1] <- "sampleID"
                 nri_persampletype[nri_persampletype$sampleID=="192",7] <- "B"
                 NRIdata_persampletype[NRIdata_persampletype$sampleID=="54",7] <- "B"
                 
                 
                 nrow(nri_persampletype) # 175
                 nrow(MPD_NRI_all_persample5000) # 110
                 
                 
                 # --> copy data from MPD_NRI_all_persample5000 to nri_persampletype for only matiching sampleID values
                 NRIdata_persampletype <- merge(nri_persampletype, MPD_NRI_all_persample5000, "sampleID")
                 NRIdata_persampletype[NRIdata_persampletype$sampleID=="54",7] <- "B"
                 NRIdata_persampletype[NRIdata_persampletype$sampleID=="118",10] <- "B14.11"
                 
                 
                 NRIdata_persampletype$Type <- factor(NRIdata_persampletype$Type, levels = c("cloacal", "broodpatch","feathers","nest material", "soil"))
                 
                 
                 # Analysis of variance of NRI 
                 # lm with Larks x Type interaction
                 
                 lm_NRI_persampletype <- lm(NRI~Larks*Type, data=NRIdata_persampletype)
                 
                 anova(lm_NRI_persampletype)
                 
                 # lm NRI without interaction
                 lm_NRI <- lm(NRI~Larks+Type, data=NRIdata_persampletype)
                 anova(lm_NRI)
                 
                 
                 # create contrast matrix for each level of exp for each level of sampling_point
                 library(multcomp)
                 
                 X <- model.matrix(~ Larks * Type, data = NRIdata_persampletype)
                 
                 Tukey <- contrMat(table(NRIdata_persampletype$Type), "Tukey")
                 K1 <- cbind(Tukey, matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)))
                 rownames(K1) <- paste(levels(NRIdata_persampletype$Larks)[1], rownames(K1), sep = ":")
                 K2 <- cbind(matrix(0, nrow = nrow(Tukey), ncol = ncol(Tukey)), Tukey)
                 rownames(K2) <- paste(levels(NRIdata_persampletype$Larks)[2], rownames(K2), sep = ":")
                 
                 K <- rbind(K1,K2)
                 colnames(K) <- c(colnames(Tukey), colnames(Tukey))  
                 
                 glht_NRI <- summary(glht(lm_NRI_persampletype, linfct = K))
                 glht_NRI
                 
                 
                 
                 ## Test if NRI values differ between species in each sample type
                 
                 # make intfac
                 NRIdata_persampletype$IntFac <- interaction(NRIdata_persampletype$Larks, NRIdata_persampletype$Type, drop=T)
                 
                 # define linear model
                 lm_NRI_persampletype_IntFac <- lm(NRI~IntFac, data=NRIdata_persampletype)
                 
                 # test linear hypothesis of the difference between mean WL and mean SL
                 library(multcomp)
                 summary(glht(lm_NRI_persampletype_IntFac, mcp(IntFac="Tukey")))
                 
                 
                 
                 # Calculate summary statistics
                 (library(doBy))
                 summary_NRI_persampletype <- summaryBy(NRI~Larks+Type,data=NRIdata_persampletype, FUN = function(x) { c(m = mean(x), s = sd(x), n=length(x))})
                 summary_MPD_persampletype <- summaryBy(c(mpd.obs) ~Larks+Type,data=NRIdata_persampletype, FUN = function(x) { c(m = mean(x), s = sd(x), n=length(x))})
                 
                 summary_NRI_persampletype[["NRI.df"]] <- summary_NRI_persampletype$NRI.n-1
                 
                 # manually find t-value based on n in inverse t table (online)
                 summary_NRI_persampletype[["t"]] <-c(2.201, 2.228, 2.179,2.179,2.201, 2.179, 2.228, 2.228, 2.447, 2.447)
                 
                 # calculate SE 
                 summary_NRI_persampletype[["NRI.se"]] <- summary_NRI_persampletype$NRI.s/sqrt(summary_NRI_persampletype$NRI.n)
                 
                 # calculate CI abs
                 summary_NRI_persampletype[["NRI.CI"]] <- summary_NRI_persampletype$NRI.se*summary_NRI_persampletype$t
                 
                 # calculate CI range based on mean (low)
                 summary_NRI_persampletype[["NRI.lowCI"]] <- summary_NRI_persampletype$NRI.m-summary_NRI_persampletype$NRI.CI
                 
                 # calculate CI range based on mean (high)
                 summary_NRI_persampletype[["NRI.highCI"]] <- summary_NRI_persampletype$NRI.m+summary_NRI_persampletype$NRI.CI
                 
                 
                 summary_addMPD <- cbind(summary_NRI_persampletype, summary_MPD_persampletype[,3:5])
                 
                 colnames(summary_addMPD)[12] <- "mpd.obs.m"
                 colnames(summary_addMPD)[13] <- "mpd.obs.s"
                 colnames(summary_addMPD)[14] <- "mpd.obs.n"
                 
                 
                 #calculate mpd SE
                 summary_addMPD[["mpd.se"]] <- summary_addMPD$mpd.obs.s/sqrt(summary_addMPD$mpd.obs.n)
                 
                 # calculate mpd CI abs
                 summary_addMPD[["mpd.CI"]] <- summary_addMPD$mpd.se*summary_addMPD$t
                 
                 
                 summary_NRI_persampletype
                 
                 
                 # add NRI column to main data
                 NRIdata_persampletype$NRI <- NRIdata_persampletype$mpd.obs.z*-1
                 
                 summary_NRI_persampletype$Type <- factor(summary_NRI_persampletype$Type, levels = c("cloacal", "broodpatch","feathers","nest material", "soil"))
                 
                 pdf("NRI_perType_perSpecies_persample5000.nullmodel4.pdf", useDingbats = F)
                 ggplot(aes(x=Type, y=NRI.m, fill=Larks), data=summary_NRI_persampletype) +
                   geom_point(data=summary_NRI_persampletype, aes(x=Type, y=NRI.m), position=position_dodge(width=0.75), size=3, color="black") +
                   geom_errorbar(aes(ymin=NRI.lowCI, ymax=NRI.highCI), width=.3, position=position_dodge(width=0.75),color="black") +
                   geom_point(data=NRIdata_persampletype, aes(x=Type, y=NRI, color=Larks), position=position_jitterdodge(), size=3, alpha=0.5) +
                   #geom_jitter(data=NRIdata_no_eggs, aes(x=Type, y=NRI, color=Larks), position=position_jitter(width=.3), size=3, alpha=0.5) +
                   geom_hline(yintercept=2) +
                   geom_hline(yintercept=-2) +
                   scale_y_continuous(limits=c(-4,6))
                 dev.off()  
                 
                 
                 
                 
                 