# selection analysis


#install.packages("mgcv")
#install.packages("mvtnorm")
### install from source: https://cran.r-project.org/src/contrib/Archive/gsg/
#install.packages("~/Downloads/gsg_2.0.tar.gz", repos = NULL, type="source")
#install.packages("AER")

library("mgcv")
library("mvtnorm")
library("gsg")
library("AER")
library("tidyverse")


#### Estimate linear and quadratic selection gradients ####

obj <- load("data/abundance_vs_phenotype.rda")
d <- abundance_vs_phenotype
hn <- subset(d, nitrogen %in% "HN")
ln <- subset(d, nitrogen %in% "LN")
dt <- as.data.frame(table(hn$ASV))
# plant traits
pt <- as.data.frame(table(d$phenotype)) #pt
mt <- as.data.frame(table(d$ASV)) #150


# Before running the model, microbiome data are mean-centered and variance standardized.

### calculate selection gradients for HN

out <- data.frame()
for(i in 1:nrow(mt)){
  sub <- subset(hn, ASV %in% mt$Var1[i] & phenotype %in% "CC_Aug12")
  m1 <- gam(blup_logrel ~ s(value), data=sub)
  
  fit <- gam.gradients(mod=m1, phenotype="value", standardized =F, se.method = "boot.case", n.boot=1000,
                       refit.smooth = T)
  
  temp <- fit$ests
  temp$type <- c("B", "G")
  temp$ASV <- as.character(mt$Var1[i])
  out <- rbind(out, temp)
}
write.table(out, file="cache/hn_selection_gredients.csv", sep=",", row.names = FALSE, quote=FALSE)
out1 <- subset(out, estimates != 0)
par(mfrow=c(2,2)) 
gam.check(m1, pch=19, cex=.3)


### calculate selection gradients for LN

out <- data.frame()
for(i in 1:nrow(mt)){
  sub <- subset(ln, ASV %in% mt$Var1[i] & phenotype %in% "CC_Aug12")
  m1 <- gam(blup_logrel ~ s(value), data=sub)
  
  fit <- gam.gradients(mod=m1, phenotype="value", standardized =F, se.method = "boot.case", n.boot=1000,
                       refit.smooth = T)
  
  temp <- fit$ests
  temp$type <- c("B", "G")
  temp$ASV <- as.character(mt$Var1[i])
  out <- rbind(out, temp)
}
out2 <- subset(out, estimates != 0)
write.table(out, file="cache/ln_selection_gredients.csv", sep=",", row.names = FALSE, quote=FALSE)


# save selection gradient data

hn <- read_csv("cache/hn_selection_gredients.csv")
hn1 <- subset(hn, estimates != 0 & P.value < 0.05)
hn1$treatment <- "HN"
ln <- read_csv("cache/ln_selection_gredients.csv")
ln1 <- subset(ln, estimates != 0 & P.value < 0.05)
ln1$treatment <- "LN"
df <- rbind(hn1, ln1)
df1 <- subset(df, type == "B")
write.table(df1, "data/results_linear_selection_gredients.csv", sep=",", row.names = FALSE, quote=FALSE)
df2 <- subset(df, type == "G")
write.table(df2, "data/quadratic_selection_gredients.csv", sep=",", row.names = FALSE, quote=FALSE)


# plot selection gradient data


sg_linear <- read_csv("data/results_linear_selection_gredients.csv")


sg_quadratic <- read_csv("data/quadratic_selection_gredients.csv")
# only 4 tax groups have significant quadratic selection gradients:
# Blastococcus, Pseudomonas umsongensis, Chthoniobacter flavus, Luteolibacter pohnpeiensis


load("data/group_data.rda")

sg_data <- sg_linear %>%
  #mutate(model = "linear") %>%
  #bind_rows(mutate(sg_quadratic, model = "quadratic")) %>%
  left_join(select(group_data, ASV, tax_group)) %>%
  arrange(-estimates) %>%
  mutate(tax_group = factor(tax_group, levels = unique(tax_group)))

p <- ggplot(sg_data, aes(x=tax_group, y=estimates, color= treatment)) +
  geom_point(size=2, alpha = 0.7, position=position_dodge(width=0.5)) +
  geom_errorbar(width=0, size=1, alpha = 0.7, position=position_dodge(width=0.5), aes(ymin=estimates-SE, ymax=estimates+SE)) +
  # + scale_color_manual("Treatment", breaks=c(1:32), values=c(rep("red",6),rep("blue",7), rep("black",13), rep("purple",6)))+
  xlab("")+
  ylab("Linear selection grient") +
  #facet_wrap(~model, ncol = 2) +
  coord_flip() +
  scale_color_manual(values = c("LN"="#C00001", "HN"="#000080")) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_bw() +
  theme(#panel.grid.major = element_blank(),
        axis.text.x = element_text(color="black", angle=90, hjust=1, vjust=0.5), panel.grid.minor = element_blank()) +
  theme(text = element_text(size = 10)) 


p





## get tax groups with significant selection gradients
sg_taxa_HN <- sg_data %>%
  filter(treatment == "HN") %>%
  select(tax_group) %>%
  unique() # 28
sg_taxa_LN <- sg_data %>%
  filter(treatment == "LN") %>%
  select(tax_group) %>%
  unique() # 46







#### GCTB analysis ####


## plot this for 58 tax groups that have significant selection gradient:
sg_taxa <- unique(as.character(sg_data$tax_group))
sg_taxa_ids <- unique(as.character(sg_data$ASV))



load("cache/gctb_selection_values.rda")


# rhizobiome traits

load("data/group_data.rda")

rhizobiome_traits <- selection_values %>%
  filter(type %in% c("microbe_from_Mike")) %>%
  rename(ASV = trait) %>%
  left_join(dplyr::select(group_data, ASV, tax_group)) %>%
  filter((sd1 >0 & sd2 >0) | (sd1 <0 &sd2 <0))

#save(rhizobiome_traits, file="cache/gctb_rhizobiome_traits.rda")

plot_rhizo_traits <- selection_values %>%
  filter(trait %in% rhizobiome_traits$ASV) %>%
  rename(ASV = trait) %>%
  left_join(dplyr::select(group_data, ASV, tax_group)) %>%
  filter(tax_group %in% sg_taxa)


ggplot(plot_rhizo_traits, aes(x=reorder(tax_group, Mean), y=Mean, color=treatment)) +
  #facet_wrap(~treatment, ncol = 2) +
  geom_point(position=position_dodge(width=0.5), size=2) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), 
                width = 0.5,
                position=position_dodge(width=0.5)) +
  coord_flip() +
  scale_color_manual(values = c("LN"="#C00001", "HN"="#000080")) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text=element_text(size=10, face="bold"),
        strip.text.x = element_text(size = 10, face = "bold")
  )




# plant traits

plant_traits <- selection_values %>%
  filter(type %in% c("leaf_nutrient_traits_Field")) %>%
  filter((sd1 >0 & sd2 >0) | (sd1 <0 &sd2 <0))

plot_traits <- selection_values %>%
  filter(trait %in% plant_traits$trait)

ggplot(plot_traits, aes(x=reorder(trait, Mean), y=Mean, color=treatment)) +
  #facet_wrap(~treatment, ncol = 2) +
  geom_point(position=position_dodge(width=0.5), size=2) +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), 
                width = 0.5,
                position=position_dodge(width=0.5)) +
  coord_flip() +
  scale_color_manual(values = c("LN"="#C00001", "HN"="#000080")) +
  ylim(-3.2, 2) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text=element_text(size=10, face="bold"),
        strip.text.x = element_text(size = 10, face = "bold")
  )




# plot nonzero SNPs

df <- read_csv("cache/gtcb_HN_LN_1873traits.csv")

nz <- df %>%
  filter(par %in% "NnzSnp" & type %in% c("leaf_nutrient_traits_Field", "microbe_from_Mike", "Aerial_trait")) %>%
  mutate(type = ifelse(type == "microbe_from_Mike", "Rhizobiome Traits", "Plant Traits")) %>%
  filter(type == "Plant Traits" | trait %in% sg_taxa_ids)


ggplot(nz, aes(x=treatment, y=log10(Mean), color=treatment)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitterdodge(), alpha = 0.5) +
  #scale_color_manual(values = c("Rhizobiome Traits"="#333333", "Plant Traits"="#333333")) +
  scale_color_manual(values = c("HN"="#000080", "LN"="#C00001")) +
  facet_wrap(~type, ncol = 2) +
  ylab("Log10(Nonzero SNPs)") +
  #geom_hline(yintercept=0, linetype="dashed", color = "red") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=10, face="bold"),
        axis.text=element_text(size=10, face="bold"),
        #strip.text.x = element_text(size = 16, face = "bold"),
        legend.position = "none"
  )



# plot effect size vs minor allele frequency
# for examples Bacillus fumarioli and f_Comamonadaceae Unknown Genus


load("data/examples_MAF_vs_EffectSize.rda")

ggplot(examples_MAF_vs_EffectSize, aes(x = A1Frq, y = A1Effect, color = treatment)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("HN"="#000080", "LN"="#C00001")) +
  facet_wrap( ~tax_group, nrow = 1, scales = "free_y") +
  ylab("Minor Allele Effect") +
  xlab("Minor Allele Frequency") +
  theme_bw() +
  theme(legend.position = "none")




