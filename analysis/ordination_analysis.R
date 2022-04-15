
library("tidyverse")
library("phyloseq")
library("vegan")
library("FactoMineR")
library("factoextra")



#### plot Shannon diversity between 4632 (ps_common) and 3728 sets (ps_bothyears) ####

load(file = "data/ps_common.rda")

rch_4632 <- estimate_richness(ps_common, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "Fisher"))
rch_4632$asv_set <- "all_ASVs"
median(rch_4632$Shannon) # 6.434921

load(file = "data/ps_bothyears.rda")

rch_3728 <- estimate_richness(ps_bothyears, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "Fisher"))
rch_3728$asv_set <- "ASVs_common_to_both_years"
median(rch_3728$Shannon) # 6.287561


diversity <- rbind(rch_4632, rch_3728)


colors <- c("#743282", "#dc9200")

shannon_plot <- ggplot(diversity, aes(x = asv_set, y = Shannon)) +
  geom_jitter(aes(color = asv_set), size=1, width=0.2) +
  geom_boxplot(outlier.shape = NA) +
  #stat_compare_means(aes(group = asv_set), label = "p.format") +
  scale_color_manual(values=colors) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

shannon_plot



#### Constrained Ordination, 4632 ASVs ####



## https://rdrr.io/rforge/vegan/man/adonis.html
# convert to relative abundance, aka Total Sum Normalization [TSS] and log transform ASV counts
ps <- transform_sample_counts(ps_common, function(x) x / sum(x) + 0.001)
sdat <- data.frame(sample_data(ps))
# exclude checks
ps <- subset_samples(ps, !genotype %in% c("CHECK", "CHECK2"))


## get distance matrix: this takes a long time!
#dm <- distance(ps, method = "wunifrac")

### replace missing values with median
#dm[is.na(dm)] <- median(dm, na.rm = TRUE)
#save(dm, file='largedata/distance_matrix_wunifrac.rda')

load("largedata/distance_matrix_wunifrac.rda")


## calculate bray-custis distance matrix
#bray.cap.whole <- capscale(as.dist(dm) ~ year + genotype + nitrogen + block + sp + spb, data = sdat, add = T, na.action = na.exclude)

load("largedata/bray.cap.whole.rda")


caps <- rownames_to_column(data.frame(scores(bray.cap.whole)$sites), var = "Sample_ID")
bray.cap.whole.axes <- left_join(sdat, caps)

## find genotypes common to both years, without check

y2018 <- sdat %>%
  filter(year == "Y2018")
y2019 <- sdat %>%
  filter(year == "Y2019")
common_genotypes <- Reduce(intersect, list(unique(y2018$genotype),unique(y2019$genotype)))
common_genotypes <- common_genotypes[-1]

bray.cap.whole.axes$color_groups <- ifelse(bray.cap.whole.axes$genotype %in% common_genotypes, "common_HN", bray.cap.whole.axes$nitrogen)
bray.cap.whole.axes$color_groups <- ifelse(bray.cap.whole.axes$genotype %in% common_genotypes & bray.cap.whole.axes$nitrogen == "-N", "common_LN", bray.cap.whole.axes$color_groups)

bray.cap.whole.axes$shape_groups <- ifelse(bray.cap.whole.axes$genotype %in% common_genotypes, "common", "unique")
bray.cap.whole.axes$shape_groups <- factor(bray.cap.whole.axes$shape_groups, levels = c("unique", "common"))

colors <- c("-N"="#ff0000","common_LN"="#C00001", "+N"="#0000ff", "common_HN"="#000080")
colors = c("#ffc30b", "#000000")


ord_plot <- ggplot(bray.cap.whole.axes, aes(x = CAP1, y = CAP2, color = shape_groups, shape = nitrogen)) +
  geom_point(size=1, alpha = 0.5) +
  facet_wrap(~ year, nrow = 1) +
  scale_color_manual(values=colors) +
  labs(x = "Constrained PCo1 (31.80%)", y = "Constrained PCo2 (26.24%)") +
  theme_bw()

ord_plot

# Permutation ANOVA (this takes 24h or so...) 
#permanova <- adonis(dm~year + genotype + nitrogen + block + sp + spb,data = sdat,add = T)
#save(permanova, file = "largedata/permanova.rda")
#load("largedata/permanova.rda")




#### Constrained Ordination, 3626 ASVs ####

load("data/ps_asv.rda")

# convert to relative abundance, aka Total Sum Normalization [TSS] and log transform ASV counts
ps <- transform_sample_counts(ps_asv, function(x) log(x / sum(x) + 0.001))
# exclude checks
ps <- subset_samples(ps, !genotype %in% c("CHECK", "CHECK2"))
sdat <- data.frame(sample_data(ps))

## get distance matrix: this takes a long time!
#dm <- distance(ps, method = "wunifrac")

### replace missing values with median
#dm[is.na(dm)] <- median(dm, na.rm = TRUE)
#save(dm, file='largedata/distance_matrix_wunifrac_3626.rda')

load("largedata/distance_matrix_wunifrac_3626.rda")

## calculate bray-curtis distance matrix

#bray.cap.whole <- capscale(as.dist(dm) ~ year + genotype + nitrogen + block + sp + spb, data = sdat, add = T, na.action = na.exclude)

#save(bray.cap.whole, file='largedata/bray.cap.whole_3626.rda')

load("largedata/bray.cap.whole_3626.rda")


caps <- rownames_to_column(data.frame(scores(bray.cap.whole)$sites), var = "Sample_ID")

bray.cap.whole.axes <- left_join(sdat, caps)


colors <- c("#c00001", "#000080")

ord_plot <- ggplot(filter(bray.cap.whole.axes, year=="Y2019"), aes(x = CAP1, y = CAP2, color = nitrogen)) +
  geom_point(size=2, alpha = 0.5, shape=18) +
  scale_color_manual(values=colors) +
  labs(x = "Constrained PCo1 (31.80%)", y = "Constrained PCo2 (26.24%)") +
  theme_bw()

ord_plot





#### PCA for high-level rhizobiome traits, using ASV BLUPs ####





## LN PCA
asvs_logrel <- read_csv("data/blup_lowN_3618_asvs.csv")

asv_table <- asvs_logrel %>%
  #replace(is.na(.), 0) %>%
  column_to_rownames(var = "genotype")


res_pca_LN <- PCA(asv_table, scale.unit = FALSE, ncp = 10, graph = FALSE)
res_pca_LN


#sdat <- read_csv("data/BG_sample_data.csv") 

#subpop <- sdat %>%
#  dplyr::select(genotype, subpopulation) %>%
#  unique()


coords <- data.frame(res_pca_LN$ind$coord)

plot_data_LN <- coords %>%
  rownames_to_column(var = "genotype") %>%
  mutate(nitrogen = "LN")
#left_join(subpop, by = "genotype")

xlab <- paste0("PC1 [", round(res_pca_LN$eig[1, c("percentage of variance")], 1),"%]" )
ylab <- paste0("PC2 [", round(res_pca_LN$eig[2, c("percentage of variance")], 1),"%]" )

pca_plot <- ggplot(plot_data_LN, aes(x = Dim.1, y = Dim.2)) + # color = subpopulation
  geom_point(alpha = 0.5) +
  xlab(xlab) +
  ylab(ylab) +
  theme_bw()

pca_plot


PCs_LN <- round(res_pca_LN$eig[1:10, c("percentage of variance")], 2)



## HN PCA
asvs_logrel <- read_csv("data/blup_stdN_3618_asvs.csv")

asv_table <- asvs_logrel %>%
  #replace(is.na(.), 0) %>%
  column_to_rownames(var = "genotype")


res_pca_HN <- PCA(asv_table, scale.unit = FALSE, ncp = 10, graph = FALSE)
res_pca_HN



coords <- data.frame(res_pca_HN$ind$coord)

plot_data_HN <- coords %>%
  rownames_to_column(var = "genotype") %>%
  mutate(nitrogen = "HN")
#left_join(subpop, by = "genotype")

xlab <- paste0("PC1 [", round(res_pca_HN$eig[1, c("percentage of variance")], 1),"%]" )
ylab <- paste0("PC2 [", round(res_pca_HN$eig[2, c("percentage of variance")], 1),"%]" )

pca_plot <- ggplot(plot_data_HN, aes(x = Dim.1, y = Dim.2)) + # color = subpopulation
  geom_point(alpha = 0.5) +
  xlab(xlab) +
  ylab(ylab) +
  theme_bw()

pca_plot


PCs_HN <- round(res_pca_HN$eig[1:10, c("percentage of variance")], 2)

## combine datasets for plots

PCs <- data.frame("PC" = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), "HN" = PCs_HN, "LN" = PCs_LN)

PCs$PC <- factor(PCs$PC, levels = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"))

## top 10 PCs, -N

ggplot(PCs, aes(x = PC, y = LN)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label= paste0(round(LN, 1), "%")), size = 3 , vjust=-0.5) +
  scale_y_continuous(limits = c(0, 25)) +
  ylab("% variance explained") +
  theme_bw() +
  theme(axis.title.x = element_blank())


## top 10 PCs, +N
ggplot(PCs, aes(x = PC, y = HN)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label= paste0(round(HN, 1), "%")), size = 3 , vjust=-0.5) +
  scale_y_continuous(limits = c(0, 25)) +
  ylab("% variance explained") +
  theme_bw() +
  theme(axis.title.x = element_blank())





load("data/asv_taxonomy.rda")


## get ASV contributions LN
asv_contributions_LN <- res_pca_LN$var$contrib

## plot top 20 ASV contributions using library factoextra

fviz_contrib(res_pca_LN, choice = "var", axes = 1, top = 20, ggtheme = theme_classic()) +
  ggtitle("Top 20 ASVs contributing to PC1")

contrib <- data.frame(asv_contributions_LN) %>%
  rownames_to_column("ASV") %>%
  rename(PC1 = Dim.1) %>%
  dplyr::select(ASV, PC1) %>%
  left_join(dplyr::select(asv_taxonomy, ASV, tax_group), by = "ASV") %>%
  group_by(tax_group) %>%
  summarize( contrib_to_PC1 = sum(PC1)) %>%
  arrange(-contrib_to_PC1) %>%
  mutate(tax_group = ifelse(contrib_to_PC1 > 1, tax_group, "other")) %>%
  group_by(tax_group) %>%
  summarize( contrib_to_PC1 = sum(contrib_to_PC1)) %>%
  arrange(-contrib_to_PC1)


lbls <- paste0(round(contrib$contrib_to_PC1, 2), "% ",contrib$tax_group)
colors = c("Ralstonia pickettii" = "blue", "Burkholderia oklahomensis" = "lightblue", "other" = "grey", "Sphingobium herbicidovorans 1" = "red", "Dyella jiangningensis" = "orange") 


pie(contrib$contrib_to_PC1 ,labels = lbls, col=colors,
    main="Variance contribution to PC1")




## get ASV contributions HN

asv_contributions_HN <- res_pca_HN$var$contrib

fviz_contrib(res_pca, choice = "var", axes = 2, top = 5, ggtheme = theme_classic()) +
  ggtitle("Top 5 ASVs contributing to PC2")

contrib <- data.frame(asv_contributions_HN) %>%
  rownames_to_column("ASV") %>%
  rename(PC1 = Dim.1) %>%
  dplyr::select(ASV, PC1) %>%
  left_join(dplyr::select(asv_taxonomy, ASV, tax_group), by = "ASV") %>%
  group_by(tax_group) %>%
  summarize( contrib_to_PC1 = sum(PC1)) %>%
  arrange(-contrib_to_PC1) %>%
  mutate(tax_group = ifelse(contrib_to_PC1 > 1, tax_group, "other")) %>%
  group_by(tax_group) %>%
  summarize( contrib_to_PC1 = sum(contrib_to_PC1)) %>%
  arrange(-contrib_to_PC1)

contrib$tax_group


lbls <- paste0(round(contrib$contrib_to_PC1, 2), "% ",contrib$tax_group)
colors = c("Massilia putida" = "#ffc500",
           "Enterobacter"   = "brown",
           "other" = "grey",
           "Ralstonia pickettii" = "blue",
           "Labrys miyagiensis" = "darkgreen",
           "Sphingobium herbicidovorans 1" = "red",
           "Acinetobacter nosocomialis"  = "purple",
           "Kribbella karoonensis" = "#00cdac",
           "Sphingobacterium siyangense-multivorum" = "lightgreen") 


pie(contrib$contrib_to_PC1 ,labels = lbls, col=colors,
    main="Variance contribution to PC1")





