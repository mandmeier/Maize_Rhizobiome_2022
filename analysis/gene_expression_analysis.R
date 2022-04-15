#gene_expression_analysis


library("tidyverse")
library("ggpubr")





#### find genes associated with 613 MAPLs (all bins with sign. SNPs) ####

# import maize reference gene data from maizegdb v4
load("data/B73_RefGen_v4_gene_ids.rda")
load("data/B73_RefGen_v4_gene_pos.rda")

load("cache/all_gwas_signals.rda")

shortlist <- all_gwas_signals %>%
  ungroup() %>%
  dplyr::select(nitrogen, chr, bin) %>%
  unique() %>%
  mutate(ext_bin_start = (bin-2)*10000+1, ext_bin_end = (bin+1)*10000)

LN <- filter(shortlist, nitrogen == "LN")
HN <- filter(shortlist, nitrogen == "HN")


### find all genes in reference genome that are within +/- 10kb of the
## 10kb regions with significant GWAS signals

shortlist$genes <- NA
gene_list <- c()
for(i in c(1:nrow(shortlist))){
  #i <- 1
  chr <- shortlist$chr[i]
  print(paste("MAPL",i))
  range_low <- (shortlist$bin[i]*10000) - 100000
  range_hi <- (shortlist$bin[i]*10000) + 100000
  gene_window <- gene_pos %>%
    filter(seqid == chr) %>%
    filter(start > range_low & end < range_hi)
  
  if (nrow(gene_window) > 0){
    genes <- c()
    for (gene in c(1:nrow(gene_window))){
      #gene <- 1
      overlap <- sum(c(gene_window$start[gene]:gene_window$end[gene]) %in% c(shortlist$ext_bin_start[i]:shortlist$ext_bin_end[i]))
      #print(overlap)
      if(overlap > 0){
        #print("gene is near bin")
        gene <- str_match(gene_window$attributes[gene], "gene_id=(.+);")[2]
        #genes <- c(genes, str_match(gene_window$attributes[gene], "gene_id=(.+);")[2])
        row <- shortlist[i,]
        row$genes <- gene
        print(gene)
        gene_list <- rbind(gene_list, row)
      }
    }
    shortlist$genes[i] <- ifelse(is.null(genes), NA, paste(genes, collapse = ';'))
  }
}



MAPL_genes <- gene_list


#save(MAPL_genes, file = "cache/MAPL_genes_613.rda")


#### find genes associated with 119 MAPLs (at least 2 sign. SNPs per bin) ####

# import maize reference gene data from maizegdb v4
load("data/B73_RefGen_v4_gene_ids.rda")
load("data/B73_RefGen_v4_gene_pos.rda")

load("cache/mapl_gwas_signals.rda")

shortlist <- mapl_gwas_signals %>%
  ungroup() %>%
  dplyr::select(nitrogen, chr, bin) %>%
  unique() %>%
  mutate(ext_bin_start = (bin-2)*10000+1, ext_bin_end = (bin+1)*10000)

LN <- filter(shortlist, nitrogen == "LN")
HN <- filter(shortlist, nitrogen == "HN")


### find all genes in reference genome that are within +/- 10kb of the
## 10kb regions with significant GWAS signals

shortlist$genes <- NA
gene_list <- c()
for(i in c(1:nrow(shortlist))){
  #i <- 1
  chr <- shortlist$chr[i]
  print(paste("MAPL",i))
  range_low <- (shortlist$bin[i]*10000) - 100000
  range_hi <- (shortlist$bin[i]*10000) + 100000
  gene_window <- gene_pos %>%
    filter(seqid == chr) %>%
    filter(start > range_low & end < range_hi)
  
  if (nrow(gene_window) > 0){
    genes <- c()
    for (gene in c(1:nrow(gene_window))){
      #gene <- 1
      overlap <- sum(c(gene_window$start[gene]:gene_window$end[gene]) %in% c(shortlist$ext_bin_start[i]:shortlist$ext_bin_end[i]))
      #print(overlap)
      if(overlap > 0){
        #print("gene is near bin")
        gene <- str_match(gene_window$attributes[gene], "gene_id=(.+);")[2]
        #genes <- c(genes, str_match(gene_window$attributes[gene], "gene_id=(.+);")[2])
        row <- shortlist[i,]
        row$genes <- gene
        print(gene)
        gene_list <- rbind(gene_list, row)
      }
    }
    shortlist$genes[i] <- ifelse(is.null(genes), NA, paste(genes, collapse = ';'))
  }
}



MAPL_genes <- gene_list


save(MAPL_genes, file = "cache/MAPL_genes_119.rda")


unique(MAPL_genes$genes)



#### compare MAPL genes vs other genes in Kremling 2018 data ####


## translate Zm gene names into AC gene names
gene_names <- read_csv("data/V3_V4.csv", col_names = FALSE)

gene_names <- gene_names %>%
  rename(gene_AC=X1, gene=X2) %>%
  dplyr::select(gene, gene_AC)

MAPL_genes <- MAPL_genes %>%
  rename(gene=genes) %>%
  left_join(gene_names)

unique(MAPL_genes$gene_AC)


## 119 MAPLs, 97 Genes within +- 10kb, 73 are in Kremling RNA seq dataset




#load gene expression data from Kremling et.al, 2018
rnadat <- read_delim(file="largedata/RNAseq/rna_seq_kremling.txt", delim="\t")




## these MAPL genes have RNA seq data available
MAPL_genes_rna <- unique(MAPL_genes$gene_AC)[unique(MAPL_genes$gene_AC) %in% colnames(rnadat)]


## select MAPL gene columns from rna dat


MAPL_rna <- rnadat %>%
  dplyr::select(TissueWODate,RawPhenotypeNames, MAPL_genes_rna) %>%
  rename(tissue = TissueWODate, genotype = RawPhenotypeNames) %>%
  pivot_longer(-c(tissue, genotype), names_to = "gene_AC", values_to = "expression") %>%
  left_join(gene_names) %>%
  group_by(tissue, genotype) %>%
  summarize(mean_expression_MAPL = mean(expression))



all_colnames <- data.frame("colname" = colnames(rnadat))


B73_genes_in_data <- all_colnames %>%
  filter(colname %in% gene_names$gene_AC)

#save(B73_genes_in_data, file="cache/B73_genes_in_data.rda")



not_mapl <- gene_names %>%
  filter(!(gene %in% unique(MAPL_genes$gene)))


## calculate mean gene expression by tissue type
ALL_rna <- rnadat %>%
  filter(Tissue != "L3Mid") %>%
  dplyr::select(TissueWODate,RawPhenotypeNames, B73_genes_in_data$colname) %>%
  #dplyr::select(TissueWODate,RawPhenotypeNames, AC148152.3_FG001:GRMZM2G009080) %>%
  #dplyr::select(TissueWODate,RawPhenotypeNames, AC148152.3_FG001:GRMZM2G123843) %>%
  rename(tissue = TissueWODate, genotype = RawPhenotypeNames) %>%
  pivot_longer(-c(tissue, genotype), names_to = "gene_AC", values_to = "expression") %>%
  left_join(gene_names) %>%
  mutate(group=ifelse(gene %in% MAPL_genes$gene,"MAPL", "other")) 



unique(filter(ALL_rna, group == "MAPL")$gene) # 73 MAPL genes
unique(filter(ALL_rna, group == "other")$gene) # 29771 other genes 


genes_in_7_tissue_data <- unique(ALL_rna$gene)

#save(genes_in_7_tissue_data, file = "cache/genes_in_7_tissue_data.rda")


mean_gene_expression <- ALL_rna %>%
  group_by(tissue, genotype, group) %>%
  summarize(mean_expression = mean(expression), sd_expression = sd(expression))


#save(mean_gene_expression, file = "data/mean_gene_expression_by_tissue_type_kremling.rda")
load("data/mean_gene_expression_by_tissue_type_kremling.rda")


mean_gene_expression$tissue <- factor(mean_gene_expression$tissue, levels = c("GRoot", "GShoot", "L3Base", "L3Tip", "LMAD", "LMAN", "Kern"))

colors <- c(MAPL="#9cd6ff", other="#fff2bd")

ggplot(mean_gene_expression, aes(x=tissue, y=mean_expression, fill=group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(group = group), label = "p.signif") +
  ylab("mean gene expression [FPKM]") +
  ggtitle("mean gene expression by tissue type") +
  scale_fill_manual(name = "", labels = c("MAPL genes", "other genes"), values = colors) +
  #scale_fill_manual(values = colors) +
  theme_bw() +
  theme(axis.title.x = element_blank(), legend.position = c(0.85, 0.7))












#### compare MAPL genes vs other genes in 4 genotypes (growth chamber) ####


gendat <- read_csv("largedata/rna_seq_4genotypes.csv")


load("cache/MAPL_genes_119.rda")
load("cache/genes_in_7_tissue_data.rda")

MAPL_genes_zm <- unique(MAPL_genes$genes)


genes_in_gen_data <- unique(c(MAPL_genes_zm, genes_in_7_tissue_data))

colnames(gendat)[1] <- "X1"

all_Zm <- sort(unique(gendat$X1[startsWith(gendat$X1, "Zm00001")]))

test <- data.frame("gen" = all_Zm)

rna_gen <- gendat %>%
  rename(gene = X1) %>%
  pivot_longer(-gene, names_to = "id", values_to = "expression") %>%
  filter(gene %in% all_Zm) %>%
  separate(id, c("nitrogen", "genotype", "tissue", "rep"), "_") %>%
  mutate(group = ifelse(gene %in% MAPL_genes_zm, "MAPL", "other"))

mean_gene_expression <- rna_gen %>%
  group_by(nitrogen, genotype, tissue, group) %>%
  summarize(mean_expression = mean(expression), sd_expression = sd(expression))


unique(filter(rna_gen, group == "MAPL")$gene) # 97 MAPL genes
unique(filter(rna_gen, group == "other")$gene) # 44049 other genes



#save(mean_gene_expression, file = "data/mean_gene_expression_by_tissue_type_4genotypes.rda")

load("data/mean_gene_expression_by_tissue_type_4genotypes.rda")

mean_gene_expression$tissue <- factor(mean_gene_expression$tissue, levels = c("Root", "leaf"))


colors <- c(MAPL="#9cd6ff", other="#fff2bd")

ggplot(mean_gene_expression, aes(x=tissue, y=mean_expression, fill=group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(group = group), label = "p.signif") +
  facet_wrap(~nitrogen, nrow = 1) +
  ylab("mean gene expression [FPKM]") +
  ggtitle("mean gene expression by tissue type") +
  scale_fill_manual(name = "", labels = c("MAPL genes", "other genes"), values = colors) +
  #scale_fill_manual(values = colors) +
  theme_bw() +
  theme(axis.title.x = element_blank())

