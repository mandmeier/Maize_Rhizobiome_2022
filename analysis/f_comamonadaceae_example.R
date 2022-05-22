# f_comamonadaceae_example, -N treatment

library("tidyverse")
library("LDheatmap")
library("genetics")
library("ggpubr")
library("ggpmisc")



#### prepare plot object for GWAS manhattan plot  ####

load("data/group_data.rda")

tgr <- "f_Comamonadaceae Unknown Genus"
asv <- "asv_000260"
gwas_dat <- read_delim("largedata/GWAS/LN/asv_000260.assoc.txt", delim = "\t")
nit <- "LN"


chr_len_10k <- c(308460000, 243680000, 238020000, 250340000, 226360000, 181360000, 185810000, 182420000, 163010000, 152440000)


nCHR <- 10
chindex <- c(0) ## this number needs to be added to position within each chromosome
for (i in c(1:nCHR)){
  #print(i)
  chindex[i+1] <- sum(chr_len_10k[1:i])
}
chindex <- chindex[1:10]
bpadd <- data.frame("chr" = c(1:nCHR), "chindex" = chindex)

options(scipen = 999)


threshold_line <- 7.187346

taxgrp <- tgr
nitr <- nit
title <- paste(taxgrp, "|", nitr)


# add bins
plot_obj <- gwas_dat %>%
  mutate(log10p = -log10(p_wald)) %>%
  left_join(bpadd, by = "chr") %>%
  mutate(psabs = chindex + ps) %>%
  mutate(ps = as.character(ps)) %>%
  mutate(bin = as.numeric(gsub('.{4}$', '', ps)) + 1) %>% # absolute bin position in 10kb units
  group_by(chr, bin) %>%
  add_count(count = sum(-log10(p_wald) >= threshold_line), name="total_sign_snps") %>%
  mutate(ASV=asv) %>%
  left_join(dplyr::select(group_data, ASV, tax_group)) %>%
  mutate(chr = as.character(chr))




## save plot object to plot LD plot later
#save(plot_obj, file="largedata/GWAS/Manhattan_plot_objects/f_Comamonadaceae_Unknown_Genus_LN.rda")


## find position for labels
axis_set <- plot_obj %>% 
  group_by(chr) %>% 
  summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)


## find range for y axis
ylim <- floor(max(plot_obj$log10p, na.rm = TRUE)) + 1


colors = c("1" = "#276FBF", "3" = "#276FBF", "5" = "#276FBF", "7" = "#276FBF", "9" = "#276FBF",
           "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059")

  
  
#### plot full manhattan plot ####
  
  

#test <- plot_obj %>%
#  filter(chr =="3")


manhplot <- ggplot(plot_obj, aes(x = psabs, y = log10p, color=chr)) +
  geom_point(alpha = 0.75) +
  #geom_hline(yintercept = -log10(5), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  #scale_x_continuous(limits = c(xmin, xmax), breaks = breaks, labels=labels) +
  scale_y_continuous(expand = c(0,0), limits = c(1, ylim), breaks = c(1:ylim)) +
  geom_hline(yintercept=threshold_line, linetype="dashed", color = "red") +
  #geom_vline(xintercept = vlines) +
  scale_color_manual(values = colors) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "-log10(p)") + 
  ggtitle(title) +
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )

manhplot

  

#### zoom into GWAS peak ####

  
## define range and labels
chdx <- chindex[10]
xmin <- 23800000 + chdx
xmax <- 24100000 + chdx
breaks <- seq(from=xmin, to=xmax, by=10000)
labels <- seq(xmin-chdx, xmax-chdx, 10000 ) /1000


manhplot_zoom <- ggplot(plot_obj, aes(x = psabs, y = log10p, color=chr)) +
  geom_point(alpha = 0.75) +
  scale_x_continuous(limits = c(xmin, xmax), breaks = breaks, labels=labels) +
  scale_y_continuous(expand = c(0,0), limits = c(1, ylim), breaks = c(1:ylim)) +
  geom_hline(yintercept=threshold_line, linetype="dashed", color = "red") +
  scale_color_manual(values = colors) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, y = "-log10(p)") + 
  ggtitle(title) +
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )

manhplot_zoom


#### LD plot in zoomed in range ####


tgr <- "f_Comamonadaceae Unknown Genus"
asv <- "asv_000260"
nit <- "LN"


## get hapmap v4 data for chromosome 10 range 23800k-24100k
haplotypes <- read.table("largedata/hapmap_chr10_23800k-24100k.txt",head=T,com="",sep="\t",row=1,na.string="NN")

haplotypes <- haplotypes %>% 
  rownames_to_column(var = "snp_id")



### subset samples with gwas log10p > 1

load("largedata/GWAS/Manhattan_plot_objects/f_Comamonadaceae_Unknown_Genus_LN.rda")

plot_obj$ps <- as.numeric(plot_obj$ps)

zoom <- plot_obj %>%
  filter(chr == "10" && (ps >= 23800000 && ps <= 24100000))


# renomly sample 200 SNPs in the window
set.seed(2021)
LD_sample <- haplotypes %>%
  filter(pos %in% zoom$ps) %>%
  sample_n(200)

## draw LD plot

fname <- "f_Comamonadaceae_Unknown_Genus_LN_chr10_23800k-24100k_200"
fpath <- paste0("figures/LD_plots/", fname, ".png")
rgb.palette <- colorRampPalette(rev(c("blue","orange" ,"red")), space = "rgb")
gene <- data.frame(t(LD_sample[,5:ncol(LD_sample)]))
gty <- makeGenotypes(gene,sep="")
png(file=fpath,res=300,width=1000,height=1000)
myld <- LDheatmap(gty,genetic.distances=LD_sample$pos,flip=TRUE, text=FALSE, color=rgb.palette(20), title=fname)
dev.off()






#### plot correlation of abundance vs canopy coverage, genotypes marked by min vs maj allele ####


load("data/abundance_vs_phenotype.rda")

load("cache/hmp_data_candidate_Comamonadaceae.rda") # haplotype information at peak SNP


cc_df_com <- abundance_vs_phenotype %>%
  filter(tax_group == "f_Comamonadaceae Unknown Genus") %>%
  filter(phenotype %in% c("CC_Aug12"))


tgr <- "f_Comamonadaceae Unknown Genus"
snp <- "S10_23954233"
chr <- "10"


allele_df_com <- hmp_data_candidates %>%
  filter(snp_id == "S10_23954233" & tax_group == "f_Comamonadaceae Unknown Genus") %>%
  ungroup() %>%
  dplyr::select(GX_name, allele) %>%
  unique()


cc_allele_df <- left_join(cc_df_com, allele_df_com)

colors <- c("major allele"="#ecb602", "minor allele"="#a45ee5", "NA"="#666666")


plot_dat_com <- cc_allele_df %>%
  mutate(allele = ifelse(is.na(allele), "NA", allele)) %>%
  filter(nitrogen == "LN")


## abundance vs canopy coverage
ggplot(plot_dat_com, aes( x = value,  y = blup_logrel)) +
  geom_point(aes(color = allele)) +
  geom_smooth(method = "lm", color="black") +
  ylab("Abundance") +
  xlab("Canopy Coverage") +
  scale_color_manual(values = colors) +
  scale_x_continuous(position = "top") +
  #facet_wrap(~nitrogen, scales = "free", nrow = 1) +
  stat_fit_glance(method = 'lm',
                  geom = 'text',
                  aes(label = paste0('p = ', round(..p.value.., 7))),
                  label.x = 15, label.y = "top") +
  theme_bw()



## abundance major vs minor allele

ggplot(filter(plot_dat_com, allele != "NA"), aes(x=nitrogen, y=blup_logrel, fill=allele)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(group = allele), label = "p.signif") +
  ylab("log(microbe relative abundance)") +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


## cc major vs minor allele

plot_dat_com$allele <- factor(plot_dat_com$allele, levels = c("minor allele", "major allele"))

ggplot(filter(plot_dat_com, allele != "NA"), aes(x=nitrogen, y=value, fill=allele)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(group = allele), label = "p.signif") +
  ylab("canopy coverage") +
  scale_fill_manual(values = colors) +
  coord_flip() +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())




### Gen gene expression data for chr 10 locus genes
#4 genotypes with N treatment


locus_genes <- c("Zm00001d023838", "Zm00001d023839", "Zm00001d023840")



fpkm <- read_csv("data/FPKM_Old_new_four_lines2.csv")
colnames(fpkm)[1] <- "gene"

N_plot_data <- fpkm %>%
  filter(gene %in% locus_genes) %>%
  ##dplyr::select(gene, contains("Root")) %>%
  pivot_longer(-gene, names_to = "Sample", values_to = "fpkm") %>%
  separate(Sample, c("nitrogen", "genotype", "tissue", "rep"), sep = "_") %>%
  mutate(nitrogen = ifelse(nitrogen == "HN", "+N", "-N"))

colors <- c("leaf"="#76ba1b", "Root" = "#dca85c")

ggplot(N_plot_data, aes(x=nitrogen, y=fpkm, fill=tissue)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge(jitter.width = 0.1)) +
  stat_compare_means(aes(group = tissue), label = "p.signif") +
  ylab("gene expression (FPKM)") +
  facet_wrap(~gene) +
  scale_fill_manual(values = colors) +
  theme_bw()



