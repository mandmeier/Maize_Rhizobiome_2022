# GWAS analysis

library("tidyverse")
library("scales")


#### plot genome overview: loci with significant associations with one or more microbial groups ####

load("data/gwas_dat_threshold7.rda")
load("data/group_data.rda")


## lengths of chromosomes
chr_len_10k <- c(308460000, 243680000, 238020000, 250340000, 226360000, 181360000, 185810000, 182420000, 163010000, 152440000)

nCHR <- length(unique(gwas_dat$chr))
chindex <- c(0) ## this number needs to be added to position within each chromosome
for (i in c(1:nCHR)){
  #print(i)
  chindex[i+1] <- sum(chr_len_10k[1:i])
}
chindex <- chindex[1:10]
bpadd <- data.frame("chr" = c(1:nCHR), "chindex" = chindex)


plot_dat <- gwas_dat %>%
  #dplyr::select(-X1) %>%
  #mutate(log10p = -log10(p_wald)) %>%
  #filter(log10p >= 5) %>%
  left_join(bpadd, by = "chr") %>%
  mutate(psabs = chindex + ps) %>%
  mutate(psabs = as.numeric(gsub('.{4}$', '', psabs))) %>% ## pos. in 10kb units
  mutate(ps = as.character(ps)) %>%
  mutate(bin = as.numeric(gsub('.{4}$', '', ps)) + 1) %>% # absolute bin position in 10kb units
  group_by(chr, bin, nitrogen) %>%
  add_count(count = sum(-log10(p_wald) >= 7.187346), name="total_sign_snps") %>%
  #summarize(total_p = sum(log10p))
  filter(total_sign_snps >= 1) %>%
  mutate(mean_p = mean(log10p)) %>%
  ##mutate(total_mean_p = mean(log10p)) %>%
  group_by(chr, bin, nitrogen, trait) %>%
  #add_count(count = sum(-log10(p_wald) >= 5), name="sign_snps") %>%
  ##mutate(mean_p = mean(log10p)) %>%
  rename(ASV=trait) %>%
  left_join(dplyr::select(group_data, ASV, tax_group)) %>%
  #mutate(color = ifelse(chr %in% c(1,3,5,7,9), "#276FBF", "#183059")) %>%
  mutate(chr = as.character(chr)) 


# get taxa per bin
taxa_per_bin <- plot_dat %>%
  dplyr::select(nitrogen, chr, bin, tax_group) %>%
  unique() %>%
  group_by(nitrogen, chr, bin) %>%
  add_tally(name="taxa_per_bin") 

plot_dat <- plot_dat %>%
  left_join(taxa_per_bin)

plot_dat$nitrogen <- factor(plot_dat$nitrogen, levels = c("LN", "HN"))


## count microbe associated plant loci MAPLs

## count all loci with GWAS signals

all_gwas_signals <- plot_dat

nrow(all_signals) # 1089 significant signals

microbial_groups_1089 <- unique(all_gwas_signals$tax_group) # 104 microbial groups

mapls_1089 <- unique(all_gwas_signals$bin) # 613 bins

save(all_gwas_signals, file = "cache/all_gwas_signals.rda")

# count loci with at least 2 significant SNPs per 10kb bin

mapl_gwas_signals <- plot_dat %>%
  filter(total_sign_snps >= 2) 

nrow(mapl_gwas_signals) # 586 significant signals

microbial_groups_586 <- unique(mapl_gwas_signals$tax_group) # 35 microbial groups

mapls_586 <- unique(mapl_gwas_signals$bin) # 119 bins

save(mapl_gwas_signals, file = "cache/mapl_gwas_signals.rda")


colors = c("1" = "#276FBF", "3" = "#276FBF", "5" = "#276FBF", "7" = "#276FBF", "9" = "#276FBF",
           "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059")

## find position for labels
axis_set <- plot_dat %>% 
  group_by(chr) %>% 
  summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)


overview_plot <-  ggplot(mapl_gwas_signals, aes(x = psabs, y = mean_p, color = chr)) +
  geom_bar(stat="identity", position="dodge") +
  geom_point(aes(size = taxa_per_bin), alpha = 0.3) +
  geom_hline(yintercept = 7.187346, color = "red", linetype = "dashed")
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  scale_y_continuous(limits = c(7.187346, 11), oob = rescale_none) +
  scale_color_manual(values = colors) +
  scale_size_continuous(breaks = c(1, 2, 3)) +
  facet_wrap(~nitrogen, nrow = 2) +
  labs(x = NULL, y = "mean -Log10(p)") + 
  geom_vline(xintercept =  c(0, 212850), color = "red", size=0.1) + # to combine plots
  geom_hline(yintercept = 7, size=0.1) + 
  theme_bw() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

overview_plot



#### GWAS for high-level traits ####


chr_len_10k <- c(308460000, 243680000, 238020000, 250340000, 226360000, 181360000, 185810000, 182420000, 163010000, 152440000)


nCHR <- 10
chindex <- c(0) ## this number needs to be added to position within each chromosome
for (i in c(1:nCHR)){
  #print(i)
  chindex[i+1] <- sum(chr_len_10k[1:i])
}
chindex <- chindex[1:10]
bpadd <- data.frame("chr" = c(1:nCHR), "chindex" = chindex)



## function to draw manhattan plot for any trait
draw_manhattan <- function(gwas_dat, trait, threshold_line = 7.187346) {
  title <- trait
  # add bins
  plot_obj <- gwas_dat %>%
    mutate(log10p = -log10(p_wald)) %>%
    left_join(bpadd, by = "chr") %>%
    mutate(psabs = chindex + ps) %>%
    mutate(ps = as.character(ps)) %>%
    mutate(bin = as.numeric(gsub('.{4}$', '', ps)) + 1) %>% # absolute bin position in 10kb units
    group_by(chr, bin) %>%
    add_count(count = sum(-log10(p_wald) >= threshold_line), name="sign_snps") %>%
    mutate(chr = as.character(chr))

  ## find position for labels
  axis_set <- plot_obj %>% 
    group_by(chr) %>% 
    summarize(center = (max(psabs) + min(psabs)) / 2) # find middle of chromosome (to place labels)
  
  ## find range for y axis
  ylim <- floor(max(plot_obj$log10p, na.rm = TRUE)) + 1
  
  if (ylim < threshold_line){
    ylim <- 8
  }
  
  colors = c("1" = "#58CCED", "3" = "#58CCED", "5" = "#58CCED", "7" = "#58CCED", "9" = "#58CCED",
             "2" = "#183059", "4" = "#183059", "6" = "#183059", "8" = "#183059", "10" = "#183059")
  
  manhplot <- ggplot(plot_obj, aes(x = psabs, y = log10p, color=chr)) +
    geom_point(alpha = 0.75) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
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
      panel.grid.minor.x = element_blank()#,
    )
  
  return(manhplot)
  
}


## example, draw one plot
#gwas_dat <- read_delim("largedata/GWAS_high_level_traits/LN_PC1.assoc.txt", delim = "\t")
#trait <- "LN_PC1"
#draw_manhattan(gwas_dat, trait, 7.187346)


## draw all manhattan plots for 2x 14 traits (this will take several hours)

traits <- c("HN_Shannon", "LN_Shannon",
            "HN_Observed", "LN_Observed",
            "HN_InvSimpson", "LN_InvSimpson",
            "HN_Fisher", "LN_Fisher",
            "HN_PC1", "LN_PC1",
            "HN_PC2", "LN_PC2",
            "HN_PC3", "LN_PC3",
            "HN_PC4", "LN_PC4",
            "HN_PC5", "LN_PC5",
            "HN_PC6", "LN_PC6",
            "HN_PC7", "LN_PC7",
            "HN_PC8", "LN_PC8",
            "HN_PC9", "LN_PC9",
            "HN_PC10", "LN_PC10")


for(trait in traits){
  
  infile <- paste0("largedata/GWAS_high_level_traits/", trait,".assoc.txt")
  print(infile)
  outfile <- paste0("figures/Manhattan_lots/", trait,".png")
  print(outfile)
  
  gwas_dat <- read_delim(infile, delim = "\t")
  
  p <- draw_manhattan(gwas_dat, trait, 7.187346)
  
  print("plotting")
  
  png(outfile)
  print(p)
  dev.off()
  
}


