#abundance_vs_phenotype



library("tidyverse")
library("ggpubr")


#load("data/group_data.rda")

#test <- read_csv("data/blup_stdN_150_tax_groups.csv")


#data_blup_LN <- read_delim(file="data/Blup_and_Heritability/data_blup_LN_result.txt", delim = "\t")



#### for each microbial group and plant trait calculate correlation between abundance and phenotype ###


load("data/abundance_vs_phenotype.rda")


corr_data <- abundance_vs_phenotype %>%
  dplyr::select(nitrogen, tax_group, phenotype) %>%
  unique()

corr_data$pearson <- 0
corr_data$p_value <- 0


get_correlation <- function(n, m, p, var = "blup_logrel"){
  #n <- "LN"
  #m <- "Burkholderia sp 2"
  #p <- "cob_length"
  #var <- "blup_logrel"
  #p <- corr_data[i,]$phenotype
  nmp <- data.frame(filter(abundance_vs_phenotype, nitrogen == n & tax_group == m & phenotype == p))
  cor <- cor.test(nmp[, var], nmp[, "value"], method=c("pearson"))
  return(cor)
}


for (i in c(1:nrow(corr_data))){
  print(i)
  cor <- get_correlation(corr_data[i,]$nitrogen, corr_data[i,]$tax_group, corr_data[i,]$phenotype)
  pears <- cor$estimate[[1]]
  p_value <- cor$p.value
  corr_data[i,]$pearson <- pears
  corr_data[i,]$p_value <- p_value
}


#save(corr_data, file="data/corr_data.rda")



#### plot correlation of microbe abundance with phenotypes for 2 N treatments and 150 microbial groups ###


load("data/corr_data.rda")

# exclude cob data (separately published)
corr_data_17 <- corr_data %>%
  filter(!(phenotype %in% c("cob_length", "cob_weight", "cob_width")))


corr_data_17$phenotype <- factor(corr_data_17$phenotype,
                                 levels = c("CC_Aug12", "ExG_Aug12",
                                            "CHL", "DW", "FW", "LA",
                                            "B", "Ca", "Cu", "Fe", "K", "Mg", "Mn", "N", "P", "S", "Zn"))

ggplot(corr_data_17, aes(x=pearson, y=-log10(p_value), color=nitrogen)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 2, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("#000080", "#C00001")) +
  facet_wrap(~phenotype+nitrogen, ncol = 10) +
  xlab("pearson correlation [r]") +
  theme_bw()




#### plot correlation of microbe abundance with canopy coverage for 2 N treatments and 150 microbial groups ###



CC <- corr_data %>%
  filter(phenotype == "CC_Aug12") %>%
  mutate(class = ifelse(p_value < 0.01, ifelse(pearson < 0, "negative", "positive"), "ns"))


ggplot(CC, aes(x=pearson, y=-log10(p_value), color=class)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 2, color = "black", linetype = "dashed") +
  scale_color_manual(values = c("ns"= "#333333", "negative" = "#FF0000", "positive" = "#00B050")) +
  facet_wrap(~nitrogen, nrow = 1) +
  xlab("pearson correlation") +
  theme_bw()






#### plot heritability vs p_value of correlation with canopy coverage ####


load("data/group_data.rda")


H2_vs_CC <- corr_data %>%
  filter(phenotype == "CC_Aug12") %>%
  mutate(class = ifelse(p_value < 0.01, ifelse(pearson < 0, "negative", "positive"), "ns")) %>%
  left_join(group_data) %>%
  #dplyr::select(nitrogen, tax_group, pearson, p_value, H2_stdN_19, H2_lowN_19) %>%
  mutate(log10p=-log10(p_value)) %>%
  mutate(heritability = ifelse(nitrogen == "HN", H2_stdN_19, H2_lowN_19)) %>%
  mutate(abundance = ifelse(nitrogen == "HN", mean_blup_stdN, mean_blup_lowN)) 



#### plot correlation
ggplot(H2_vs_CC, aes( x = heritability, y = log10p, color = class)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 2, color = "black", linetype = "dashed") +
  geom_smooth(aes(x = heritability, y = log10p), method='lm', formula= y~x, inherit.aes = FALSE, color = "#333333") +
  scale_color_manual(values = c("ns"= "#333333", "negative" = "#FF0000", "positive" = "#00B050")) +
  facet_wrap(~nitrogen, nrow = 1) +
  ylab("-log10(p_value)") +
  theme_bw()








