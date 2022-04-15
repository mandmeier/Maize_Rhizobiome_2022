# heritability analysis

library("tidyverse")
library("phyloseq")
library("lme4")
library("ggpubr")
library("ggrepel")


#### calculate heritability ####

load("data/ps_grp.rda")

# calculate mean microbe abundance for each maize genotype

### filter asv table for 2019 data without checks
ps_19 <- subset_samples(ps_grp, year == "Y2019" & genotype != "CHECK")
asv_table <- rownames_to_column(data.frame(otu_table(ps_19)), var = "Sample_ID")
counts <- left_join(data.frame(sample_data(ps_19)), asv_table)


### calculate means of ASV counts across all subsamples by group
mean_counts <- counts %>%
  group_by(genotype, nitrogen, block, sp, spb) %>%
  summarize_each(funs(mean), -Sample_ID, -raw_seq_count, -filename, -year, -subsample, -collected_by_person, -pedigree, -subpopulation)


### calculate log relative abundance
asvtab <- mean_counts[, 7:ncol(mean_counts)]
logrel <- t(apply(asvtab, 1, function(x) log(x/sum(x) + 0.001)))
mean_counts_logrel <- cbind(data.frame(mean_counts[, 1:6]), logrel)



### split stdN, lowN samples, retain only genotypes with at least 2 reps (both blocks)
# (note that genotypes 38-11, A214N, and CI90C were planted twice and have 4 reps)

h2dat_stdN <- mean_counts_logrel %>%
  filter(nitrogen == "+N") %>%
  group_by(genotype) %>%
  add_tally(name="count") %>%
  filter(count >= 2) %>%
  dplyr::select(-count)


#save(h2dat_stdN, file = "cache/h2dat_stdN.rda")

# unique(as.character(h2dat_stdN$genotype))
# 206 genotypes

h2dat_lowN <- mean_counts_logrel %>%
  filter(nitrogen == "-N") %>%
  group_by(genotype) %>%
  add_tally(name="count") %>%
  filter(count >= 2) %>%
  dplyr::select(-count)

#save(h2dat_lowN, file = "cache/h2dat_lowN.rda")

# unique(as.character(h2dat_lowN$genotype))
# 206 genotypes


# calculate heritability using formula H2 = Vg/(Vg + Ve/6)
# Vg : variance of the genotype
# Ve : error variance
# divide by 6 for 6 reps (2 blocks x 3 subsamples)



## calculate H2 for +N
df <- h2dat_stdN
traits <- colnames(df)[7:ncol(df)]

H2_stdN <- data.frame()

for(trait in traits){
  #trait <- "asv_000013"
  f <- formula(paste0(trait, ' ~ (1|genotype) + (1|block)'))
  fit <- lmer(f, data=df)
  v <- as.data.frame(VarCorr(fit))
  Vg <- v$vcov[v$grp == "genotype"]
  Ve <- v$vcov[v$grp == "Residual"]
  H2 <- round(Vg/(Vg + Ve/6), 6)
  H2_stdN <- rbind(H2_stdN, data.frame("ASV"=trait, "H2_stdN_19"=H2))
}

H2_stdN




## calculate H2 for -N
df <- h2dat_lowN
traits <- colnames(df)[7:ncol(df)]

H2_lowN <- data.frame()

for(trait in traits){
  #trait <- "asv_000013"
  f <- formula(paste0(trait, ' ~ (1|genotype) + (1|block)'))
  fit <- lmer(f, data=df)
  v <- as.data.frame(VarCorr(fit))
  Vg <- v$vcov[v$grp == "genotype"]
  Ve <- v$vcov[v$grp == "Residual"]
  H2 <- round(Vg/(Vg + Ve/6), 6)
  H2_lowN <- rbind(H2_lowN, data.frame("ASV"=trait, "H2_lowN_19"=H2))
}

H2_lowN



## add heritability to group data
load("data/group_data.rda")
#group_data <- group_data %>%
# cleft_join(H2_stdN) %>%
#  left_join(H2_lowN)
#save(group_data, file = "data/group_data.rda")



# Compare heritability between HN and LN, paired T test

load("data/group_data.rda")
h2dat <- group_data %>%
  select(ASV, H2_lowN_19, H2_stdN_19) %>%
  filter(H2_lowN_19 > 0 & H2_stdN_19 > 0)

t.test(h2dat$H2_lowN_19, h2dat$H2_stdN_19, paired = TRUE, alternative = "greater" )
# p-value = 0.02147





#### plot correlation of microbe abundance vs heritability ####


HN <- group_data %>%
  dplyr::select(tax_group, mean_blup_stdN, H2_stdN_19) %>%
  #pivot_longer(starts_with("mean"), names_to = "mean_blup", values_to = "logrel") %>%
  mutate(nitrogen = "HN") %>%
  rename("abundance" = mean_blup_stdN, "heritability" = H2_stdN_19)

LN <- group_data %>%
  dplyr::select(tax_group, mean_blup_lowN, H2_lowN_19) %>%
  #pivot_longer(starts_with("mean"), names_to = "mean_blup", values_to = "logrel") %>%
  mutate(nitrogen = "LN") %>%
  rename("abundance" = mean_blup_lowN, "heritability" = H2_lowN_19)

plot_data <- rbind(HN, LN)



ggscatter(plot_data, x = "heritability", y = "abundance",
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.method = "pearson", #color = "nitrogen",
          xlab = "Heritability", ylab = "Mean Log(relative abundance)", alpha=0.5) +
  #scale_color_manual(values = c("#000080", "#C00001")) +
  #geom_hline(yintercept = 2, color = "black", linetype = "dashed") +
  facet_wrap(~nitrogen, nrow = 1) +
  theme_bw()


#### plot heritability +N vs -N ####


# permutation_test
# do a permutation test  to determine a threshold to distinguish heritable vs non heritable traits
# the threshold is the heritability score you would expect by chance

# 1) for each trait and N treatment: calculate "observed" H2
# 2) shuffle each genotype's ASV counts of this trait 1000 times
# 3) calculate 1000 permutation H2 scores
# 4) how many of the 1000 are larger than the observed H2?
# 5) p-value = #larger/1000
# 6) this yields 300 p values for each trait and treatment
# 7) mark quadrants: check if p < 0.05 for +N and -N

load("cache/h2dat_stdN.rda")
load("cache/h2dat_lowN.rda")



# Function to calculate H2

getH2 <- function(h2dat) {
  #h2dat <- h2dat_stdN
  traits <- colnames(h2dat)[7:ncol(h2dat)]
  #print(traits)
  H2_df <- data.frame()
  
  for(trait in traits){
    #trait <- "asv_000013"
    f <- formula(paste0(trait, ' ~ (1|genotype) + (1|block)'))
    fit <- suppressMessages(lmer(f, data=h2dat))
    v <- as.data.frame(VarCorr(fit))
    Vg <- v$vcov[v$grp == "genotype"]
    Ve <- v$vcov[v$grp == "Residual"]
    H2 <- round(Vg/(Vg + Ve/6), 6)
    H2_df <- rbind(H2_df, data.frame("ASV"=trait, "H2"=H2))
  }
  
  return(H2_df)
  
}

# function to shuffle abundances
shuffle <- function(h2dat) {
  h2dat$genotype <- sample(h2dat$genotype)
  return(h2dat)
}



#### calculate permutation p_values stdN

set.seed(2021)
H2_permutations_stdN <-  getH2(h2dat_stdN)
colnames(H2_permutations_stdN)[2] <- "obs"

for( i in c(1:1000)){
  print(paste("permutation", i))
  perm <-  getH2(shuffle(h2dat_stdN))
  colnames(perm)[2] <- paste0("p",i)
  H2_permutations_stdN <- suppressMessages(left_join(H2_permutations_stdN, perm))
}

p_values_stdN <- H2_permutations_stdN %>%
  rowwise() %>%
  mutate(perm_p_stdN = (sum(c_across(3:length(H2_permutations_stdN))> obs)+1)/(length(H2_permutations_stdN)-1)) %>%
  dplyr::select(ASV, perm_p_stdN)



#### calculate permutation p_values lowN

set.seed(2021)
H2_permutations_lowN <-  getH2(h2dat_lowN)
colnames(H2_permutations_lowN)[2] <- "obs"

for( i in c(1:1000)){
  print(paste("permutation", i))
  perm <-  getH2(shuffle(h2dat_lowN))
  colnames(perm)[2] <- paste0("p",i)
  H2_permutations_lowN <- suppressMessages(left_join(H2_permutations_lowN, perm))
}

p_values_lowN <- H2_permutations_lowN %>%
  rowwise() %>%
  mutate(perm_p_lowN = (sum(c_across(3:length(H2_permutations_lowN))> obs)+1)/(length(H2_permutations_lowN)-1)) %>%
  dplyr::select(ASV, perm_p_lowN)


### save permutation data

#save(H2_permutations_stdN, file="cache/H2_permutations_stdN.rda")
#save(H2_permutations_lowN, file="cache/H2_permutations_lowN.rda")


### add pvalues to group data

#load("data/group_data.rda")

#group_data <- group_data %>%
#  left_join(p_values_stdN) %>%
#  left_join(p_values_lowN)


#load("data/group_data.rda")


group_data <- group_data %>%
  mutate(diffab_group = ifelse(log2FoldChange < 0, "negative", "positive")) %>%
  mutate(diffab_group = ifelse(padj >= 0.05, "n.s.", diffab_group))



group_data$sign_lowN <- ifelse(group_data$perm_p_lowN < 0.05, "sign", "n.s.")
group_data$sign_stdN <- ifelse(group_data$perm_p_stdN < 0.05, "sign", "n.s.")

group_data$sign_group <- "none"
group_data$sign_group <- ifelse(group_data$sign_lowN == "sign", "lowN", group_data$sign_group)
group_data$sign_group <- ifelse(group_data$sign_stdN == "sign", "stdN", group_data$sign_group)
group_data$sign_group <- ifelse(group_data$sign_stdN == "sign" & group_data$sign_lowN == "sign", "both", group_data$sign_group)
group_data$sign_group <- factor(group_data$sign_group, levels = c("stdN", "lowN", "both", "none"))

group_data <- group_data %>%
  dplyr::select(-sign_lowN, -sign_stdN)

#save(group_data, file = "data/group_data.rda")



ggplot(group_data, aes(x=H2_lowN_19, y=H2_stdN_19, color = sign_group)) +
  geom_density_2d(color="#cccccc") +
  geom_point(alpha=0.7) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_smooth(method='lm', se=TRUE, color="#008000") +
  scale_color_manual(values= c(lowN="#C00001", stdN="#000080", none="#999999", both="#800080"), labels=c(lowN="heritable under -N", stdN="heritable under +N", both="heritable under both treatments", none="heritability not significant")) +
  theme_classic() +
  ylab("Heritability under +N") +
  xlab("Heritability under -N")





#### plot heritability +N vs -N, label most heritable microbial groups ####


top <- group_data %>%
  filter(H2_lowN_19 > 0.6 & H2_stdN_19 > 0.6)


ggplot(group_data, aes(x=H2_lowN_19, y=H2_stdN_19, color = sign_group)) +
  geom_density_2d(color="#cccccc") +
  geom_point(alpha=0.7) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_vline(xintercept = 0.6) +
  geom_hline(yintercept = 0.6) +
  geom_text_repel(data=top, aes(label=tax_group)) +
  geom_smooth(method='lm', se=TRUE, color="#008000") +
  scale_color_manual(values= c(lowN="#C00001", stdN="#000080", none="#999999", both="#800080"), labels=c(lowN="heritable under -N", stdN="heritable under +N", both="heritable under both treatments", none="heritability not significant")) +
  theme_classic() +
  ylab("Heritability under +N") +
  xlab("Heritability under -N")


## table for supplementary figure
View(top)









