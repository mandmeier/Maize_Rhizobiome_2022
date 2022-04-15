# BLUP calculation

library("tidyverse")
library("phyloseq")
library("lme4")

devtools::install_github("jyanglab/g3tools")
library("g3tools")


#### BLUP calculation, ASV level ####

load("data/ps_asv.rda")


### use 2019 data only, exclude checks
ps19 <- subset_samples(ps_asv, year == "Y2019" & genotype != "CHECK")
asv_table <- rownames_to_column(data.frame(otu_table(ps19)), var = "Sample_ID")
counts <- left_join(data.frame(sample_data(ps19)), asv_table)

### calculate means of ASV counts of all subsamples by group
mean_counts <- counts %>%
  group_by(genotype, nitrogen, block, sp, spb) %>%
  summarize_each(funs(mean), -year, -Sample_ID, -raw_seq_count, -filename, -subsample, -collected_by_person, -pedigree, -subpopulation)


### calculate log relative abundance
asvtab <- mean_counts[, 7:ncol(mean_counts)]
logrel <- t(apply(asvtab, 1, function(x) log(x/sum(x) + 0.001)))
mean_counts_logrel <- cbind(data.frame(mean_counts[, 1:6]), logrel)



### BLUP +N

if(!(exists("blup_stdN"))){
  blup_stdN <- data.frame(genotype = sort(unique(df$genotype)))
}

for (asvid in colnames(asvtab)){
  
  df <- mean_counts_logrel %>%
    filter(nitrogen == "+N") %>%
    select_("genotype", "block", "sp", "spb", "row", asvid)
  
  colnames(df)[ncol(df)] <- "asv_abundance"
  
  get_BLUP(data = df, model = asv_abundance ~ (1 | genotype) + (1 | block) + (1 | sp) + (1 | spb), which.factor = "genotype",
           outfile = "cache/BLUP/tmp_blup.csv")
  
  blup <- read.csv("cache/BLUP/tmp_blup.csv")
  colnames(blup)[1] <- "genotype"
  colnames(blup)[2] <- asvid
  
  blup_stdN <- left_join(blup_stdN, blup)
  
}



### save data
write_csv(blup_stdN, file="data/blup_stdN_3618_asvs.csv")




### BLUP -N

if(!(exists("blup_lowN"))){
  blup_lowN <- data.frame(genotype = sort(unique(df$genotype)))
}

for (asvid in colnames(asvtab)){
  
  df <- mean_counts_logrel %>%
    filter(nitrogen == "-N") %>%
    select_("genotype", "block", "sp", "spb", "row", asvid)
  
  colnames(df)[ncol(df)] <- "asv_abundance"
  
  get_BLUP(data = df, model = asv_abundance ~ (1 | genotype) + (1 | block) + (1 | sp) + (1 | spb), which.factor = "genotype",
           outfile = "cache/BLUP/tmp_blup.csv")
  
  blup <- read.csv("cache/BLUP/tmp_blup.csv")
  colnames(blup)[1] <- "genotype"
  colnames(blup)[2] <- asvid
  
  blup_lowN <- left_join(blup_lowN, blup)
  
}


### save data
write_csv(blup_lowN, file="data/blup_lowN_3618_asvs.csv")




#### BLUP calculation, microbial group level ####


load("data/ps_grp.rda")


### use 2019 data only, exclude checks
ps19 <- subset_samples(ps_grp, year == "Y2019" & genotype != "CHECK")
asv_table <- rownames_to_column(data.frame(otu_table(ps19)), var = "Sample_ID")
counts <- left_join(data.frame(sample_data(ps19)), asv_table)

### calculate means of ASV counts of all subsamples by group
mean_counts <- counts %>%
  group_by(genotype, nitrogen, block, sp, spb) %>%
  summarize_each(funs(mean), -year, -Sample_ID, -raw_seq_count, -filename, -subsample, -collected_by_person, -pedigree, -subpopulation)


### calculate log relative abundance
asvtab <- mean_counts[, 7:ncol(mean_counts)]
logrel <- t(apply(asvtab, 1, function(x) log(x/sum(x) + 0.001)))
mean_counts_logrel <- cbind(data.frame(mean_counts[, 1:6]), logrel)



### BLUP std N

if(!(exists("blup_stdN"))){
  blup_stdN <- data.frame(genotype = sort(unique(df$genotype)))
}

for (asvid in colnames(asvtab)){
  
  df <- mean_counts_logrel %>%
    filter(nitrogen == "+N") %>%
    select_("genotype", "block", "sp", "spb", "row", asvid)
  
  colnames(df)[ncol(df)] <- "asv_abundance"
  
  get_BLUP(data = df, model = asv_abundance ~ (1 | genotype) + (1 | block) + (1 | sp) + (1 | spb), which.factor = "genotype",
           outfile = "cache/BLUP/tmp_blup.csv")
  
  blup <- read.csv("cache/BLUP/tmp_blup.csv")
  colnames(blup)[1] <- "genotype"
  colnames(blup)[2] <- asvid
  
  blup_stdN <- left_join(blup_stdN, blup)
  
}



### BLUP low N

if(!(exists("blup_lowN"))){
  blup_lowN <- data.frame(genotype = sort(unique(df$genotype)))
}

for (asvid in colnames(asvtab)){
  
  df <- mean_counts_logrel %>%
    filter(nitrogen == "-N") %>%
    select_("genotype", "block", "sp", "spb", "row", asvid)
  
  colnames(df)[ncol(df)] <- "asv_abundance"
  
  get_BLUP(data = df, model = asv_abundance ~ (1 | genotype) + (1 | block) + (1 | sp) + (1 | spb), which.factor = "genotype",
           outfile = "cache/BLUP/tmp_blup.csv")
  
  blup <- read.csv("cache/BLUP/tmp_blup.csv")
  colnames(blup)[1] <- "genotype"
  colnames(blup)[2] <- asvid
  
  blup_lowN <- left_join(blup_lowN, blup)
  
}






### use T1-T150 for colnames (to avoid confusing tax groups with ASVs)
rep_asvs <- data.frame("ASV" = colnames(blup_stdN)[-1])
traits <- left_join(rep_asvs, select(group_data, ASV, trait))
newnames <- c("genotype", traits$trait)
colnames(blup_stdN) <- newnames
colnames(blup_lowN) <- newnames

### save data

write_csv(blup_stdN, file="data/blup_stdN_150_tax_groups.csv")
write_csv(blup_lowN, file="data/blup_lowN_150_tax_groups.csv")



### calculae mean BLUP

mean_blup_stdN <- blup_stdN %>%
  select(-1) %>%
  colMeans(na.rm = TRUE)
mean_blup_stdN <- data.frame("ASV" = names(mean_blup_stdN), "mean_blup_stdN" = mean_blup_stdN)

mean_blup_lowN <- blup_lowN %>%
  select(-1) %>%
  colMeans(na.rm = TRUE)
mean_blup_lowN <- data.frame("ASV" = names(mean_blup_lowN), "mean_blup_lowN" = mean_blup_lowN)


mean_blup <- left_join(mean_blup_stdN, mean_blup_lowN)



### add mean blup to group data
load("data/group_data.rda")
group_data <- left_join(group_data, mean_blup)
#save(group_data, file = "data/group_data.rda")





