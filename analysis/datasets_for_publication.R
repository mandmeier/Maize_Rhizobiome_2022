## Datasets for publication

library("tidyverse")
library("phyloseq")
library("stringr")


## Dataset 1

load("cache/ps_core5t.rda")


sdat <- data.frame(sample_data(ps_core5t))

sdat <- sdat %>%
  dplyr::select(-collected_by_person, -raw_seq_count) %>%
  mutate(nitrogen = ifelse(nitrogen == "+N", "HN", "LN")) %>%
  rename(nitrogen_treatment = nitrogen) %>%
  rename(SRA_filename = filename) %>%
  rename(split_plot = sp) %>%
  rename(split_plot_block = spb) %>%
  rename(subplot_id = row) %>%
  mutate(quadrant = str_sub(subplot_id, 1, 1)) %>%
  mutate(range = str_sub(subplot_id, 2, 2)) %>%
  mutate(row = str_sub(subplot_id, -2, -1)) %>%
  mutate(genotype = ifelse(genotype == "CHECK", "Check_B73xMo17", as.character(genotype))) %>%
  mutate(genotype = ifelse(genotype == "CHECK2", "Check_B37xMo17", as.character(genotype))) %>%
  rename(GRIN_Accession = pedigree) %>%
  mutate(GRIN_Accession = ifelse(GRIN_Accession == "B73xMo17", NA, as.character(GRIN_Accession))) %>%
  mutate(GRIN_Accession = ifelse(GRIN_Accession == "B37xMo17", NA, as.character(GRIN_Accession))) %>%
  dplyr::select(Sample_ID, SRA_filename, year, genotype, GRIN_Accession, subpopulation, nitrogen_treatment,
                block, split_plot, split_plot_block, subplot_id, quadrant, range, row, subsample)




asvtab <- data.frame(otu_table(ps_core5t))

asvtab <- asvtab %>%
  rownames_to_column(var = "Sample_ID")


dataset1 <- left_join(sdat, asvtab)

write_csv(dataset1, file = "data/Datasets/dataset1.csv", na = "")



## Dataset 2


taxtab <- data.frame(tax_table(ps_core5t))

taxtab <- taxtab %>%
  rownames_to_column(var = "ASV") %>%
  dplyr::select(ASV, Phylum, Class, Order, Family, Genus, Species, tax_group) %>%
  rename(Functional_Tax_Group = tax_group)


load("cache/sequences.rda")

head(names(sequences))


seq_3626 <- sequences[sequences %in% taxtab$ASV]


seqs <- data.frame("ASV" = seq_3626, "Sequence" = names(seq_3626))


dataset2 <- left_join(taxtab, seqs)


write_csv(dataset2, file = "data/Datasets/dataset2.csv", na = "")


## Dataset 3


load("data/group_data.rda")

grpdat <- group_data %>%
  dplyr::select(ASV, tax_group, unique_ASVs, count_stdN, count_lowN, log2FoldChange, lfcSE, padj, H2_stdN_19, H2_lowN_19)


load("cache/gctb_selection_values.rda")


svalues_HN <- selection_values %>%
  filter(grepl('asv_', trait)) %>%
  filter(treatment == "HN") %>%
  dplyr::select(trait, Mean, SD)

colnames(svalues_HN) <- c("ASV", "Selection_coeff_HN", "Selection_coeff_HN_sd")



svalues_LN <- selection_values %>%
  filter(grepl('asv_', trait)) %>%
  filter(treatment == "LN") %>%
  dplyr::select(trait, Mean, SD)

colnames(svalues_LN) <- c("ASV", "Selection_coeff_LN", "Selection_coeff_LN_sd")


grpdat <- grpdat %>%
  left_join(svalues_HN) %>%
  left_join(svalues_LN) %>%
  dplyr::select(-ASV)



colnames(grpdat) <- c("Microbial_Group", "unique_ASVs", "total_ASV_count_HN", "total_ASV_count_LN",
                      "deseq2_log2FoldChange_HNvsLN", "deseq2_lfc_SE", "deseq2_lfc_padj", "Heritability_HN", "Heritability_LN",
                      "Selection_coeff_HN", "Selection_coeff_HN_sd", "Selection_coeff_LN", "Selection_coeff_LN_sd")


grpdat_HN <- grpdat %>%
  mutate(nitrogen_treatment = "HN") %>%
  dplyr::select(nitrogen_treatment, Microbial_Group, unique_ASVs, total_ASV_count_HN, deseq2_log2FoldChange_HNvsLN, deseq2_lfc_SE, deseq2_lfc_padj, Heritability_HN, Selection_coeff_HN, Selection_coeff_HN_sd)
colnames(grpdat_HN) <- c("nitrogen_treatment", "Microbial_Group", "unique_ASVs_in_group", "total_ASV_count", "deseq2_log2FoldChange_HNvsLN", "deseq2_lfc_SE", "deseq2_lfc_padj", "Heritability", "Selection_coeff", "Selection_coeff_sd")


grpdat_LN <- grpdat %>%
  mutate(nitrogen_treatment = "LN") %>%
  dplyr::select(nitrogen_treatment, Microbial_Group, unique_ASVs, total_ASV_count_LN, deseq2_log2FoldChange_HNvsLN, deseq2_lfc_SE, deseq2_lfc_padj, Heritability_LN, Selection_coeff_LN, Selection_coeff_LN_sd)
colnames(grpdat_LN) <- c("nitrogen_treatment", "Microbial_Group", "unique_ASVs_in_group", "total_ASV_count", "deseq2_log2FoldChange_HNvsLN", "deseq2_lfc_SE", "deseq2_lfc_padj", "Heritability", "Selection_coeff", "Selection_coeff_sd")


dataset3 <- rbind(grpdat_HN, grpdat_LN)

dataset3 <- dataset3 %>%
  rename( total_read_count = total_ASV_count)

## canopy coverage correlation


load("data/corr_data.rda")

corrdat <- corr_data %>%
  filter(!(phenotype %in%c("cob_length", "cob_width", "cob_weight"))) %>%
  mutate(phenotype2 = phenotype)

corr <- corrdat %>%
  dplyr::select(nitrogen, tax_group, phenotype, pearson) %>%
  pivot_wider(names_from = phenotype, values_from = pearson)

colnames(corr) <- paste0("corr_", colnames(corr),"_r")
colnames(corr)[1] <- "nitrogen_treatment"
colnames(corr)[2] <- "Microbial_Group"


p <- corrdat %>%
  dplyr::select(nitrogen, tax_group, phenotype, p_value) %>%
  pivot_wider(names_from = phenotype, values_from = p_value)


colnames(p) <- paste0("corr_", colnames(p),"_p")
colnames(p)[1] <- "nitrogen_treatment"
colnames(p)[2] <- "Microbial_Group"

carrelations <- left_join(corr, p) %>%
  dplyr::select(nitrogen_treatment, Microbial_Group,
                corr_CC_Aug12_r, corr_CC_Aug12_p,
                corr_ExG_Aug12_r, corr_ExG_Aug12_p,
                corr_B_r, corr_B_p,
                corr_Ca_r, corr_Ca_p,
                corr_CHL_r, corr_CHL_p,
                corr_Cu_r, corr_Cu_p,
                corr_DW_r, corr_DW_p,
                corr_Fe_r, corr_Fe_p,
                corr_FW_r, corr_FW_p,
                corr_K_r, corr_K_p,
                corr_LA_r, corr_LA_p,
                corr_Mg_r, corr_Mg_p,
                corr_Mn_r, corr_Mn_p,
                corr_N_r, corr_N_p,
                corr_P_r, corr_P_p,
                corr_S_r, corr_S_p,
                corr_Zn_r, corr_Zn_p)

dataset3 <- left_join(dataset3, carrelations)



write_csv(dataset3, file = "data/Datasets/dataset3.csv", na = "")




## Dataset 4


blup_LN <- read_csv("data/blup_lowN_150_tax_groups.csv")

blup_HN <- read_csv("data/blup_lowN_150_tax_groups.csv")

names <- group_data %>%
  dplyr::select(trait, tax_group)

blup_HN <- blup_HN %>%
  pivot_longer(-genotype, names_to = "trait", values_to = "value") %>%
  left_join(names) %>%
  dplyr::select(genotype, tax_group, value) %>%
  rename(trait = tax_group) %>%
  mutate(unit = "BLUP of ln(relative microbe abundance)") %>%
  mutate(type = "rhizobiome trait") %>%
  mutate(nitrogen_treatment = "HN")



blup_LN <- read_csv("data/blup_lowN_150_tax_groups.csv")

blup_LN <- blup_LN %>%
  pivot_longer(-genotype, names_to = "trait", values_to = "value") %>%
  left_join(names) %>%
  dplyr::select(genotype, tax_group, value) %>%
  rename(trait = tax_group) %>%
  mutate(unit = "BLUP of ln(relative microbe abundance)") %>%
  mutate(type = "rhizobiome trait") %>%
  mutate(nitrogen_treatment = "LN")


rhizobiome_traits <- rbind(blup_HN, blup_LN)



## Semra yield component traits, Eric imaging traits
load("data/yield_analysis/microbe_counts_and_yield.rda")
load("data/group_data.rda")

yield_traits <- microbe_counts_and_yield %>%
  rename(MM_name=genotype) %>%
  mutate(nitrogen = ifelse(nitrogen == "+N", "HN", "LN")) %>%
  left_join(dplyr::select(group_data, ASV, tax_group))



unique(yield_traits$MM_name)


## yufeng nutrient traits
nutrient_traits <- read_delim(file="data/yield_analysis/phenotype.txt", delim="\t")
nutrient_traits_HN <- nutrient_traits %>%
  dplyr::select(ID, starts_with("HN")) %>%
  mutate(nitrogen="HN")
colnames(nutrient_traits_HN) <- c("GX_name", "B", "Ca", "CHL", "Cu", "DW", "Fe", "FW", "K", "LA", "Mg", "Mn", "N", "P", "S", "Zn", "nitrogen")

nutrient_traits_LN <- nutrient_traits %>%
  dplyr::select(ID, starts_with("LN")) %>%
  mutate(nitrogen="LN")

colnames(nutrient_traits_LN) <- c("GX_name", "B", "Ca", "CHL", "Cu", "DW", "Fe", "FW", "K", "LA", "Mg", "Mn", "N", "P", "S", "Zn", "nitrogen")

nutrient_traits <- rbind(nutrient_traits_LN, nutrient_traits_HN)


yield <- yield_traits %>%
  dplyr::select(MM_name, nitrogen, CC_Aug12, ExG_Aug12) %>%
  unique()


names <- read_csv("data/BG_MM_Gen_names.csv")

nutrient <- nutrient_traits %>%
  left_join(dplyr::select(names, MM_name, GX_name)) %>%
  unique()



info <- data.frame("trait" = c("CC_Aug12", "ExG_Aug12", "CHL", "DW", "FW", "LA", "B","Ca", "Cu", "Fe", "K", "Mg", "Mn", "N", "P", "S", "Zn"),
                   "unit" = c("canopy coverage (%)", "excess green index", "chlorophyll content (µmol/m2)", "leaf dry weight (g)", "leaf fresh weight (g)", "leaf area (m2)",
                              "micronutrient (%)", "micronutrient (%)", "micronutrient (%)", "micronutrient (%)", "micronutrient (%)", "micronutrient (%)", "micronutrient (%)", "micronutrient (%)", "micronutrient (%)", "micronutrient (%)", "micronutrient (%)"))


plant_traits <- yield %>%
  left_join(nutrient) %>%
  rename(genotype = MM_name) %>%
  rename(nitrogen_treatment = nitrogen) %>%
  dplyr::select(genotype, nitrogen_treatment, CC_Aug12, ExG_Aug12, CHL, DW, FW, LA, B, Ca, Cu, Fe, K, Mg, Mn, N, P, S, Zn) %>%
  pivot_longer(c(-genotype, -nitrogen_treatment), names_to = "trait", values_to = "value") %>%
  left_join(info) %>%
  mutate(type = "plant trait")






rhizobiome_traits <- rhizobiome_traits %>%
  dplyr::select(nitrogen_treatment, genotype, trait, value, unit, type)

plant_traits <- plant_traits %>%
  dplyr::select(nitrogen_treatment, genotype, trait, value, unit, type)


dataset4 <- rbind(plant_traits, rhizobiome_traits)



dataset4 <- dataset4 %>%
  filter(!(is.na(value)))



write_csv(dataset4, file = "data/Datasets/dataset4.csv", na = "")




## Dataset 5


load("cache/all_gwas_signals.rda")

load("cache/MAPL_genes_613.rda")


MAPL_traits <- all_gwas_signals %>%
  ungroup() %>%
  dplyr::select(nitrogen, chr, bin, total_sign_snps, tax_group, taxa_per_bin) %>%
  mutate(MAPL_id = paste0(nitrogen,"_chr",chr,"_",bin)) %>%
  mutate(start_pos = (bin-1)*10000+1, end_pos = (bin)*10000) %>%
  select(MAPL_id, nitrogen, chr, bin, start_pos, end_pos, total_sign_snps, tax_group) %>%
  unique()

unique(MAPL_traits$MAPL_id) # 622 MAPLs

MAPL_gns <- MAPL_genes %>%
  mutate(MAPL_id = paste0(nitrogen,"_chr",chr,"_",bin)) %>%
  select(MAPL_id, genes)


MAPLs <- MAPL_traits %>%
  dplyr::select(MAPL_id,  nitrogen, chr, start_pos, end_pos, total_sign_snps) %>%
  rename(nitrogen_treatment = nitrogen) %>%
  unique() %>%
  mutate(assoc_traits = NA, assoc_genes = NA)



for(i in 1:nrow(MAPLs)){
  mapl <- as.character(MAPLs[i, "MAPL_id"][1,1])
  
  traits <- MAPL_traits %>%
    dplyr::filter(MAPL_id == mapl) %>%
    arrange(desc(total_sign_snps))
  #print(traits)
  
  if(nrow(traits) > 0){
    str <- c()
    for(x in 1:nrow(traits)){
      s <- paste0(traits[x,"tax_group"],"(",traits[x,"total_sign_snps"],")")
      str <- c(str, s)
    } 
    str <- paste(str, collapse = ";")
    print(str)
    MAPLs[i, "assoc_traits"] <- str
  } else {
    MAPLs[i, "assoc_traits"] <- NA
  }
  
  genes <- MAPL_gns %>%
    dplyr::filter(MAPL_id == mapl) %>%
    arrange(genes)
  if(nrow(genes) > 0){
    gn <- paste(genes$genes, collapse = ";")
    print(gn)
    MAPLs[i, "assoc_genes"] <- gn
  } else {
    MAPLs[i, "assoc_genes"] <- NA
  }
}




signsnps <- all_gwas_signals %>%
  dplyr::select(chr, rs, ps, p_wald, tax_group, bin) %>%
  ungroup() %>%
  group_by(chr, bin) %>%
  summarize(lowest_p = min(p_wald, na.rm = TRUE), .groups = "keep")

most_sign_snp <- signsnps %>%
  left_join(all_gwas_signals) %>%
  filter(p_wald == lowest_p) %>%
  group_by(chr, bin) %>%
  sample_n(1) %>%
  ungroup() %>%
  mutate(MAPL_id = paste0(nitrogen,"_chr",chr,"_",bin)) %>%
  select(MAPL_id, rs, ps, p_wald, tax_group) %>%
  rename(most_sign_snp_rs = rs, most_sign_snp_pos = ps, most_sign_snp_p_wald = p_wald, trait_with_most_sign_snp = tax_group )



dataset5 <- left_join(MAPLs, most_sign_snp)



write_csv(dataset5, file = "data/Datasets/dataset5.csv", na = "")

