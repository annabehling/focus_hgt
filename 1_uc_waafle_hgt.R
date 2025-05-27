## load libraries

library(tidyverse)
library(ggbeeswarm)
library(rstatix)
library(DescTools)
library(lme4)
library(lmerTest)
library(openxlsx)

## load functions

get_contigs <- function(results_dir, pattern){
  ## returns a list of dataframes (one dataframe per sample) containing contigs with classification matching 'pattern'
  # results_dir : path to directory containing WAAFLE output files (str)
  # pattern : file extension matching contig classification (".lgt.tsv" OR ".no_lgt.tsv" OR ".unclassified.tsv") (str)
  fnames <- list.files(path = results_dir, pattern = paste0("\\", pattern)) #list file names
  root <- sub("\\..+", "", fnames) #get file basename
  file_paths <- paste0(results_dir, root, pattern) #get full file path
  df_list <- lapply(file_paths,
                    FUN = function(files) {
                      read.table(files, header = TRUE, sep = "\t")
                    }) #read in files and make list of dataframes - one df per sample
  names(df_list) <- root # name list of dataframes, can extract later
  
  df_list # return
}

## format waafle output

# make lists of dataframes containing contigs with HGT / no HGT / unclassified
#hgt_dfs_uc <- get_contigs("waafle_data/", ".lgt.tsv")
#no_hgt_dfs_uc <- get_contigs("waafle_data/", ".no_lgt.tsv")
#unclassified_dfs_uc <- get_contigs("waafle_data/", ".unclassified.tsv")

# save files
#save(hgt_dfs_uc, file = "hgt_dfs_uc.RData")
#save(no_hgt_dfs_uc, file = "no_hgt_dfs_uc.RData")
#save(unclassified_dfs_uc, file = "unclassified_dfs_uc.RData")

## load data

# lists of dataframes containing of contigs for each sample (created from WAAFLE output)
load("hgt_dfs_uc.RData") # contigs with HGT event
load("no_hgt_dfs_uc.RData") # contigs with no HGT event
load("unclassified_dfs_uc.RData") # unclassified contigs

load("uc_metadata_subset.RData") # uc_metadata

load("metaphlan.Rdata") # metaphlan species data

## load palettes

palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "slategrey", "#CC79A7", "darkslateblue", "#D55E00")

## running

# check all names lists same
all(names(hgt_dfs_uc) == names(no_hgt_dfs_uc)) # TRUE
all(names(hgt_dfs_uc) == names(unclassified_dfs_uc)) # TRUE

# extract sample IDs
sample_IDs <- names(hgt_dfs_uc)
all(sample_IDs == uc_metadata %>% pull(Sample_ID)) # TRUE

# get count of each type of contig for each sample
nrow_hgt_dfs_uc <- lapply(hgt_dfs_uc, nrow) # HGT
nrow_no_hgt_dfs_uc <- lapply(no_hgt_dfs_uc, nrow) # no HGT
nrow_unclassified_dfs_uc <- lapply(unclassified_dfs_uc, nrow) # unclassified
nrow_total_uc <- mapply(sum, nrow_hgt_dfs_uc, nrow_no_hgt_dfs_uc, nrow_unclassified_dfs_uc) # total

# get percentage of each type of contig for each sample
pct_hgt_uc <- unlist(nrow_hgt_dfs_uc)/nrow_total_uc*100 # HGT
pct_no_hgt_uc <- unlist(nrow_no_hgt_dfs_uc)/nrow_total_uc*100 # no HGT
pct_unclassified_uc <- unlist(nrow_unclassified_dfs_uc)/nrow_total_uc*100 # unclassified

# make dataframes of contig counts
contig_count_data <- data.frame(sample_IDs, unlist(nrow_hgt_dfs_uc), unlist(nrow_no_hgt_dfs_uc), unlist(nrow_unclassified_dfs_uc))
contig_count_data_individuals <- contig_count_data %>%
  left_join(uc_metadata, by = c("sample_IDs" = "Sample_ID")) %>% # join with metadata
  filter(Group != "Donor_batch") %>% # remove donor batch samples
  select(-c(Participant_ID, Group, Timepoint)) # deselect unneeded columns
colnames(contig_count_data_individuals) <- c("Sample_ID", "n_HGT_contigs", "n_no_HGT_contigs", "n_unclassified_contigs")

# make dataframes of contig percentages
contig_pct_data <- data.frame(sample_IDs, pct_hgt_uc, pct_no_hgt_uc, pct_unclassified_uc)
contig_pct_data_individuals <- contig_pct_data %>%
  left_join(uc_metadata, by = c("sample_IDs" = "Sample_ID")) %>% # join with metadata
  filter(Group != "Donor_batch") %>% # remove donor batch samples
  select(-c(Participant_ID, Group, Timepoint)) # deselect unneeded columns
colnames(contig_pct_data_individuals) <- c("Sample_ID", "pct_HGT_contigs", "pct_no_HGT_contigs", "pct_unclassified_contigs")

# count individual samples
individual_samples <- uc_metadata %>%
  filter(Group != "Donor_batch") # remove donor batch samples
nrow(individual_samples) # 132
nrow(contig_count_data_individuals) # 132
nrow(contig_pct_data_individuals) # 132

nrow(uc_metadata %>% filter(Group == "FMT")) # 64 total FMT samples
nrow(uc_metadata %>% filter(Group == "Placebo")) # 41 total placebo samples

# get mean, min, max percentage values for results (individual donors and recipients only)
# HGT contigs
mean(contig_pct_data_individuals$pct_HGT_contigs) #0.1192644; gutbugs:0.1528122
min(contig_pct_data_individuals$pct_HGT_contigs) #0; gutbugs:0.07412074
max(contig_pct_data_individuals$pct_HGT_contigs) #0.3496226; gutbugs:0.2556661
# no HGT contigs
mean(contig_pct_data_individuals$pct_no_HGT_contigs) #71.14581; gutbugs:58.6406
min(contig_pct_data_individuals$pct_no_HGT_contigs) #50; gutbugs:38.1893
max(contig_pct_data_individuals$pct_no_HGT_contigs) #99.08257; gutbugs:80.33072
# unclassified contigs
mean(contig_pct_data_individuals$pct_unclassified_contigs) #28.73493; gutbugs:41.20659
min(contig_pct_data_individuals$pct_unclassified_contigs) #0.9174312; gutbugs:19.54895
max(contig_pct_data_individuals$pct_unclassified_contigs) #50; gutbugs:61.72985

# plot contig percentages for each sample
contig_pct_plot <- contig_pct_data_individuals %>% # already filtered for only recipient or individual donor data
  pivot_longer(cols = `pct_HGT_contigs`:`pct_unclassified_contigs`, # reshape data to long format
               names_to = "Contig_type",
               values_to = "Percentage_value") %>%
  ggplot(aes(x = factor(Contig_type, level = c("pct_HGT_contigs", "pct_no_HGT_contigs", "pct_unclassified_contigs")), 
             y = Percentage_value, 
             fill = Contig_type)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 1.7) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA) +
  scale_fill_manual(values = palette, name = NULL) +
  scale_colour_manual(values = palette, name = NULL) +
  scale_y_log10() +
  scale_x_discrete(labels=c("pct_HGT_contigs" = "With HGT", "pct_no_HGT_contigs" = "Without HGT", "pct_unclassified_contigs" = "Unclassified")) +
  xlab("Contig type") + ylab("Percentage of sample contigs (log base 10)") + 
  theme_bw() +
  theme(axis.text.x = element_text(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")
# warning message: Transformation introduced infinite values in continuous y-axis (will be due to 0 values for some HGT contigs) 


# compare FMT and placebo HGT events
metaphlan_spp <- rowSums(species != 0) # number of species in each sample (species richness)
hgt_contigs_uc <- unlist(nrow_hgt_dfs_uc) # number of contigs in each sample
hgt_cont_per_spp_uc <- hgt_contigs_uc/metaphlan_spp # number of contigs per species
richness_df_uc <- data.frame(sample_IDs, hgt_cont_per_spp_uc) # number of contigs per species for each sample
names(richness_df_uc)[1] <- "Sample_ID" # name sample ID column
richness_meta_uc <- inner_join(richness_df_uc, uc_metadata, by = "Sample_ID") # merge with metadata

# plot HGT events per species in each recipient sample
normalised_hgt_plot_donor_mean <- richness_meta_uc %>% # take raw data (all samples including donor batches)
  filter(Group == "FMT" | Group == "Placebo") %>%
  ggplot(aes(x = factor(Timepoint, levels = c("BL", "wk8"), labels = c("Baseline", "Week 8")), 
             y = hgt_cont_per_spp_uc, fill = Group)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 1.7) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c( "#b2182b", "#2166ac"), labels = c("FMT", "Placebo")) +
  scale_colour_manual(values = c("#b2182b", "#2166ac"), name = NULL) +
  scale_y_continuous(limits = c(0, 3)) +
  xlab(NULL) + ylab("HGT events/species") + 
  theme_bw() +
  theme(axis.text.x = element_text(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# prepare data for linear mixed effect model
full_dataset_uc <- richness_meta_uc %>% # full dataset (include all participants)
  filter(Group == "FMT" | Group == "Placebo") %>% # select only recipient data
  select(c(Participant_ID, Timepoint, hgt_cont_per_spp_uc)) %>% # select relevant columns
  reshape(idvar = "Participant_ID", timevar = "Timepoint", direction = "wide") %>% # make wide format
  inner_join(richness_meta_uc %>% select(c(Participant_ID, Group))) %>% # join with selected metadata
  distinct() %>% # select only distinct rows (duplicates after inner_join)
  rename("Baseline" = "hgt_cont_per_spp_uc.BL", "Week 8" = "hgt_cont_per_spp_uc.wk8") %>% # rename columns
  select(c(Participant_ID, Group, Baseline, `Week 8`)) %>%
  gather(key = "Timepoint", value = "hgt_cont_per_spp_uc", Baseline, `Week 8`) %>% # gather data back to long format
  mutate(across(c(Participant_ID, Group, Timepoint), as.factor)) %>%
  na.omit() # remove samples with missing data
full_dataset_uc$Group <- factor(full_dataset_uc$Group, levels = c("Placebo", "FMT")) # reorder Group levels

# fit full model
full_lmer <- lmer(hgt_cont_per_spp_uc ~ Group * Timepoint+(1|Participant_ID),data = full_dataset_uc) #full dataset ok for lmm
summary(full_lmer)

# fit reduced model (no interactions between fixed effects)
reduced_lmer <- lmer(hgt_cont_per_spp_uc ~ Group + Timepoint+(1|Participant_ID),data =full_dataset_uc)
summary(reduced_lmer)
# Group p = 0.244
# Timepoint p = 0.226

# run anova to compare full and reduced models
anova(full_lmer, reduced_lmer) # p = 0.7917
# reduced model fits better - non-significant p value and lower AIC value

# get confidence interval of reduced model fixed effect coefficients
confint(reduced_lmer)

## supplementary file
suppl_file_uc <- uc_metadata %>%
  inner_join(contig_count_data_individuals) %>% # contig counts
  rename(`Number of contigs with HGT` = n_HGT_contigs,
         `Number of contigs without HGT` = n_no_HGT_contigs,
         `Number of unclassified contigs` = n_unclassified_contigs)
#write.xlsx(suppl_file_uc, "hgt_contig_counts.xlsx")
