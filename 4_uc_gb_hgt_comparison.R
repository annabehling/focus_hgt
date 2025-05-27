## code summary:
# compare the rate of HGT in FOCUS Trial and Gut Bugs Trial samples
# normalise HGT events by the number of contigs in each sample

## load libraries

library(tidyverse)
library(ggbeeswarm)
library(purrr)

## load data

# metadata
load("uc_metadata_subset.RData") # FOCUS
load("metadata_clean.RData") # Gut Bugs

# lists of dataframes containing of contigs for each sample (created from WAAFLE output)
# FOCUS
load("hgt_dfs_uc.RData") # contigs with HGT event
load("no_hgt_dfs_uc.RData") # contigs with no HGT event
load("unclassified_dfs_uc.RData") # unclassified contigs
# Gut Bugs
load("hgt_dfs.RData") # contigs with HGT event
load("no_hgt_dfs.RData") # contigs with no HGT event
load("unclassified_dfs.RData") # unclassified contigs

## running

## FOCUS Trial

sample_IDs <- names(hgt_dfs_uc) # extract sample IDs from waafle data

# get count of each type of contig for each sample
nrow_hgt_dfs_uc <- lapply(hgt_dfs_uc, nrow) # HGT
nrow_no_hgt_dfs_uc <- lapply(no_hgt_dfs_uc, nrow) # no HGT
nrow_unclassified_dfs_uc <- lapply(unclassified_dfs_uc, nrow) # unclassified
nrow_total_uc <- mapply(sum, nrow_hgt_dfs_uc, nrow_no_hgt_dfs_uc, nrow_unclassified_dfs_uc) # total

# make dataframes of contig counts
contig_count_data_uc <- data.frame(sample_IDs, unlist(nrow_hgt_dfs_uc), unlist(nrow_total_uc)) %>%
  left_join(uc_metadata, by = c("sample_IDs" = "Sample_ID")) %>% # join with metadata
  mutate(Trial = "FOCUS") %>% # mutate trial column
  rename(Sample_ID = sample_IDs,
         n_HGT_contigs = unlist.nrow_hgt_dfs_uc.,
         n_total_contigs = unlist.nrow_total_uc.) %>%
  mutate(pct_hgt_contigs = n_HGT_contigs/n_total_contigs*100)

## Gut Bugs Trial

Sample_ID <- reduced_meta %>% pull(Sample_ID) # extract sample IDs from matching metadata as Gut Bugs HGT dataframe was unnamed

# get count of each type of contig for each sample
nrow_hgt_dfs_gb <- lapply(hgt_dfs, nrow) # HGT
nrow_no_hgt_dfs_gb <- lapply(no_hgt_dfs, nrow) # no HGT
nrow_unclassified_dfs_gb <- lapply(unclassified_dfs, nrow) # unclassified
nrow_total_gb <- mapply(sum, nrow_hgt_dfs_gb, nrow_no_hgt_dfs_gb, nrow_unclassified_dfs_gb) # total

# make dataframes of contig counts
contig_count_data_gb <- data.frame(Sample_ID, unlist(nrow_hgt_dfs_gb), unlist(nrow_total_gb)) %>%
  left_join(reduced_meta, by = "Sample_ID") %>% # join with metadata
  mutate(Trial = "Gut Bugs") %>% # mutate trial column
  select(-Sex) %>% # data not in FOCUS Trial
  relocate(Group, .before = Timepoint) %>%
  rename(n_HGT_contigs = unlist.nrow_hgt_dfs_gb.,
         n_total_contigs = unlist.nrow_total_gb.) %>%
  mutate(pct_hgt_contigs = n_HGT_contigs/n_total_contigs*100)

## join data from two trials

all(colnames(contig_count_data_uc) == colnames(contig_count_data_gb)) # TRUE

contig_hgt_joined <- bind_rows(contig_count_data_uc, contig_count_data_gb) %>%
  filter(Group != "Donor_batch") %>% # remove donor batch samples
  filter(Timepoint != "Week 12" & Timepoint != "Week 26") %>% # remove additional timepoints from Gut Bugs
  mutate(Timepoint_general = case_when(Timepoint == "BL" & Group == "FMT" ~ "Pre-FMT", # make generalised timepoints
                                       Timepoint == "Baseline" & Group == "FMT" ~ "Pre-FMT",
                                       Timepoint == "wk8" & Group == "FMT" ~ "Post-FMT",
                                       Timepoint == "Week 6" & Group == "FMT" ~ "Post-FMT",
                                       Timepoint == "BL" & Group == "Placebo" ~ "Pre-placebo",
                                       Timepoint == "Baseline" & Group == "Placebo" ~ "Pre-placebo",
                                       Timepoint == "wk8" & Group == "Placebo" ~ "Post-placebo",
                                       Timepoint == "Week 6" & Group == "Placebo" ~ "Post-placebo",
                                       Timepoint == "Donor" ~ "Donor")) %>%
  group_by(Participant_ID, Timepoint_general, Trial) %>%
  mutate(mean_value = mean(pct_hgt_contigs)) %>% # get mean values for donors with multiple samples
  ungroup() %>%
  select(mean_value, Timepoint_general, Trial, Participant_ID) %>%
  distinct()
nrow(contig_hgt_joined) # 298 samples total

# reorder timepoint levels for plotting
contig_hgt_joined$Timepoint_general <- factor(contig_hgt_joined$Timepoint_general, levels = c("Donor", "Pre-FMT", "Post-FMT", "Pre-placebo", "Post-placebo"))

## plot data

contig_hgt_plot <- contig_hgt_joined %>%
  ggplot(aes(x = Timepoint_general, y = mean_value, fill = Trial)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 1.7) +
  geom_boxplot(alpha = 0.7, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) +
  ylim(0.0, 0.4) +
  xlab(NULL) + ylab("Normalised HGT rate (HGT events/contigs)") +
  theme_bw() +
  theme(axis.text.x = element_text(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

## run stats

# normality testing of data
contig_hgt_joined %>% filter(Trial == "FOCUS" & Timepoint_general == "Donor") %>% pull(mean_value) %>% shapiro.test() # 0.1436 normal
contig_hgt_joined %>% filter(Trial == "FOCUS" & Timepoint_general == "Pre-FMT") %>% pull(mean_value) %>% shapiro.test() # 0.0003093 not normal
# data for at least one cohort is not normal - run Wilcoxon test to compare all cohort pairs

# define cohorts for comparison
cohorts_to_compare <- c("Donor", "Pre-FMT", "Post-FMT", "Pre-placebo", "Post-placebo")

# loop through each pair of cohorts and run Wilcoxon test
test_results <- cohorts_to_compare %>%
  set_names() %>%
  map(function(cohort_name) {
    df_sub <- contig_hgt_joined %>%
      filter(Timepoint_general == cohort_name)
    wilcox.test(mean_value ~ Trial, data = df_sub)
  })

# extract p values
p_values <- sapply(test_results, function(x) x$p.value)

# fdr correction of p values using BH method
adjusted_p_values <- p.adjust(p_values, method = "fdr")

# summarise results
final_results <- data.frame(
  cohort = names(p_values),
  p_value = p_values,
  adj_p_value = adjusted_p_values
)
final_results

#                   cohort      p_value  adj_p_value
#Donor               Donor 8.568387e-03 8.568387e-03 **
#Pre-FMT           Pre-FMT 3.822668e-06 1.911334e-05 ***
#Post-FMT         Post-FMT 1.433340e-03 1.802753e-03 **
#Pre-placebo   Pre-placebo 5.500447e-05 1.375112e-04 ***
#Post-placebo Post-placebo 1.442202e-03 1.802753e-03 **