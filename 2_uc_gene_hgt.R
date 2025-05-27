## load libraries

library(tidyverse)

## load data

# metadata
load("uc_metadata_subset.RData")

# high quality MAGs, high quality genes, high quality gene clusters
load("HQ-gene-mapping-filtered.RData")

# true donor-recipient pairings
load("true_pairings_formatted.RData")

## load functions

# fmt hgt events
tentative_hgts_fmt <- function(HQ_genes){
  
  recipients <- 
    HQ_genes %>%
    filter(Group == "FMT") %>% 
    distinct(Participant_ID) %>% 
    pull(Participant_ID)
  
  donor_genes <- 
    HQ_genes %>%
    filter(Group == "Donor_individual") %>%
    rename(donor_gene = Gene_ID,
           donor_species = Species) %>%
    select(Cluster_ID, donor_gene, donor_species)
  
  tentative_hgt_events_fmt <- tibble()
  for (recipient in recipients) {
    
    # make a table of recipient gene clusters
    # only include necessary columns to make operations quicker
    recipient_genes <- 
      HQ_genes %>%
      filter(Participant_ID == recipient) %>%
      rename(recipient_gene = Gene_ID,
             recipient_species = Species) %>%
      select(Cluster_ID, recipient_gene, recipient_species, Timepoint)
    
    # find gene clusters present at baseline
    baseline_genes <- 
      recipient_genes %>%
      filter(Timepoint == "BL") %>%
      pull(Cluster_ID)
    
    # filter out gene clusters present at baseline
    recipient_genes_filtered <-
      recipient_genes %>%
      filter(!(Cluster_ID %in% baseline_genes),
             Timepoint == "wk8") %>%
      select(-Timepoint)
    
    # find joint gene clusters with donors
    # filter out cases where the donor and recipient species match
    tentative_hgt_events_fmt <-
      inner_join(donor_genes, 
                 recipient_genes_filtered,
                 relationship = "many-to-many") %>%
      filter(donor_species != recipient_species) %>% 
      bind_rows(tentative_hgt_events_fmt)
    
  }
  tentative_hgt_events_fmt
}

# placebo hgt events
tentative_hgts_placebo <- function(HQ_genes){
  
  recipients <- 
    HQ_genes %>%
    filter(Group == "Placebo") %>% 
    distinct(Participant_ID) %>% 
    pull(Participant_ID)
  
  donor_genes <- 
    HQ_genes %>%
    filter(Group == "Donor_individual") %>%
    rename(donor_gene = Gene_ID,
           donor_species = Species) %>%
    select(Cluster_ID, donor_gene, donor_species)
  
  tentative_hgt_events_fmt <- tibble()
  for (recipient in recipients) {
    
    # make a table of recipient gene clusters
    # only include necessary columns to make operations quicker
    recipient_genes <- 
      HQ_genes %>%
      filter(Participant_ID == recipient) %>%
      rename(recipient_gene = Gene_ID,
             recipient_species = Species) %>%
      select(Cluster_ID, recipient_gene, recipient_species, Timepoint)
    
    # find gene clusters present at baseline
    baseline_genes <- 
      recipient_genes %>%
      filter(Timepoint == "BL") %>%
      pull(Cluster_ID)
    
    # filter out gene clusters present at baseline
    recipient_genes_filtered <-
      recipient_genes %>%
      filter(!(Cluster_ID %in% baseline_genes),
             Timepoint == "wk8") %>%
      select(-Timepoint)
    
    # find joint gene clusters with donors
    # filter out cases where the donor and recipient species match
    tentative_hgt_events_fmt <-
      inner_join(donor_genes, 
                 recipient_genes_filtered,
                 relationship = "many-to-many") %>%
      filter(donor_species != recipient_species) %>% 
      bind_rows(tentative_hgt_events_fmt)
    
  }
  tentative_hgt_events_fmt
}

# randomise donor assignment for each placebo recipient to average HGT genes and clusters
randomise_donor_htgcs <- function(placebo_ids, donor_ids, donor_probs, uc_placebo_hgts_tmp){
  # placebo_ids: participant IDs of placebo recipients with HQ genes (chr vector)
  # donor_ids: participant IDs of individual donors with HQ genes (chr vector)
  # donor_probs: frequency of donor usage ( /1) across all FMT recipients
  # uc_placebo_hgts_tmp: output of tentative_hgts_placebo() - donor-specific HTGCs matching all donors
  
  # set seed
  set.seed(1234)
  
  # initialise empty tibble
  uc_placebo_hgts_randomised <- tibble()
  
  # loop 1000 times
  for(i in 1:1000){
    
    # sample 4-7 donors for each placebo recipient with HQ genes (n=14) based on the FMT donor usage frequency
    sampled_donors <- lapply(1:14, function(x) sample(x = donor_ids, sample(4:7, 1), replace = FALSE, prob = donor_probs)) # run sampling
    
    # format as dataframe with sampled donors for each placebo recipient
    dummy_placebo_donors <- data.frame( # format result into dataframe
      Recipient_ID = rep(placebo_ids, times = sapply(sampled_donors, length)),
      Donor_ID = unlist(sampled_donors)) %>%
      mutate(Pairings = paste0(Donor_ID, "_", Recipient_ID)) %>% # make pairings column
      select(Pairings)
    
    # get clean hgt events
    uc_placebo_hgts_randomised <- uc_placebo_hgts_tmp %>% # load temp HGT events object
      mutate(Recipient_sample_ID = str_split(recipient_gene, "_", simplify = TRUE)[,1]) %>% # extract recipient sample ID from gene ID
      mutate(Donor_sample_ID = str_split(donor_gene, "_", simplify = TRUE)[,1]) %>% # extract donor sample ID from gene ID
      # join recipient sample ID with metadata
      left_join(uc_metadata %>% select(-c(Timepoint, Group)), by = c("Recipient_sample_ID" = "Sample_ID")) %>% 
      rename(Recipient_ID = Participant_ID) %>% # get recipient IDs
      # join donor sample ID with metadata
      left_join(uc_metadata, by = c("Donor_sample_ID" = "Sample_ID")) %>% 
      rename(Donor_ID = Participant_ID) %>% # get donor IDs
      
      mutate(Pairings = paste0(Donor_ID, "_", Recipient_ID)) %>% # make pairings column
      inner_join(dummy_placebo_donors) %>% # join with dummy pairings data by pairing
      
      group_by(Recipient_ID) %>% # group by recipient ID
      summarise(hgt_genes_per_indiv = n_distinct(recipient_gene), # get number of genes involved in HGT for each individual
                hgt_clusters_per_indiv = n_distinct(Cluster_ID)) %>%
      bind_rows(uc_placebo_hgts_randomised) # bind with previous results
  }
  
  uc_placebo_hgts_randomised_mean <- uc_placebo_hgts_randomised %>% # overall binded data df
    group_by(Recipient_ID) %>% 
    summarise(mean_hgt_genes_per_indiv = mean(hgt_genes_per_indiv),
              mean_hgt_clusters_per_indiv = mean(hgt_clusters_per_indiv))
  
  uc_placebo_hgts_randomised_mean # return dataframe
}


## running

# filter the HQ MAG data to get just MAG ID and classification
mag_spp <- 
  HQ_mags %>%
  select(MAG_ID, classification)

# join the mag species data with the genes on high quality MAGs, extract sample ID and species data
HQ_genes_spp <- 
  HQ_genes %>%
  inner_join(mag_spp) %>%
  mutate(Species = str_split(classification, "s__", simplify = TRUE)[,2]) %>% # extract species of MAG
  mutate(across(where(is.character), ~ na_if(.,""))) %>% # replace empty cells (no species) with NA
  select(-c(classification, MAG_ID)) %>%
  mutate(Sample_ID = str_split(Gene_ID, "_", simplify = TRUE)[,1]) # extract sample ID from gene ID

# join metadata with the HQ gene and MAG data
HQ_genes_meta <-
  HQ_genes_spp %>%
  left_join(uc_metadata)

# filter out rows with missing species data and donor batch samples
HQ_genes_meta_spp <-
  HQ_genes_meta %>%
  filter(!is.na(Species)) %>%
  filter(Group != "Donor_batch")

# only include recipients with HQ genes at both baseline and week 8, individual donors with HQ genes
HQ_genes_samples <- HQ_genes_meta_spp %>%
  select(Sample_ID, Participant_ID, Group, Timepoint) %>%
  distinct() # 110 before subset

# donors
HQ_genes_donor_samples <- HQ_genes_samples %>%
  filter(Group == "Donor_individual") # 27/27

# FMT recipients 
HQ_genes_FMT_samples <- HQ_genes_samples %>%
  filter(Group == "FMT") %>%
  group_by(Participant_ID) %>%
  filter(all(c("BL", "wk8") %in% Timepoint)) %>%
  ungroup() # 18/32 participants @ wk8

# placebo recipients 
HQ_genes_placebo_samples <- HQ_genes_samples %>%
  filter(Group == "Placebo") %>%
  group_by(Participant_ID) %>%
  filter(all(c("BL", "wk8") %in% Timepoint)) %>%
  ungroup() # 14/21 participants @ wk8
  
# detail samples with missing HQ genes (not considering donor batches) for supplementary
HQ_genes_missing_samples <- uc_metadata %>%
  filter(Group != "Donor_batch") %>%
  anti_join(HQ_genes_samples) # 22

# participants to keep in analysis
HQ_genes_samples_keep <- bind_rows(HQ_genes_donor_samples, HQ_genes_FMT_samples, HQ_genes_placebo_samples)

# remove underscores from genera and species, and filter for recipients with both timepoints
HQ_genes_clean <- HQ_genes_meta_spp %>%
  mutate(Species_clean = str_replace_all(Species, "_.+?", "")) %>% # '?' for non-greedy matching
  inner_join(HQ_genes_samples_keep, by = c("Sample_ID", "Participant_ID", "Group", "Timepoint")) %>% # keep samples with needed timepoints
  select(-c(Species, Contig_ID, Sample_ID)) %>% # deselect columns
  rename(Species = Species_clean)
nrow(HQ_genes_clean) # 665578


## FMT cohort

# get hgt events (temporary - need to join with pairings data)
uc_fmt_hgts_tmp <- tentative_hgts_fmt(HQ_genes_clean)

# get FMT donor batch pairings
batch_pairings <- true_pairings_formatted %>%
  mutate(Pairings = paste0(Donor_ID, "_", Recipient_ID)) %>% # make pairings column
  select(Pairings)

# get clean hgt events
uc_fmt_hgts_clean <- uc_fmt_hgts_tmp %>% # load temp HGT events object
  mutate(Recipient_sample_ID = str_split(recipient_gene, "_", simplify = TRUE)[,1]) %>% # extract recipient sample ID from gene ID
  mutate(Donor_sample_ID = str_split(donor_gene, "_", simplify = TRUE)[,1]) %>% # extract donor sample ID from gene ID
  # join recipient sample ID with metadata
  left_join(uc_metadata %>% select(-c(Timepoint, Group)), by = c("Recipient_sample_ID" = "Sample_ID")) %>% 
  rename(Recipient_ID = Participant_ID) %>% # get recipient IDs
  # join donor sample ID with metadata
  left_join(uc_metadata %>% select(-c(Timepoint, Group)), by = c("Donor_sample_ID" = "Sample_ID")) %>% 
  rename(Donor_ID = Participant_ID) %>% # get donor IDs
  
  mutate(Pairings = paste0(Donor_ID, "_", Recipient_ID)) %>% # make pairings column
  inner_join(batch_pairings) # join with true pairings data by pairing

length(table(batch_pairings$Pairings)) # total pairings used in FMT batch treatment = 173
length(table(uc_fmt_hgts_clean$Recipient_ID)) # represented FMT recipients = 16/32

#nrow(uc_fmt_hgts_clean) # 8242 hgt events
#nrow(uc_fmt_hgts_clean %>% distinct(recipient_gene)) # 2707 distinct recipient genes involved

# count HQ genes, HQ clusters per individual
HQ_genes_fmt <- HQ_genes_clean %>% # load all HQ genes
  filter(Group == "FMT" & Timepoint == "wk8") # filter for FMT week 8 
nrow(HQ_genes_fmt) # 128232

HQ_genes_fmt_individual <- HQ_genes_fmt %>% # load HQ genes for FMT week 8 (note: HQ genes only available for 18/32 FMT recipients - see above)
  group_by(Participant_ID) %>% # group by participant ID
  summarise(HQ_genes_per_indiv = n(), # get number of HQ genes for each individual
            HQ_clusters_per_indiv = n_distinct(Cluster_ID)) # get number of distinct gene clusters for each individual

# find relative proportions of genes involved in hgt across each fmt recipient (week 8)
hgt_genes_fmt_individual <- uc_fmt_hgts_clean %>% # already joined with meta data
  group_by(Recipient_ID) %>% # group by recipient ID
  summarise(hgt_genes_per_indiv = n_distinct(recipient_gene), # get number of genes involved in HGT for each individual
            hgt_clusters_per_indiv = n_distinct(Cluster_ID)) %>% # get number of distinct gene clusters involved in HGT for each individual
  right_join(HQ_genes_fmt_individual, by = c("Recipient_ID" = "Participant_ID")) %>% # join with the "total" data from HQ genes above
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% # recipients with HQ genes but none involved in HGT
  mutate(pct_hgt_genes = hgt_genes_per_indiv/HQ_genes_per_indiv*100) %>% # calculate percentage of genes at week 8 involved in HGT
  mutate(pct_hgt_clusters = hgt_clusters_per_indiv/HQ_clusters_per_indiv*100) %>% # calculate percentage of gene clusters at week 8 involved in HGT
  left_join(uc_metadata, by = c("Recipient_ID" = "Participant_ID")) %>% # rejoin with meta data
  filter(Group == "FMT" & Timepoint == "wk8") %>% # filter for relevant data
  select(-Sample_ID) # deselect irrelevant columns

length(table(uc_fmt_hgts_clean$Pairings)) # represented pairings in HGT events = 67

# recipient 1037 has very high proportions of genes involved in HGT
HGT_events_1037 <- uc_fmt_hgts_clean %>% filter(Recipient_ID == "1037")
table(HGT_events_1037$Donor_ID)
#D031 D033 D043 D052 D054 D055 
#1773 3514   11    5    2 1818

HGT_events_1037_high_donors <- HGT_events_1037 %>%
  filter(Donor_ID %in% c("D031", "D033", "D055"))
table(HGT_events_1037_high_donors$donor_species) # mostly B. dorei
table(HGT_events_1037_high_donors$recipient_species) # mostly B. vulgatus
# note: possible misclassification of 1037 B. dorei as B. vulgatus?


## placebo cohort

# get hgt events (temporary - need to join with pairings data)
uc_placebo_hgts_tmp <- tentative_hgts_placebo(HQ_genes_clean)

# need to simulate random combinations of 4-7 FMT donors for each placebo recipient
range(true_pairings_formatted$all_donors) # 4-7 donors used for each FMT recipient

nrow(HQ_genes_placebo_samples %>% select(Participant_ID) %>% distinct()) # 14 placebo recipients have HQ genes at baseline and week 8
nrow(HQ_genes_donor_samples %>% select(Participant_ID) %>% distinct()) # 14 (all) donors had HQ genes

# specify probabilities for donors based on their use across all FMT recipients
donor_probabilities <- true_pairings_formatted %>%
  count(Donor_ID) %>% # count number of times each donor was used across all FMT participants
  mutate(prob = n/32) # donor use frequency = probability

donor_ids <- donor_probabilities %>% pull(Donor_ID) # donor IDs
donor_probs <- donor_probabilities %>% pull(prob) # frequency ( /1) of each donor's use across all FMT recipients in FOCUS trial

# placebo recipients with HQ genes (analysis requires HQ genes)
placebo_ids <- HQ_genes_placebo_samples %>% select(Participant_ID) %>% distinct() %>% pull(Participant_ID)

# randomise donor assignment for each placebo recipient to average HGT genes and clusters
uc_placebo_hgts_randomised_mean <- randomise_donor_htgcs(placebo_ids, donor_ids, donor_probs, uc_placebo_hgts_tmp)

# count HQ genes, HQ clusters per individual
HQ_genes_placebo <- HQ_genes_clean %>% # load all HQ genes
  filter(Group == "Placebo" & Timepoint == "wk8") # filter for placebo week 8 
nrow(HQ_genes_placebo) # 96283

HQ_genes_placebo_individual <- HQ_genes_placebo %>% # load HQ genes for placebo week 8 (note: HQ genes only available for 14/21 placebo recipients - see above)
  group_by(Participant_ID) %>% # group by participant ID
  summarise(HQ_genes_per_indiv = n(), # get number of HQ genes for each individual
            HQ_clusters_per_indiv = n_distinct(Cluster_ID)) # get number of distinct gene clusters for each individual

# find relative proportions of (mean randomised) genes involved in hgt across each placebo recipient (week 8)
hgt_genes_placebo_individual <- uc_placebo_hgts_randomised_mean %>% # already joined with meta data
  right_join(HQ_genes_placebo_individual, by = c("Recipient_ID" = "Participant_ID")) %>% # join with the "total" data from HQ genes above
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% # recipients with HQ genes but none involved in HGT
  mutate(pct_hgt_genes = mean_hgt_genes_per_indiv/HQ_genes_per_indiv*100) %>% # calculate percentage of genes at week 8 involved in HGT
  mutate(pct_hgt_clusters = mean_hgt_clusters_per_indiv/HQ_clusters_per_indiv*100) %>% # calculate percentage of gene clusters at week 8 involved in HGT
  left_join(uc_metadata, by = c("Recipient_ID" = "Participant_ID")) %>% # rejoin with meta data
  filter(Group == "Placebo" & Timepoint == "wk8") %>% # filter for relevant data
  select(-Sample_ID) # deselect irrelevant columns

# don't quantify hgt events and distinct recipient genes involved overall as it depends on donor sampling

# save data
#save(HQ_genes_clean, file = "HQ_genes_clean.RData") # FMT and placebo HQ genes
#save(uc_fmt_hgts_clean, file = "uc_fmt_hgts_clean.RData") # note: no placebo equivalent
#save(hgt_genes_fmt_individual, file = "hgt_genes_fmt_individual.RData")
#save(hgt_genes_placebo_individual, file = "hgt_genes_placebo_individual.RData")

# bind data for further analysis
all_hgt_genes_individual <- hgt_genes_placebo_individual %>%
  rename(hgt_genes_per_indiv = mean_hgt_genes_per_indiv, # rename placebo cols
         hgt_clusters_per_indiv = mean_hgt_clusters_per_indiv) %>%
  bind_rows(hgt_genes_fmt_individual) # bind rows
nrow(all_hgt_genes_individual) # 32

# plot percentage of gene clusters involved in HGT

# remove outliers
remove_outliers <- function(x) {
  x[x > quantile(x, .25) - 1.5*IQR(x) & x < quantile(x, .75) + 1.5*IQR(x)]
}

# plot
pct_htgcs_plot <- 
  all_hgt_genes_individual %>% # load combined data
  filter(pct_hgt_clusters %in% remove_outliers(pct_hgt_clusters)) %>%
  ggplot(aes(x = factor(Timepoint, labels = "Week 8"), y = pct_hgt_clusters, fill = Group)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 1.7) +
  geom_boxplot(alpha = 0.6, outlier.colour = NA) +
  scale_fill_manual(values = c("#b2182b", "#2166ac")) + 
  scale_colour_manual(values = c("#b2182b", "#2166ac"), name = NULL) +
  scale_y_continuous(limits = c(0, 2)) +
  xlab(NULL) + ylab("Percentage of gene clusters in HGT") + 
  theme_bw() +
  theme(#legend.position = "none",
        axis.text.x = element_text(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4)) 

# run statistics to find impact of FMT on percentage of HTGCs (increase/ decrease?) from placebo vs FMT distributions
# check normality
all_hgt_genes_individual %>% ggplot(aes(x= pct_hgt_clusters)) + geom_histogram() # not normal
pull(all_hgt_genes_individual, pct_hgt_clusters) %>% shapiro.test() # p = 1.005e-10 (significant) == not normal
# run wilcoxon test
all_box_data <- all_hgt_genes_individual %>% 
  spread(Group, pct_hgt_clusters)
wilcox.test(all_box_data$FMT, all_box_data$Placebo)  # p = 0.2961
