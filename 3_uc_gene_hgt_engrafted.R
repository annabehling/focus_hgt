## load libraries

library(tidyverse)
library(ggalluvial)
library(chisq.posthoc.test)

## load data

# novel donor matching strains at week 8 from optimised strain engraftment analysis
load("sea_output_optimised.RData") 

# formatted HQ gene data from 0_uc_gene_hgt.R
load("HQ_genes_clean.RData")

# HGT events data from 0_uc_gene_hgt.R
load("uc_fmt_hgts_clean.RData")

# high quality MAGs, high quality genes, high quality gene clusters
load("HQ-gene-mapping-filtered.RData")

# metadata
load("uc_metadata_subset.RData") # uc_metadata

## load palettes

large_colour_palette_3 <- c("hotpink", "maroon", "orange", "pink2" ,"#1B9E77", "turquoise", "khaki2", "cornflowerblue", "#3B9AB2", "#EBCC2A",
                            "lightblue", "#BEAED4", "lightgreen", "slategray3")

## running

# identify single or multiple donor origin of novel donor-matching strains in FMT recipients at week 8
novel_strain_source <- sea_output_best %>%
  select(Species, Recipient_ID, Donor_ID) %>%
  group_by(Species, Recipient_ID) %>%
  summarise(n_donors = n_distinct(Donor_ID)) %>%
  ungroup() %>%
  mutate(source = case_when(n_donors == 1 ~ "single_donor", 
                            n_donors > 1 ~ "multi_donor")) %>%
  left_join(sea_output_best %>% select(Species, Recipient_ID, Donor_ID), by = c("Species", "Recipient_ID")) %>%
  distinct() %>% # removes strain matches to the same donor twice (two samples)
  mutate(Source = ifelse(grepl("multi_donor", source), "multiple_donors", Donor_ID)) %>%
  select(Species, Recipient_ID, Source)

# get engrafted donor hgt strains at week 8 that horizontally transferred genes to recipient species
donor_engrafted_strains <- novel_strain_source %>%
  filter(Source != "multiple_donors") %>% # specific donor strain
  filter(grepl("s__", Species)) %>% # filter for species
  mutate(Species = str_split(Species, "s__", simplify = TRUE)[,2]) %>% # remove 's__' prefix
  mutate(Species = str_replace_all(Species, "_", " ")) # replace '_' with ' '

uc_engrafted_hgts <-
  uc_fmt_hgts_clean %>%
  select(c(Cluster_ID, Donor_ID, donor_species, Recipient_ID, recipient_species)) %>%
  # join with engrafted strain data by Donor ID, Recipient ID and donor species
  inner_join(donor_engrafted_strains, by = c("Donor_ID" = "Source", "Recipient_ID", "donor_species" = "Species"))

# plot engraftment-specific HGT events
uc_engrafted_hgts_alluvial <- uc_engrafted_hgts %>%
  select(Donor_ID, Recipient_ID) %>%
  group_by(Donor_ID, Recipient_ID) %>% # group by donor / recipient pairing
  mutate(n_transfers = n()) %>% # calculate number of treatment-specific transfers between fmt donor / recipient pairings
  distinct() %>%
  ggplot(aes(axis1 = factor(Donor_ID), axis2 = factor(Recipient_ID), y = n_transfers)) +
  geom_alluvium(aes(fill = Donor_ID), alpha = 0.6) + 
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual("Donor ID", drop = FALSE, values = large_colour_palette_3[c(3, 8)]) +
  scale_x_continuous(breaks=c(1, 2), 
                     labels=c("Donor", "FMT recipient")) +
  ylab("HGT from engrafted species") +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# investigate species involved in engraftment-dependent HGT for each FMT donor-recipient pairing

D031_1037 <- uc_engrafted_hgts %>%
  filter(Donor_ID == "D031" & Recipient_ID == "1037")
table(D031_1037$donor_species) # 3 Bacteroides uniformis
table(D031_1037$recipient_species) # 3 Bacteroides vulgatus

D043_1038 <- uc_engrafted_hgts %>%
  filter(Donor_ID == "D043" & Recipient_ID == "1038")
table(D043_1038$donor_species) # 54 Bacteroides vulgatus
table(D043_1038$recipient_species) # 54 Prevotella copri

D043_1041 <- uc_engrafted_hgts %>%
  filter(Donor_ID == "D043" & Recipient_ID == "1041")
table(D043_1041$donor_species) # 15 Prevotella copri
table(D043_1041$recipient_species) # 15 Butyrivibrio crossotus 
  
D043_2219 <- uc_engrafted_hgts %>%
  filter(Donor_ID == "D043" & Recipient_ID == "2219")
table(D043_2219$donor_species) # 38 Barnesiella intestinihominis
table(D043_2219$recipient_species) # 38 Bacteroides plebeius


# identify inter-/ intra-phylum HGT events between engrafted donor species and host species

# make lookup table for clean (no underscores) species and phyla names
phylum_species <- HQ_mags %>%
  select(classification) %>%
  mutate(Species = str_split(classification, "s__", simplify = TRUE)[,2]) %>% # extract species of MAG
  mutate(Species_clean = str_replace_all(Species, "_.+?", "")) %>% # '?' for non-greedy matching
  select(-Species) %>% # deselect columns
  rename(Species = Species_clean) %>%
  mutate(tmp_Phylum = str_split(classification, "p__", simplify = TRUE)[,2]) %>%
  mutate(Phylum = str_split(tmp_Phylum, ";", simplify = TRUE)[,1]) %>%
  select(c(Species, Phylum)) %>%
  mutate(across(where(is.character), ~ na_if(.,""))) %>% # replace empty cells with NA
  filter(!is.na(Species) & !is.na(Phylum)) %>%
  distinct()

# add alternative phylum names
phylum_alt_species <- phylum_species %>%
  distinct(Phylum) %>%
  mutate(Phylum_alt = c("Firmicutes", "Firmicutes", "Bacteroidetes", "Proteobacteria", 
                        "Actinobacteria", "Verrucomicrobia", "Firmicutes")) %>% # entering manually based on uniprot synonyms
  right_join(phylum_species)

# plot phyla heatmap
engrafted_phyla_hgt <- uc_engrafted_hgts %>% # engrafted fmt-specific hgt data
  group_by(donor_species, recipient_species) %>%
  count() %>% # count combinations of donor and recipient species
  left_join(phylum_alt_species, by = c("donor_species" = "Species")) %>% # join to get phylum names for donor species
  rename("donor_phylum" = "Phylum_alt") %>%
  select(-Phylum) %>%
  left_join(phylum_alt_species, by = c("recipient_species" = "Species")) %>% # join to get phylum names for recipient species
  rename("recipient_phylum" = "Phylum_alt") %>%
  select(-Phylum)

engrafted_phyla_hgt_plot <- engrafted_phyla_hgt %>%
  ggplot(aes(x = donor_species, y = recipient_species, fill = n)) +
  geom_tile() +
  facet_grid (recipient_phylum~ donor_phylum, scales = "free", space = "free") + # rows ~ columns
  xlab("Engrafted donor species") + ylab("Recipient species") +
  scale_y_discrete(limits = rev) + # descending alphabetical order
  geom_vline(xintercept = seq(1.5, length(unique(engrafted_phyla_hgt$donor_species))-0.5, 1), linewidth = 1, colour = "grey") + # add vertical lines
  geom_hline(yintercept = seq(1.5, length(unique(engrafted_phyla_hgt$recipient_species))-0.5, 1), linewidth = 1, colour = "grey") + # add horizontal lines
  scale_fill_gradient(name = "Number of\nHGT events", low = "#f3a5ae", high = "#b2182b") +
  theme_bw() +
  theme(strip.text.x = element_text(angle = 90), # rotate x axis facel labels vertical
        strip.text.y = element_text(angle = 0), # rotate y axis facet labels horizontal 
        axis.text.x = element_text(angle = 45, hjust=1), # rotate x axis labels 45 degrees
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
engrafted_phyla_hgt_sum <- engrafted_phyla_hgt %>% group_by(donor_phylum) %>% summarise(sum(n)) # quantify HGT donors

# find relative abundance of background phyla (phyla represented by HQ donor genes)
HQ_donor_genes_phyla <- HQ_genes_clean %>%
  left_join(phylum_alt_species) %>% # join to get phylum names for species
  filter(Timepoint == "Donor") %>% # donor phyla
  group_by(Phylum_alt) %>% # group data
  summarise(Freq = n()) %>% # find frequencies of phyla
  mutate(Rel_abundance = Freq/sum(Freq)) %>% # find relative abundance
  mutate(Type = "background donor phyla") # make a type column

# find relative abundance of phyla in engrafted FMT donor strains with HGT
stackedbar_engrafted_hgt_data <- uc_engrafted_hgts %>% # load engrafted HGT data
  select(donor_species) %>%
  left_join(phylum_alt_species, by = c("donor_species" = "Species")) %>% # join with species:phylum lookup table
  group_by(Phylum_alt) %>% # group by phylum
  summarise(Freq = n()) %>% # find frequencies of phyla
  mutate(Rel_abundance = Freq/sum(Freq)) %>% # find relative abundance
  mutate(Type = "engrafted phyla with HGT") # make a type column

# find rel abundance of phyla in engrafted FMT donor strains
stackedbar_engrafted_all_data <- donor_engrafted_strains %>% # load week 8 engrafted strains in FMT recipients
  #left_join(uc_metadata) %>% # join with meta data
  select(Species) %>%
  left_join(phylum_alt_species, by = c("Species" = "Species")) %>% # join with species:phylum lookup table
  # manually include missing phylum data, looking up on UniProt taxonomy
  mutate(new_Phylum_alt = case_when(grepl("Bacteroides", Species) ~ "Bacteroidetes",
                                    grepl("Ruminococcus", Species) ~ "Firmicutes", # Bacillota
                                    grepl("Eubacterium", Species) ~ "Firmicutes", # Bacillota
                                    grepl("Parabacteroides", Species) ~ "Bacteroidetes" # Bacteroidota
  )) %>% 
  mutate(Phylum_alt_replacement = coalesce(Phylum_alt, new_Phylum_alt)) %>%
  #filter(is.na(Phylum_alt_replacement)) %>% # check remaining NAs
  select(-c(Phylum_alt, new_Phylum_alt)) %>%
  rename("Phylum_alt" = "Phylum_alt_replacement") %>%
  # format
  group_by(Phylum_alt) %>% # group by phylum
  summarise(Freq = n()) %>% # find frequencies of phyla
  mutate(Rel_abundance = Freq/sum(Freq)) %>% # find relative abundance
  mutate(Type = "engrafted phyla") %>% # make a type column
  bind_rows(stackedbar_engrafted_hgt_data)%>% # bind with engrafted HGT data
  bind_rows(HQ_donor_genes_phyla) # bind with all HQ donor genes phyla

# sort colouring for stacked bar plot
phylum_hexcode <- phylum_alt_species %>%
  select(Phylum_alt) %>%
  distinct() %>%
  mutate(Hexcode = c("#44AA99", "#88CCEE", "#cc7fbf", "#DDCC77", "#6699CC"))

all_engrafted_phyla_meta_alt <- stackedbar_engrafted_all_data %>%
  left_join(phylum_hexcode)

unique_all_phyla <- unique(all_engrafted_phyla_meta_alt$Phylum_alt)
my_phylum_hexcodes <- phylum_hexcode[phylum_hexcode$Phylum_alt %in% unique_all_phyla, ]
my_phylum_hexcodes_distinct <- my_phylum_hexcodes %>% arrange(Phylum_alt) %>% distinct()

# plot stacked bar
hgt_phyla_stacked_plot <- stackedbar_engrafted_all_data %>%
  ggplot(aes(x = Type, y = Rel_abundance, fill = Phylum_alt)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual("Phylum", values = my_phylum_hexcodes_distinct$Hexcode, 
                    labels = my_phylum_hexcodes_distinct$Phylum_alt) +
  xlab(NULL) + ylab("Phylum relative abundance") + 
  scale_x_discrete(labels = c("Background", "Engrafted", "Engrafted\nwith HGT")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4)) 

# statistics to compare relative abundances of phyla in 'engrafted' and 'engrafted with HGT' distributions
uc_engrafted_hgt_data <- 
  stackedbar_engrafted_all_data %>% 
  select(-Rel_abundance) %>% 
  filter(str_detect(Type, "engrafted")) %>% 
  pivot_wider(names_from = Type, values_from = Freq, values_fill = 0) %>% 
  ungroup() 

uc_names <- 
  uc_engrafted_hgt_data %>% 
  select(Phylum_alt)%>% 
  mutate(Dimension = as.character(row_number()))

uc_chisq <- 
  uc_engrafted_hgt_data %>% 
  select(-Phylum_alt) 

chisq.test(uc_chisq)  #p-value = 4.019e-09

uc_results <- chisq.posthoc.test(uc_chisq, method = "BH", round = 17)

uc_results <- 
  left_join(uc_names, uc_results)
uc_pvalues <- 
  uc_results %>% 
  filter(Value == "p values", if_any(.cols = everything(), ~ . <= 0.05))
