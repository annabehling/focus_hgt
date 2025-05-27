# focus_hgt

## The FOCUS Trial

The FOCUS Trial was a double-blinded randomised placebo-controlled trial that assessed the efficacy of FMT to treat adult (aged 18-75 years) ulcerative colitis (UC). FMT from 4-7 donors was administered via an initial colonoscopic infusion, followed by enemas at a frequency of 5 days per week for 8 weeks. Stool samples were obtained from individual donors and donor batches, as well as recipients every 4 weeks during the blinded treatment phase. A subset of 157 metagenomic sequencing files were selected for analysis.

Trial paper: https://doi.org/10.1016/S0140-6736(17)30182-4

Metagenomic data: https://doi.org/10.1053/j.gastro.2018.12.001

## Horizontal gene transfer analysis

This repository contains the following R scripts used to analyse metagenomic data from the FOCUS Trial for evidence of horizontal gene transfer (HGT).

- 1_uc_waafle_hgt.R : analysis of [WAAFLE](https://github.com/biobakery/waafle) HGT output
- 2_uc_gene_hgt.R : gene-based detection of HGT events
- 3_uc_gene_hgt_engrafted.R : identification of HGT from engrafted donor strains
- 4_uc_gb_hgt_comparison.R : comparison of HGT (WAAFLE) in the FOCUS Trial and Gut Bugs Trial
