#!/usr/bin/env Rscript

library("tximeta")
library("DESeq2")

library("BiocParallel")

library(tictoc)
library(readr)
library(here)
library(tidyverse)
library(readxl)

register(MulticoreParam(128))

tic();

well_to_treatment = read_excel(here('Sequencing_results/MB_Array1_Layout.xlsx')) %>%
  dplyr::select(Dest384Well,Target,Vector) %>%
  mutate('well_letter' = str_extract(Dest384Well,'[[:alpha:]]'),
         'well_digit' = str_extract(Dest384Well,'[[:digit:]]{2}$')) %>%
  mutate('well_id' = paste0(well_letter,well_digit)) %>%
  dplyr::select(-Dest384Well,-well_letter,-well_digit) %>% 
  
  #Adding a column that indicates the repeat number
  group_by(Target, Vector) %>% 
  mutate(rep = 1:n()) %>% 
  ungroup()

#There was a problem with the barcode on well P11, so the reads from P11 are
#actually assigned to a well named p11_new1 in the naming scheme I used to save
#out the salmon results. Replacing that entry here:
P11_row = well_to_treatment %>% 
  filter(well_id == "P11")
stopifnot(dim(P11_row)[1] == 1)
P11_row$well_id = "p11_new1"
well_to_treatment = well_to_treatment %>%
  filter(well_id != "P11") %>%
  add_row(P11_row)

exp_info = data.frame(files = Sys.glob(here('Sequencing_results/combined_sequencing/salmon_alignment/*/quant.sf'))) %>%
  mutate(salmon_folder = basename(dirname(as.character(files)))) %>%
  mutate(well_id = str_extract(salmon_folder,"[^\\.]*")) %>%
  left_join(well_to_treatment) %>%
  mutate(names = paste0(Target,"_",Vector,"_",rep)) %>%
  #ensure that DMSO is the first level, forcing the fold change values to be
  #calculated as compound_treatment/DMSO
  mutate(Target = relevel(as.factor(Target), "DMSO","NonTarget"),
         rep = as.factor(rep)) %>%
  #Remove wells not mapped to a target
  filter(!is.na(Target)) %>% 
  mutate(group = paste0(Target,"_", Vector))

exp_info_CRISPRi = exp_info %>%
  #Remove the drug treatments
  filter(! Target %in% c("DMSO","Trametinib","JQ1", "Tram + JQ1")) %>%
  
  #Ensure that all the NonTarget runs get grouped for comparison to the CRISPRi
  #treatments
  mutate(group = ifelse(Target == "NonTarget", "NonTarget", group)) %>%
  mutate(group = relevel(as.factor(group), "NonTarget"),
         Target = relevel(as.factor(Target), "NonTarget")) %>%
  identity()

exp_info_compound = exp_info %>%
  filter(Target %in% c("DMSO","Trametinib","JQ1", "Tram + JQ1")) %>%
  mutate(Target = ifelse(Target == "Tram + JQ1","Tram_JQ1",as.character(Target))) %>%
  mutate(Target = relevel(as.factor(Target), "DMSO")) %>%
  identity()


drug_perturb_exp = summarizeToGene(tximeta(exp_info_CRISPRi))

# dds <- DESeqDataSet(drug_perturb_exp, design = ~group)
# keep <- rowSums(counts(dds)) > 1
# dds <- dds[keep,]
# 
# dds_analysis <- DESeq(dds)
# write_rds(dds_analysis, here('DESeq_results/DEseq_analysis.rds'))

dds <- DESeqDataSet(drug_perturb_exp, design = ~Target)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

dds_analysis <- DESeq(dds)
write_rds(dds_analysis, here('DESeq_results/DEseq_analysis_per_gene.rds'))

toc();
