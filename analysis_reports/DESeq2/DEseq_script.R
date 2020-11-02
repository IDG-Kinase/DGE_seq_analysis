#!/usr/bin/env Rscript

library("tximeta")
library("DESeq2")

library("BiocParallel")

library(tictoc)
library(readr)
library(here)

register(MulticoreParam(10))

tic();

drug_perturb_exp = summarizeToGene(tximeta(read_csv(here('analysis_reports/DESeq2/CRISPRi_exp_info.csv'))))

# dds <- DESeqDataSet(drug_perturb_exp, design = ~group)
# keep <- rowSums(counts(dds)) > 1
# dds <- dds[keep,]
# 
# dds_analysis <- DESeq(dds)
# write_rds(dds_analysis, here('DESeq_results/DEseq_analysis.rds'))
toc();
