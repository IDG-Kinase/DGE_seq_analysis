# DGE-seq Analysis

This repo has code related to processing and quantifying the multiplexed DGE-seq data related to the CRISPRi and kinase inhibitor treatments. The processing pipeline involves a few simple scripts that take the sequencing results and split them into the sequencing reads from each well (see split_FASTQ_file.pl). Then takes each of the split sets of reads and gathers transcript read counts using salmon (see run_salmon_over_split.pl). Both of these scripts are a little bit rough, but have worked through this first pass at the data.

As for analyzing the results, I've put together three reports:

 * [DESeq2 Based Analysis](https://github.com/IDG-Kinase/DGE_seq_analysis/blob/master/analysis_reports/DESeq2/DESeq2_analysis.md) - Using the R package DESeq2, I've compared the reads from each of the CRISPRi targets to the non-target wells.
 * [Transcriptome Correlation](https://github.com/IDG-Kinase/DGE_seq_analysis/blob/master/analysis_reports/Correlation/DGE_seq_correlation.md) - This report covers calculating the correlation between the transcriptome in each of the wells.
 * [Bulk RNAseq Analysis](https://github.com/IDG-Kinase/DGE_seq_analysis/blob/master/analysis_reports/Bulk_compound_treatment/Bulk_DEseq2_analysis.md) - There is also some previously collected data related to the drug treatments tried in the DGE-seq, so I've reprocessed that sequencing data to produce the same type of report.
