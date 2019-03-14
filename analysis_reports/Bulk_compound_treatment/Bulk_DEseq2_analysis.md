Bulk Compound Treatment DESeq2
================
Matthew Berginski

``` r
salmon_files = list.files(path = '../../Sequencing_results/Bulk_RNAseq_drug_treatments/salmon_alignment/',
                          pattern='quant.sf',all.files = T,recursive = T, full.names = T)

samples = tibble(Treatment = c(rep("Trametinib",3), 
                               rep("Trametinib_JQ1",3),
                               rep("JQ1",3),
                               rep("DMSO",3)),
                 files=salmon_files)

ENST_to_HGNC = read_csv(here('ENST_to_HGNC.csv'))
ENST_to_ENST = ENST_to_HGNC
ENST_to_ENST$hgnc_symbol = ENST_to_ENST$ensembl_transcript_id_version

Bulk_seq_hits = c();
for (this_drug in c("Trametinib","Trametinib_JQ1","JQ1")) {
  
  these_samples = samples %>% filter(Treatment %in% c("DMSO",this_drug))

  salmon_import = tximport(these_samples$files,type="salmon",tx2gene=ENST_to_ENST)
  
  ddsTxi <- DESeqDataSetFromTximport(salmon_import,
                                     colData = these_samples,
                                     design = ~ Treatment)
  
  dds <- DESeq(ddsTxi)
  non_target_vs_results <- results(dds)
  
  hits <- as.data.frame(non_target_vs_results[which(non_target_vs_results$padj < 0.05),]) %>%
    rownames_to_column(var = 'ensembl_transcript_id_version') %>%
    left_join(ENST_to_HGNC) %>%
    mutate(treatment = this_drug) %>%
    select(hgnc_symbol,treatment,log2FoldChange,ensembl_transcript_id_version,everything())
  
  Bulk_seq_hits = rbind(Bulk_seq_hits,hits)
}
write_csv(Bulk_seq_hits,here('Bulk_RNAseq_DESeq2_hits.csv'))
```

``` r
hit_counts = Bulk_seq_hits %>%
  group_by(treatment) %>%
  summarise(hit_count = n())

kable(hit_counts)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
treatment
</th>
<th style="text-align:right;">
hit\_count
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
JQ1
</td>
<td style="text-align:right;">
1987
</td>
</tr>
<tr>
<td style="text-align:left;">
Trametinib
</td>
<td style="text-align:right;">
7370
</td>
</tr>
<tr>
<td style="text-align:left;">
Trametinib\_JQ1
</td>
<td style="text-align:right;">
10121
</td>
</tr>
</tbody>
</table>
.

Overlaps with DGE-seq
=====================

``` r
DGE_seq_hits = read_csv(here('DESeq2_CRISPR_drug_hits_combined.csv'))
```

    ## Parsed with column specification:
    ## cols(
    ##   hgnc_symbol = col_character(),
    ##   compound = col_character(),
    ##   log2FoldChange = col_double(),
    ##   ensembl_transcript_id_version = col_character(),
    ##   baseMean = col_double(),
    ##   lfcSE = col_double(),
    ##   stat = col_double(),
    ##   pvalue = col_double(),
    ##   padj = col_double()
    ## )

``` r
Trametinib_bulk = Bulk_seq_hits %>% filter(treatment == "Trametinib")
Trametinib_DGE = DGE_seq_hits %>% filter(compound == "Trametinib")
Trametinib_overlap = Trametinib_DGE %>% 
  filter(hgnc_symbol %in% Trametinib_bulk$hgnc_symbol) %>%
  select(hgnc_symbol,compound,log2FoldChange,ensembl_transcript_id_version,everything())

kable(Trametinib_overlap)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
hgnc\_symbol
</th>
<th style="text-align:left;">
compound
</th>
<th style="text-align:right;">
log2FoldChange
</th>
<th style="text-align:left;">
ensembl\_transcript\_id\_version
</th>
<th style="text-align:right;">
baseMean
</th>
<th style="text-align:right;">
lfcSE
</th>
<th style="text-align:right;">
stat
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
padj
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
CBX3
</td>
<td style="text-align:left;">
Trametinib
</td>
<td style="text-align:right;">
-22.087521
</td>
<td style="text-align:left;">
ENST00000337620.8
</td>
<td style="text-align:right;">
24.04756
</td>
<td style="text-align:right;">
3.193074
</td>
<td style="text-align:right;">
-6.917323
</td>
<td style="text-align:right;">
0.00e+00
</td>
<td style="text-align:right;">
0.0000000
</td>
</tr>
<tr>
<td style="text-align:left;">
CNBP
</td>
<td style="text-align:left;">
Trametinib
</td>
<td style="text-align:right;">
-9.508608
</td>
<td style="text-align:left;">
ENST00000500450.6
</td>
<td style="text-align:right;">
68.29618
</td>
<td style="text-align:right;">
2.285941
</td>
<td style="text-align:right;">
-4.159603
</td>
<td style="text-align:right;">
3.19e-05
</td>
<td style="text-align:right;">
0.0325381
</td>
</tr>
<tr>
<td style="text-align:left;">
RSL1D1
</td>
<td style="text-align:left;">
Trametinib
</td>
<td style="text-align:right;">
-22.355768
</td>
<td style="text-align:left;">
ENST00000573618.5
</td>
<td style="text-align:right;">
29.64767
</td>
<td style="text-align:right;">
3.191343
</td>
<td style="text-align:right;">
-7.005128
</td>
<td style="text-align:right;">
0.00e+00
</td>
<td style="text-align:right;">
0.0000000
</td>
</tr>
<tr>
<td style="text-align:left;">
GTF2I
</td>
<td style="text-align:left;">
Trametinib
</td>
<td style="text-align:right;">
-22.293531
</td>
<td style="text-align:left;">
ENST00000621040.4
</td>
<td style="text-align:right;">
28.00058
</td>
<td style="text-align:right;">
3.191485
</td>
<td style="text-align:right;">
-6.985317
</td>
<td style="text-align:right;">
0.00e+00
</td>
<td style="text-align:right;">
0.0000000
</td>
</tr>
</tbody>
</table>
``` r
JQ1_bulk = Bulk_seq_hits %>% filter(treatment == "JQ1")
JQ1_DGE = DGE_seq_hits %>% filter(compound == "JQ1")
JQ1_overlap = JQ1_DGE %>% 
  filter(hgnc_symbol %in% JQ1_bulk$hgnc_symbol) %>%
  select(hgnc_symbol,compound,log2FoldChange,ensembl_transcript_id_version,everything())

kable(JQ1_overlap)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
hgnc\_symbol
</th>
<th style="text-align:left;">
compound
</th>
<th style="text-align:right;">
log2FoldChange
</th>
<th style="text-align:left;">
ensembl\_transcript\_id\_version
</th>
<th style="text-align:right;">
baseMean
</th>
<th style="text-align:right;">
lfcSE
</th>
<th style="text-align:right;">
stat
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
padj
</th>
</tr>
</thead>
<tbody>
<tr>
</tr>
</tbody>
</table>
``` r
Trametinib_JQ1_bulk = Bulk_seq_hits %>% filter(treatment == "Trametinib_JQ1")
Trametinib_JQ1_DGE = DGE_seq_hits %>% filter(compound == "Tram + JQ1")
Trametinib_JQ1_overlap = Trametinib_JQ1_DGE %>% 
  filter(hgnc_symbol %in% Trametinib_JQ1_bulk$hgnc_symbol) %>%
  select(hgnc_symbol,compound,log2FoldChange,ensembl_transcript_id_version,everything())

kable(Trametinib_JQ1_overlap)
```

<table>
<thead>
<tr>
<th style="text-align:left;">
hgnc\_symbol
</th>
<th style="text-align:left;">
compound
</th>
<th style="text-align:right;">
log2FoldChange
</th>
<th style="text-align:left;">
ensembl\_transcript\_id\_version
</th>
<th style="text-align:right;">
baseMean
</th>
<th style="text-align:right;">
lfcSE
</th>
<th style="text-align:right;">
stat
</th>
<th style="text-align:right;">
pvalue
</th>
<th style="text-align:right;">
padj
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
NDUFC1
</td>
<td style="text-align:left;">
Tram + JQ1
</td>
<td style="text-align:right;">
-10.728116
</td>
<td style="text-align:left;">
ENST00000394225.6
</td>
<td style="text-align:right;">
161.73668
</td>
<td style="text-align:right;">
2.300238
</td>
<td style="text-align:right;">
-4.663916
</td>
<td style="text-align:right;">
3.10e-06
</td>
<td style="text-align:right;">
0.0046620
</td>
</tr>
<tr>
<td style="text-align:left;">
NDUFC1
</td>
<td style="text-align:left;">
Tram + JQ1
</td>
<td style="text-align:right;">
9.869606
</td>
<td style="text-align:left;">
ENST00000503453.5
</td>
<td style="text-align:right;">
83.24251
</td>
<td style="text-align:right;">
2.318026
</td>
<td style="text-align:right;">
4.257763
</td>
<td style="text-align:right;">
2.06e-05
</td>
<td style="text-align:right;">
0.0206850
</td>
</tr>
<tr>
<td style="text-align:left;">
CLU
</td>
<td style="text-align:left;">
Tram + JQ1
</td>
<td style="text-align:right;">
10.112438
</td>
<td style="text-align:left;">
ENST00000523500.5
</td>
<td style="text-align:right;">
98.21545
</td>
<td style="text-align:right;">
2.311465
</td>
<td style="text-align:right;">
4.374904
</td>
<td style="text-align:right;">
1.21e-05
</td>
<td style="text-align:right;">
0.0136914
</td>
</tr>
<tr>
<td style="text-align:left;">
TXNRD1
</td>
<td style="text-align:left;">
Tram + JQ1
</td>
<td style="text-align:right;">
10.394335
</td>
<td style="text-align:left;">
ENST00000526691.5
</td>
<td style="text-align:right;">
119.41084
</td>
<td style="text-align:right;">
2.319913
</td>
<td style="text-align:right;">
4.480485
</td>
<td style="text-align:right;">
7.40e-06
</td>
<td style="text-align:right;">
0.0095922
</td>
</tr>
</tbody>
</table>
