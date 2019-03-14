library(glue)
library(readxl)
library(tidyverse)
library(readr)
library(here)

plate_indexes = read_excel(here('salmon_indexes/Copy of DGE well barcodes_384_set1_newp11.xlsx')) %>%
  select(-plate_id)

write_csv(plate_indexes,'plate_ids.csv')