# Get all sequences that do not correspond to a single mutant
# Nathan Lubock

library(magrittr)
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

# load everything (~2 mins)
wt.key <-
    read_delim(args[1], delim = ' ')

counts.per.idx <-
    read_delim(args[2], delim = ' ', col_names = c('Index', 'Barcode', 'Count'))

single.mutants <-
    read_delim(args[3], delim = ' ', col_names = c('Barcode', 'Count', 'Variant'))

var.map <-
    read_csv(args[4], col_names = c('Barcode', 'Seq', 'Count'))

# filter out barcodes that are already mapped to single mutants
# (speeds up search ~2-3 mins)
counts.per.idx %>%
    select(Barcode) %>%
    distinct() %>%
    anti_join(wt.key, by = 'Barcode') %>%
    anti_join(single.mutants, by = 'Barcode') %>%
    inner_join(var.map, by = 'Barcode') %>%
    select(Barcode, Seq) %>%
    format_tsv(col_names = FALSE) %>%
    cat()
