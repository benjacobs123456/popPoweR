
<!-- README.md is generated from README.Rmd. Please edit that file -->

# popPoweR

<!-- badges: start -->

<!-- badges: end -->

GWAS in populations of non-European ancestry is often not powered to
detect variants at a stringent genome-wide significance level. With
smaller sample sizes, the goal is often to assess whether risk variants
identified in European-ancestry populations show evidence of replication
in other ancestries. Power varies according to allele frequency, with
rarer variants being harder to detect with a fixed sample size. Allele
frequencies vary according to ancestral population. This package
provides a quick way of integrating GWAS summary stats with allele
frequencies from 1,000 genomes and using these frequencies to estimate
the power for replicating each signal in non-European ancestries.

## Installation

You can install the package from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("benjacobs123456/popPoweR")
```

## Workflow

Here’s the basic workflow.

``` r
library(popPoweR)

## 1 - read in the European-ancestry GWAS
gwas = readr::read_tsv("../gwas_sumstats.tsv")

## 2 - format the sumstats (rename columns)
df = format_gwas(gwas = gwas,
                 chr = "chr",
                 effect_allele = "A1",
                 other_allele = "A2",
                 snp = "rsid",
                 eaf = "A1freq",
                 beta = "BETA",
                 n_case = "cases",
                 n_control = "controls")

## 3 - combine gwas with 1kg frequencies
kg_combo = combine_gwas_with_1kg_freqs(df)

## 4 - and work out power across every SNP
outputdf = power_calc_per_pop(kg_combo,n_perm = 1000,alpha=0.05,pop="AFR")
```
