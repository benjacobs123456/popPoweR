% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{combine_gwas_with_1kg_freqs}
\alias{combine_gwas_with_1kg_freqs}
\title{Combine GWAS with 1,000 genomes frequencies}
\usage{
combine_gwas_with_1kg_freqs(gwas)
}
\arguments{
\item{gwas}{Formatted GWAS summary statistics with 1kg frequencies (output of format_gwas)}
}
\value{
GWAS summary statistics with 1kg allele frequencies (and aligned alleles)
}
\description{
Combine GWAS with 1,000 genomes frequencies
}
\examples{
df = format_gwas(gwas = gwas_sumstats)
kg_combo = combine_gwas_with_1kg_freqs(df)
}
