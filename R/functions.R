#' Calculate power for GWAS of a binary trait
#'
#' @param eaf effect allele frequency
#' @param beta_alt beta under the alternate hypothesis (i.e. per-allele log odds ratio in your discovery GWAS)
#' @param n_case number of cases
#' @param n_control number of controls
#' @param alpha desired alpha
#' @param n_perm number of permutations for power calculations
#'
#' @return the estimated power for each SNP given the allele frequency, n, alpha, effect size, and number of permutations. Note that this is for a one-tailed test.
#' @export

gwas_power = function(eaf = 0.5,
                      beta_alt = 0.5,
                      n_case = 1000,
                      n_control = 1000,
                      alpha = 5e-8,
                      n_perm = 1000){

  # guess se of SNP effects and assume it's constant across SNPs
  test_genos = rbinom(n = n_control + n_case, size = 2, prob = eaf)
  test_phenos = rbinom(n = n_control + n_case, size = 1, prob = n_case / (n_case + n_control))

  # make model
  model = glm(test_phenos ~ test_genos,family=binomial(link="logit"))
  se = summary(model)$coefficients[2,2]

  # simulate SNP effects from a normal using this SE
  betas = rnorm(n = n_perm, mean = beta_alt, sd = se)
  z_scores = betas/se


  pval = if(beta_alt>0){
    (1-pnorm(z_scores))/2
  } else {
    pnorm(z_scores)/2
  }

  power = sum(pval<alpha) / length(pval)
  return(power)
}


#' Convert GWAS summary statistics into required format
#'
#' @param gwas GWAS summary statistics
#' @param chr Column name for chromosome - a string
#' @param snp Column name for rsid of SNP - a string
#' @param effect_allele Column name for effect allele (i.e. the allele for the signed statistic beta) - a string
#' @param other_allele Column name for the non-effect allele - a string
#' @param eaf Column name for effect allele frequency - a string
#' @param beta Column name for per-allele log odds ratio - a string
#' @param n_case Column name for number of cases - a string
#' @param n_control Column name for number of controls - a string
#'
#' @return a data frame with formatted GWAS summary statistics
#' @export
#'
#' @examples
#' df = format_gwas(gwas = gwas_sumstats,chr = "CHR",effect_allele = "A1",other_allele = "A2",snp = "SNP",eaf = "eaf",beta = "beta",n_case = "n_case",n_control = "n_control")
#' df = format_gwas(gwas = gwas_sumstats)


format_gwas = function(gwas,
                       chr = "CHR",
                       snp = "SNP",
                       effect_allele = "A1",
                       other_allele = "A2",
                       eaf = "eaf",
                       beta = "beta",
                       n_case = "N_case",
                       n_control = "N_cont")
{

  if(!is.character(c(chr,snp,effect_allele,other_allele,eaf,beta,n_case,n_control))){
    stop("Please provide column names in your GWAS")
  }

  if(!is.data.frame(gwas)){
    stop("GWAS must be a data frame")
  }

  output_gwas = data.frame(
    CHR = gwas[[chr]],
    SNP = gwas[[snp]],
    effect_allele = gwas[[effect_allele]],
    other_allele = gwas[[other_allele]],
    eaf = gwas[[eaf]],
    beta = gwas[[beta]],
    n_case = gwas[[n_case]],
    n_control = gwas[[n_control]]
  )
  return(output_gwas)
}




#' Combine GWAS with 1,000 genomes frequencies
#'
#' @param gwas Formatted GWAS summary statistics with 1kg frequencies (output of format_gwas)
#'
#' @return GWAS summary statistics with 1kg allele frequencies (and aligned alleles)
#' @export
#'
#' @examples
#' df = format_gwas(gwas = gwas_sumstats)
#' kg_combo = combine_gwas_with_1kg_freqs(df)

combine_gwas_with_1kg_freqs = function(gwas){
  # load 1kg freqs
  load("R/sysdata.rda")

  # filter input gwas to just 1kg snps
  message("Filtering ",nrow(gwas)," SNPs in input GWAS")
  gwas = gwas %>% dplyr::filter(SNP %in% kg_freqs$SNP)
  message(nrow(gwas)," SNPs from input GWAS present in 1kg")

  # and vice-versa
  kg_freqs = kg_freqs %>% dplyr::filter(SNP %in% gwas$SNP)

  # now join by rsid
  gwas = gwas %>% dplyr::left_join(kg_freqs,by=c("SNP","CHR"))

  # check alleles are aligned
  matching_alleles = gwas %>% dplyr::filter(effect_allele == A1 && other_allele == A2)
  message("There are ",nrow(matching_alleles)," matching SNPs")

  flipped_alleles = gwas %>% dplyr::filter(effect_allele == A2 && other_allele == A1)
  message("There are ",nrow(flipped_alleles)," SNPs with flipped alleles")
  flipped_alleles = flipped_alleles %>%
    dplyr::mutate(beta = beta * -1) %>%
    dplyr::mutate(eaf = 1 - eaf) %>%
    dplyr::mutate(effect_allele_1 = other_allele) %>%
    dplyr::mutate(other_allele_1 = effect_allele) %>%
    dplyr::select(-effect_allele,-other_allele) %>%
    dplyr::rename(other_allele = other_allele_1) %>%
    dplyr::rename(effect_allele = effect_allele_1)

  # print number of non-matching alleles
  binned_snps = gwas %>% dplyr::filter(!(SNP %in% matching_alleles$SNP || SNP %in% flipped_alleles$SNP))
  message("Excluded ",nrow(binned_snps)," SNPs with non-matching alleles")

  # combine
  gwas = dplyr::bind_rows(flipped_alleles,matching_alleles) %>% select(-A1,-A2)
  message(nrow(gwas)," SNPs remaining")
  return(gwas)
}


#' Power calculation
#'
#' @param gwas GWAS summary statistics with 1kg allele frequencies (output of combine_gwas_with_1kg_freqs)
#' @param n_perm Number of permutations, defaults to 100.
#' @param alpha Alpha, defaults to 0.05 (replication)
#' @param pop A string which must be one of "AFR","SAS","AMR","EAS",or "EUR" specifying the 1kg superpopulation of interest.
#'
#' @return a data frame with your input summary statistics plus a new column listing the power at each SNP in the superpopulation of interest at the desired alpha.
#' @export
#'
#' @examples
#'
#' df = format_gwas(gwas = gwas_sumstats)
#' kg_combo = combine_gwas_with_1kg_freqs(df)
#' power_stats = power_calc_per_pop(kg_combo,n_perm = 1000,alpha=0.05,pop="AFR") # for AFR
#' power_stats = power_calc_per_pop(kg_combo,n_perm = 1000,alpha=0.05,pop="SAS") # for SAS

power_calc_per_pop = function(gwas,n_perm,alpha,pop){
  if(!is.data.frame(gwas)){
    stop("GWAS input must be a data frame with 1kg allele frequencies (i.e. output of combine_gwas_with_1kg_freqs)")
  }

  if(n_perm != round(n_perm) | n_perm == 0 | n_perm > 100000){
    stop("Number of permutations must be an integer between 0 and 100000")
  }

  if(!is.numeric(alpha) | alpha == 0 | alpha > 0.5){
    stop("Alpha must be a number between 0 and 0.5")
  }

  if(!pop %in% c("AFR","SAS","EAS","EUR","AMR")){
    stop("Population must be a character string corresponding to a 1kg superpopulation, i.e. one of one of AFR, SAS, AMR, EAS, or EUR")
  }

  power_estimates = list()

  # remove SNPs where AF is missing
  gwas = gwas[!is.na(gwas[[pop]]),]

  # iterate through GWAS and calculate power for each SNP
  for(i in 1:nrow(gwas)){
    message("Estimating power at SNP ",i," of ",nrow(gwas))
    power = gwas_power(
      eaf = gwas[[pop]][i],
      beta_alt = gwas$beta[i],
      n_case = gwas$n_case[i],
      n_control = gwas$n_control[i],
      alpha = alpha,
      n_perm = n_perm
    )
    power_estimates[[i]] = power
  }

  gwas[[paste0(pop,"_alpha_",alpha,"_power")]] = unlist(power_estimates)
  return(gwas)
}
