#!/usr/bin/env Rscript

#####################################
## name: extract_reported_quantities.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## read in mcmc chains
## and output quantities that
## we report in the paper
##
####################################

suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(bayesplot)) # for rhat, neff_ratio

## read command line args
args <- commandArgs(trailingOnly=TRUE)
n_models <- length(args) - 2
mcmc_fit_paths <- args[1:n_models]
data_path <- args[n_models + 1]
outpath <- args[n_models + 2]

## read data
dat <- read.csv(data_path)

## useful functions
get_experiments_by_class <- function(class_name,
                                     dat=dat){
    return(unique(
        dat$ExperimentNumber[dat$experimentClass == class_name]))
}

get_correctness_stat <- function(chains,
                                 experimentClass,
                                 dat=dat,
                                 stat=mean,
                                 margin=1){
    experiments <- get_experiments_by_class(experimentClass,
                                            dat=dat)
    return(
        apply(chains$unbiased_p_correct_experienced[, experiments],
              margin,
              stat))
}

list_of_quantities <- list()

## do main model chain quantity extraction
cat("\n Reading in main model...\n")
fit <- readRDS(mcmc_fit_paths[1])
chains <- extract(fit)
cat(paste0("Model read successfully, ",
           "extracting reported ",
           "quantities...\n"))

null_means <- get_correctness_stat(
    chains,
    'null',
    dat=dat,
    stat=mean)

spec_means <- get_correctness_stat(
  chains,
  'spectral',
  dat=dat,
  stat=mean)

nonspec_means <- get_correctness_stat(
    chains,
    'nonspectral',
    dat=dat,
    stat=mean)

list_of_quantities['main_null_mean_mean_correctness'] <-
    mean(null_means)

list_of_quantities['main_spec_mean_mean_correctness'] <-
  mean(spec_means)

list_of_quantities['main_nonspec_mean_mean_correctness'] <-
    mean(nonspec_means)

list_of_quantities['main_nonspec_q01_mean_correctness'] <-
    quantile(nonspec_means, 0.01)

list_of_quantities['main_nonspec_q05_mean_correctness'] <-
    quantile(nonspec_means, 0.05)

list_of_quantities['main_spec_q01_mean_correctness'] <-
  quantile(spec_means, 0.01)

list_of_quantities['main_spec_q05_mean_correctness'] <-
  quantile(spec_means, 0.05)

list_of_quantities['main_null_q95_mean_correctness'] <-
    quantile(null_means, 0.95)

list_of_quantities['main_null_q99_mean_correctness'] <-
    quantile(null_means, 0.99)

list_of_quantities['main_nonspec_min_q50_correctness'] <-
    min(get_correctness_stat(
        chains,
        'nonspectral',
        dat=dat,
        margin=2,
        stat=function(x){quantile(x, 0.5)}))

list_of_quantities['main_nonspec_min_mean_correctness'] <-
    min(get_correctness_stat(
        chains,
        'nonspectral',
        dat=dat,
        margin=2,
        stat=mean))

list_of_quantities['main_spec_min_q50_correctness'] <-
  min(get_correctness_stat(
    chains,
    'spectral',
    dat=dat,
    margin=2,
    stat=function(x){quantile(x, 0.5)}))

list_of_quantities['main_spec_min_mean_correctness'] <-
  min(get_correctness_stat(
    chains,
    'spectral',
    dat=dat,
    margin=2,
    stat=mean))

list_of_quantities['main_null_max_q50_correctness'] <-
    max(get_correctness_stat(
        chains,
        'null',
        dat=dat,
        margin=2,
        stat=function(x){quantile(x, 0.5)}))

list_of_quantities['main_mean_sigma_error'] <-
    mean(chains$sd_error)

list_of_quantities['main_q95_sigma_error'] <-
    quantile(chains$sd_error, 0.95)


## key-value to dataframe
df <- stack(list_of_quantities)
names(df) <- c("value", "quantity")

## reorder for better display
df <- df[c("quantity", "value")]
## save to outpath

## output to file
cat(sprintf("saving to %s\n",
            outpath))
write.csv(df, outpath)
