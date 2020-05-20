#!/usr/bin/env Rscript

#####################################
## name: fit_stan_model.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## fit the specified stan model
## to the specified clean dataset,
## and save the result as an .Rds
## file
##
####################################

suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(parallel))


## read command line args
args <- commandArgs(trailingOnly=TRUE)
model_src_path <-args[1]
datapath <- args[2]
outpath <- args[3]


## set stan options
n_cores <- parallel::detectCores()
options(mc.cores = n_cores)
rstan_options(auto_write = TRUE)
niter <- 2000
thin <- 1
nchains = n_cores
adapt_d= 0.85
max_tree = 15
fixed_seed = 21

## load data
dat <- read.csv(datapath)
dat$dtStart <- as.POSIXct(dat$dtStart)
dat$dtEnd <- as.POSIXct(dat$dtEnd)

    
#############################
## set hyperparameters

main_hyperprior_list <- list(
    min_learning_effect = 0, # birds not worse than random
    hyper_mean_sd_left_bias = 0,
    hyper_sd_sd_left_bias = 0.5, # need to rule out error sds > 1.5 --> 3 sigma
    hyper_mean_sd_error = 0,
    hyper_sd_sd_error = 0.5, # need to rule out error sds > 1.5 --> 3 sigma
    hyper_alpha_equilibrium_E = 1, # flat beta dist
    hyper_beta_equilibrium_E = 1, 
    hyper_mean_experience_half_saturation_time = 100,
    hyper_sd_experience_half_saturation_time = 200,
    hyper_mode_p_correct_experienced = 0.5,
    hyper_certainty_p_correct_experienced = 2, # flat beta
    hyper_mode_loc_left_bias = 0.5,
    hyper_certainty_loc_left_bias = 6) # weakly informative beta, mode 0.5

multilevel_hyperprior_list <- list(
    hyper_mode_treatment_mode = 0.5,
    hyper_certainty_treatment_mode = 2,
    hyper_mean_treatment_certainty = 5) # don't allow too much pooling

#############################
## format data for model
############################

## handle regression model on uv/green and uv/red from
## same script as others
rg_uv_only <- grepl("regression", model_src_path)

if (rg_uv_only){
    cat("\nOnly considering red-uv and green-uv axis data...\n")
    dat <- dat[dat$axis %in% c(1, 2),]

    ## renumber experiments...
    dat$ExperimentNumber <- with(dat,
                                 match(ExperimentNumber,
                                       unique(ExperimentNumber)))
}


## get number of experiments
n_experiments <- max(dat$ExperimentNumber)

# enumerate locations, encoding year
dat$Location <- as.factor(format(dat$dtStart, "%Y")) 
# change location levels to 1...n_locations
levels(dat$Location) <- seq(1, nlevels(dat$Location))
# convert to number 
dat$Location <- as.numeric(dat$Location)

## make lookup table for experiment locations
locations <- unique(dat[c('ExperimentNumber', 'Location')])
## check that 1 location per experiment
stopifnot(length(locations$Location) == n_experiments)
locations <- locations$Location
n_locations <- max(locations)

# enumerate treatment IDs from experiment classes
dat$treatmentID <- as.numeric(
  factor(dat$experimentClass, 
         levels=unique(dat$experimentClass)))

treatments <- unique(dat[c('ExperimentNumber', 'treatmentID')])
stopifnot(length(treatments$treatmentID) == n_experiments)
treatment_id = treatments$treatmentID
n_treatments = max(treatment_id)

euclids <- unique(dat[c('ExperimentNumber', 'Euclid')])
stopifnot(length(euclids$Euclid) == n_experiments)
euclid <-  euclids$Euclid

## note which axis (red/uv, green/uv or neither)
## an observation row lies along
axes <- unique(dat[c('ExperimentNumber', 'axisID')])
axis_id = axes$axisID
axis_id[is.na(axis_id)] <- 4
n_axes = max(axis_id)

## calculate total visits / trials

## remove trials that provide no information
dat <- dat[(dat$cumulativeHitsAtTrialStart +
            dat$ChoseSugar +
            dat$ChoseWater) > 0,]

n_trials_total <- dim(dat)[1]
n_trials_by_hit <- sum(dat$stoppingCondition == 0)
n_trials_by_time <- sum(dat$stoppingCondition == 1)
stopifnot(n_trials_total == n_trials_by_hit + n_trials_by_time)


## make data into list
data_list <- list(
    n_experiments = n_experiments,
    n_trials_total = n_trials_total,
    n_locations = n_locations,
    n_axes = n_axes,
    n_treatments = n_treatments,
    trial_duration = dat$trialDuration,
    cumulative_duration_at_trial_start = dat$cumulativeDurationAtTrialStart,
    cumulative_hits_at_trial_start = dat$cumulativeHitsAtTrialStart,
    n_successes = dat$ChoseSugar,
    n_failures = dat$ChoseWater,
    correct_side = dat$CorrectSide,
    experiment_id = dat$ExperimentNumber,
    euclid = euclid,
    location = locations,
    treatment_id = treatment_id,
    axis_id = axis_id,
    stopping_condition = dat$stoppingCondition,
    n_trials_by_hit = n_trials_by_hit)
    

# concatenate with prior choices
stan_data = c(
    main_hyperprior_list,
    multilevel_hyperprior_list,
    data_list)

###############################
## Compile, fit, and save model
###############################
fit <- stan(
    model_src_path,
    data = stan_data,
    iter = niter,
    seed = fixed_seed,
    chains = nchains,
    thin = thin,
    control = list(max_treedepth = max_tree,
                   adapt_delta = adapt_d))

cat(sprintf("saving posterior samples to file %s...\n",
            outpath))

saveRDS(fit, outpath)

warnings()
