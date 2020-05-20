#!/usr/bin/env Rscript

#####################################
## name: figure_pp_check.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## plot posterior predictive checks 
## for all models
##
####################################

suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(bayesplot))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(readr))

## read command line args
args <- commandArgs(trailingOnly=TRUE)
n_models <- length(args) -2
data_path <- args[n_models+1]
outpath <- args[n_models+2]
mcmc_fit_paths <- args[1:n_models]

mcmc_fit_names <- c(
    'main model',
    'multilevel model',
    'no error model')


## get original data
dat <- read_csv(data_path,
                col_types = cols())

## remove trials that carry no info
dat <- dat[(dat$cumulativeHitsAtTrialStart +
            dat$ChoseSugar +
            dat$ChoseWater) > 0,]

## get mcmc chains

set.seed(543522) # reproducibility

n_to_plot = 50 # only plot a random subsample, for efficiency


#####################################
## styling parameters
#####################################

## define plot geometry
if(n_models > 3){
  co <- 2
  ro <- 2
  left_plots <- c(1, 4, 6)
  bottom_plots <- c(5, 6, 7)
  
} else{
    co <- 1
    ro <- n_models
    left_plots <- 1:3
    bottom_plots <- n_models
}

## set xlimits/coordinate system
pp_xlim <- c(0, 25)
pp_coords <- coord_cartesian(xlim=pp_xlim)

## define plot styling
pp_blue <- "#28c8ff"  # color for predictive draw lines

cat("\n\nPlotting posterior predictive checks...\n")
cat("Ignore coordinate system and colour warnings;",
    "these are expected behavior\n")


pp_check_plots <- list()
for(k in 1:length(mcmc_fit_paths)){
    model_name <- mcmc_fit_names[k]
    cat(sprintf("\nPlotting pp checks for %s...\n\n",
                model_name))


    fit <- readRDS(mcmc_fit_paths[k])
    chains <- extract(fit)
    print(names(chains))


    model_dat <- dat
    rg_uv_only <- grepl("regression", model_name)
    no_error <- grepl("no error", model_name)
    multilevel <- grepl("multilevel", model_name)

    if (rg_uv_only){
        cat("\nOnly considering red-uv and green-uv axis data...\n\n")
        model_dat <- model_dat[model_dat$axisID %in% c(1, 2)]
    }

    if(no_error) {
        prob_array <- chains$p_correct_mean
    } else {
        prob_array <- chains$p_correct_overall
    }
    
    n_checks = length(prob_array[,1])
    print(n_checks)
    
    sample_checks = sample(1:n_checks, n_to_plot)

    real_succs <- model_dat$ChoseSugar
    n_hits <- model_dat$ChoseSugar + model_dat$ChoseWater
    n_trials_total <- dim(model_dat)[1]
    yreps <- apply(
        prob_array[sample_checks,],
        MARGIN = 1,
        FUN = function(x){rbinom(n_trials_total,
                                 size = n_hits,
                                 p = x)})
    yreps <- t(yreps)
    
    plot_k <- pp_check(real_succs,
                       yrep = yreps,
                       fun = ppc_dens_overlay,
                       alpha = 0.2) +
        ggtitle(model_name) +
        pp_coords + 
        xlab(element_blank()) +
        ylab(element_blank()) +
        scale_color_manual(name="", 
                           labels = c("real data",
                                      "posterior predictive draws"),
                           values = c("black",
                                      pp_blue)) +
        theme_cowplot()

    if(k %in% bottom_plots)
        plot_k <- plot_k + xlab("rewarded feeder visits")
    if(k %in% left_plots)
        plot_k <- plot_k + ylab("density")

    pp_check_plots[[k]] <- plot_k

}

lab_letters <- c('a', 'b', 'c', 'd', 'e', 'f', 'g')

pp_check_plot <- plot_grid(plotlist=pp_check_plots,
                           nrow = ro,
                           ncol = co,
                           labels=lab_letters)

## save plot to outpath
cat("\nSaving plot...\n\n")
save_plot(outpath,
          pp_check_plot,
          base_width = 8,
          nrow = ro,
          ncol = co)
