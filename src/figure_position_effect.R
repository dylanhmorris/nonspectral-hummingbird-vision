#!/usr/bin/env Rscript

#####################################
## name: figure_position_effect.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## visualize posteriors for position effect
##
##
####################################

suppressPackageStartupMessages(library(ggplot2))   # for general plotting
suppressPackageStartupMessages(library(magrittr))  # for pipe operator 
suppressPackageStartupMessages(library(dplyr))     # for ungroup()
suppressPackageStartupMessages(library(tidybayes)) # for spread_draws()
suppressPackageStartupMessages(library(ggridges))  # for theme_ridges()
suppressPackageStartupMessages(library(cowplot))   # for publication-ready
                                                   # graphics
suppressPackageStartupMessages(library(readr))     # for read_csv

args <- commandArgs(trailingOnly=TRUE)
mcmc_fit_path <-args[1]
data_path <- args[2]
outpath <- args[3]

dat <- read_csv(data_path,
                col_types = cols())
fit <- readRDS(mcmc_fit_path)


## order factors in dat correctly
dat$experimentClass =
    factor(dat$experimentClass,
           levels=c(
               "spectral",
               "non-spectral",
               "null"))

tidyfit <- fit %>% spread_draws(loc_true_bias[Location]) %>% ungroup()

quants = tidyfit %>% group_by(Location) %>%
    do(tibble(bias_quant=quantile(.$loc_true_bias, ppoints(100))))


loc_year = c(
    "2016",
    "2017",
    "2018")

quants$loc_year = loc_year[quants$Location]

pointcolor = "#56c9ff"

xlab_string <- "left bias
(probability that a naive bird
chooses lefthand feeder)"
ylab_string <- "posterior frequency"


position_fig <- quants %>%
    ungroup() %>%
    ggplot(aes(x=bias_quant)) +
    geom_vline(xintercept=0.5) + 
    geom_dotplot(
        fill=pointcolor,
        alpha=1,
        binwidth=0.009) +
    scale_fill_distiller() +
    scale_y_continuous(breaks = NULL) +
    facet_grid(rows=vars(loc_year)) +
    coord_cartesian(xlim=c(0.25, 0.75)) +
    xlab(xlab_string) +
    ylab(ylab_string) +
    theme_classic() +
    theme_ridges() +
    panel_border()


save_plot(outpath,
          position_fig,
          base_height=9,
          base_aspect_ratio=0.75)


