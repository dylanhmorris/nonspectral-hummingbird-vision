#!/usr/bin/env Rscript

#####################################
## name: figure_experience.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## visualize posteriors for
## becoming experienced and related
## quantities (visitation rates)
##
##
####################################

script_packages <- c(
    'readr',     # for read_csv()
    'ggplot2',   # for plotting
    'cowplot',   # for publication-ready ggplot
    'magrittr',  # for pipe operator %>%
    'tidybayes', # for spread_draws(),
    'dplyr',     # for sql-style manipulations
    'tidyr')     # for cartesian product via crossing(),

## silently load packages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}


#####################
## define functions
#####################

experience_func <- function(equilibrium_E,
                            experience_saturation_rate,
                            time){
    return(-equilibrium_E * exp(-experience_saturation_rate * time) +
           equilibrium_E)
}

################################
## process script call inputs
## and read in data
################################


args <- commandArgs(trailingOnly=TRUE)
mcmc_fit_path <-args[1]
data_path <- args[2]
outpath <- args[3]

set.seed(501325) # reproducible!

dat <- read_csv(data_path,
                col_types = cols())
fit <- readRDS(mcmc_fit_path)

##############################
## Set aesthetic/styling parameters
##############################
darkblue <- "#00356B"

## various colors, sizes, transparencies
posterior_pointfill <- "#56c9ff"
hitrate_pointfill   <- "#56c9ff"

pointborder        <- 'black'
pointsize          <- 2
pointstroke        <- 0.5
linecolor          <- darkblue
func_form_alpha    <- 0.1
func_form_linesize <- 1.25
hitrate_alpha      <- 0.5 

## number of lines for plotting functional
## form posterior 
n_lines <- 100

## jitter for hitrate plots
jitwid <- 0.01
jith   <- 0.0

## overall plot styling
plot_width_cm <- 18.3
plot_width_in <- plot_width_cm * 0.393701
theme_set(
    theme_classic() +
    background_grid(major = "xy",
                    minor = "none",
                    size.major = 0.5))




##############################
## plot posterior dotplots
##############################

eq_E_draws <- fit %>%
    spread_draws(equilibrium_E) %>%
    ungroup()

e_quants = eq_E_draws %>%
    do(tibble(E_quant=quantile(.$equilibrium_E, ppoints(100))))


n_bins <- 25
e_max  <- 0.6
e_min  <- 0.25
e_range <- e_max - e_min

e_fig <- e_quants %>%
    ungroup() %>%
    ggplot(aes(x=E_quant)) +
    coord_cartesian(xlim=c(e_min, e_max)) + 
    geom_dotplot(
        fill=posterior_pointfill,
        alpha=1,
        binwidth=e_range/n_bins) +
    scale_fill_distiller() +
    scale_y_continuous(breaks = NULL) +
    xlab("equilibrium fraction of visits\nfrom experienced birds (E)") +
    ylab("posterior frequency")

half_sat_draws <- fit %>%
    spread_draws(experience_half_saturation_time) %>% ungroup()

half_sat_quants = half_sat_draws %>%
    do(tibble(half_sat_quant=quantile(.$experience_half_saturation_time,
                                      ppoints(100))))

h_max <- 60
h_min <- 0
h_range <- h_max - h_min

h_fig <- half_sat_quants %>%
    ungroup() %>%
    ggplot(aes(x=half_sat_quant)) +
    coord_cartesian(xlim=c(h_min, h_max)) + 
    geom_dotplot(
        fill=posterior_pointfill,
        alpha=1,
        binwidth=h_range/n_bins) +
    scale_fill_distiller() +
    scale_y_continuous(breaks = NULL) +
    xlab("experience level half-saturation wait\n(in feeder visits)") +
    ylab(element_blank())

#############################################
## plot distribution of experience functions
## (as superimposed semi-transparent lines)
#############################################

exp_xs <- tibble(cumul_hits = seq(0, 250, 0.5))

func_samples <- fit %>%
    spread_draws(equilibrium_E, experience_saturation_rate) %>%
    sample_n(n_lines) # get 100 random posterior lines

to_plot <- func_samples %>% crossing(exp_xs)


#############################################
## plot empirical distribution of hit rates
## as a function of euclid and correctness
#############################################

exp_funcs_fig <- to_plot %>%
    ggplot(aes(x = cumul_hits,
               y = experience_func(equilibrium_E,
                                   experience_saturation_rate,
                                   cumul_hits),
               group = .draw)) +
    geom_line(alpha=func_form_alpha,
              color=linecolor,
              size=func_form_linesize) +
    xlab("cumulative feeder visits") +
    ylab("fraction of experienced birds")

dat$hit_rate <- (dat$ChoseSugar + dat$ChoseWater) / dat$trialDuration
dat$correctness <- dat$ChoseSugar/ (dat$ChoseSugar + dat$ChoseWater)

hits_by_euclid <- ggplot(data=dat) +
    geom_point(aes(x=Euclid, y=hit_rate),
               shape = 21,
               colour=pointborder,
               fill=hitrate_pointfill,
               alpha=hitrate_alpha,
               size=pointsize,
               stroke=pointstroke,
               position=position_jitter(width=jitwid,
                                        height=jith,
                                        seed = 5)) +
    labs(x = "Euclidean color distance",
         y = "feeder visits per minute") 

##  exclude trials with <5 hits.
datsub<-subset(dat, dat$nHits>5)

hits_by_correctness <- ggplot(data=datsub) +
    geom_point(aes(x=correctness, y=hit_rate),
               shape = 21,
               colour=pointborder,
               fill=hitrate_pointfill,
               alpha=hitrate_alpha,
               size=pointsize,
               stroke=pointstroke,
               position=position_jitter(width=jitwid,
                                        height=jith,
                                        seed = 5)) +
    labs(x = "percent of visits to reward",
         y = element_blank()) +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank())



########################################
# assemble and save full paneled figure
########################################

top_row <- plot_grid(
    e_fig,
    h_fig,
    labels=c('a', 'b'),
    align='h')

bottom_row <- plot_grid(
    hits_by_euclid,
    hits_by_correctness,
    labels=c('d', 'e'),
    align='h')


full_fig <- plot_grid(
    top_row,
    exp_funcs_fig,
    bottom_row,
    labels = c('', 'c', ''),
    ncol=1)

save_plot(outpath,
          full_fig,
          base_width=plot_width_in,
          base_height=plot_width_in * 1.5)

