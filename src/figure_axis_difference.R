#!/usr/bin/env Rscript

#####################################
## name: figure_axis_difference.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## plot change in estimated correctness
## as a function of euclidean distance
## for uvred and uvgreen axes
## 
####################################

#######################
# load packages
#######################

script_packages <- c(
    'ggplot2',    # for plotting
    'magrittr',   # for (ceci n'est pas une) pipe operator %>%
    'cowplot',    # for publication-ready ggplot
    'ggridges',   # for theme_ridges()
    'tidybayes',  # for spread_draws()
    'readr',      # for read_csv()
    'dplyr',      # for sample_n()
    'tidyr',      # for cartesian product via crossing()
    'tibble')     # for add_column()

for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}


############################
## styling/aesthetic parameters
############################

plot_width_cm <- 18.3
plot_width_in <- plot_width_cm * 0.393701
fineness <- 100
n_lines  <- 100

violin_width <- 0.075
violin_alpha <- 0.25

palette <- list()
palette["UV"] <- "#EFE8F3" 
palette["Red"] <- "#ED2024"
palette["Blue"] <- "#3953A4"
palette["Green"] <- "#69BD45"
palette["Yellow"] <- "#F6EB14"
palette["Purple"] <- "#B9529F"
palette["White"] <- "#FFFFFF"

############################
## get args, read data
############################

args <- commandArgs(trailingOnly=TRUE)
mcmc_fit_path <-args[1]
data_path <- args[2]
outpath <- args[3]

dat <- read_csv(data_path,
                col_types = cols())
fit <- readRDS(mcmc_fit_path)


############################
## data prep
############################

## subset data to just look at uvred and uvgreen axes
dat <- dat[dat$axisID %in% c(1, 2, 3), ]

dat <-  rbind(
    transform(subset(dat, axisID %in% c(1, 3)), axisID=1),
    transform(subset(dat, axisID %in% c(2, 3)), axisID=2))

euclid_table <- unique(dat[c('ExperimentNumber',
                             'ExperimentDisplayName',
                             'Euclid',
                             'axisID')])

euclid_table <- euclid_table[order(euclid_table$ExperimentNumber), ]

euclid_table$axis_name <- ifelse(
    euclid_table$axisID == 1,
    "UV+Green axis",
    "UV+Red axis")

############################
## plotting 
############################
    
min_quant <- 0.025
max_quant <- 0.975

learning_effects <- fit %>%
    spread_draws(unbiased_p_correct_experienced[exp_id]) %>%
    group_by(exp_id) %>%
    filter(unbiased_p_correct_experienced > 
           quantile(unbiased_p_correct_experienced,
                    min_quant)) %>%
    filter(unbiased_p_correct_experienced <
           quantile(unbiased_p_correct_experienced,
                    max_quant)) %>%
    inner_join(euclid_table,
               by=c('exp_id' = 'ExperimentNumber'))


fig <- learning_effects %>%
    ggplot(aes(x=Euclid,
               y=unbiased_p_correct_experienced,
               group=exp_id,
               fill=axis_name)) + 
    geom_violin(position="identity",
                color="black",
                draw_quantiles=c(0.5),
                alpha=violin_alpha,
                width=violin_width) +
    scale_fill_manual(values=c(palette[['Green']],
                               palette[['Red']])) +
    facet_wrap(vars(axis_name), ncol=1) +
    coord_cartesian(xlim=c(0, 1),
                    ylim=c(0.5, 1)) +
    xlab("Euclidean distance in avian colorspace") + 
    ylab("Probability of choosing rewarded color") +
    theme_classic(base_size=22) +
    background_grid(major = "xy",
                    minor = "none",
                    size.major = 0.5) +
    theme(legend.position = "none") + 
    panel_border()


save_plot(outpath,
          fig,
          base_width=plot_width_in,
          base_height=plot_width_in * 2)

 
