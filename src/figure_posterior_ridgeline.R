#!/usr/bin/env Rscript

#####################################
## name: figure_posterior_ridgeline.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## visualize posteriors of interest
##
##
####################################

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(tidybayes))
suppressPackageStartupMessages(library(forcats))  ## for fct_rev()
suppressPackageStartupMessages(library(RColorBrewer))
## for scale_fill_distiller
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(cowplot))
# for metadata generation
suppressPackageStartupMessages(library(stringr))

###############################
## styling pre-definitions / setup
###############################
palette <- list()
palette["UV"] <- "#EFE8F3"
palette["Red"] <- "#ED2024"
palette["Blue"] <- "#3953A4"
palette["Green"] <- "#69BD45"
palette["Yellow"] <- "#F6EB14"
palette["Purple"] <- "#B9529F"
palette["White"] <- "#FFFFFF"

den_scale <- 1 ## den_scale <- 1 makes ridgeline densities not overlap
rel_min_height <- 0.01 ## minimum relative density to plot

x_light <- 0.32 ## where to put the indicators of light color
light_pad <- 0.07 ## how much to pad them
light_wid <- 0.05 ## how wide they should be

lower_ridgecolor <- "black"
upper_ridgecolor <- "#909aaf"

###############################
## plotting functions
###############################

label_fmt <- function(value){
    pattern = "\\^([[:digit:]]{2})"
    replacement = "_{\\1}"
    val <- gsub(pattern, replacement, value)
    return(TeX(val))
}

ggGradientCircle <- function(color1,
                             color2=NULL,
                             fraction1=1,
                             fineness=500,
                             blend_region=0.2,
                             bordersize=1,
                             midcolor=rgb(240/255,233/255,255/255)){
    if (is.null(color2) || color2 == color1) {
        color2 <- color1
        blend_region <- 0
    }

    heatmap_dat <- matrix(0, nrow=fineness, ncol=fineness)

    border <- floor(fraction1 * fineness)
    heatmap_dat[0:border,] <- 1

    if (blend_region > 0){
        n_blend <- floor(blend_region * fineness / 2)
        blend_low <- max(border - n_blend, 1)
        blend_high <- min(border + n_blend, fineness - 1)

        linear_blend <- matrix(
            rep(seq(1, 0, length.out=blend_high - blend_low + 1),
                fineness),
            nrow=blend_high - blend_low + 1,
            ncol=fineness,
            byrow=FALSE)

        heatmap_dat[blend_low:blend_high, ] <- linear_blend
    } else { midcolor <- color1 }
    
    heatmap_tidy <- melt(heatmap_dat,
                         varnames=c('x','y'))

    ## get rid of anything outside circle
    rad <- fineness / 2 
    outside <- (heatmap_tidy$x - rad)^2 + (heatmap_tidy$y - rad)^2 > rad^2

    heatmap_tidy$value[outside] <- NA

    ## make path for circular border
    t <- (0:6300)/1000
    circle <- data.frame(x=rad * cos(t) + rad,
                         y=rad * sin(t) + rad, t)


    plot <- heatmap_tidy %>%
        ggplot(aes(x=x,
                   y=y,
                   fill=value)) +
        geom_raster(interpolate = TRUE) +
        scale_fill_gradientn(colors=c(color1,
                                      midcolor,
                                      color2),
                             na.value='transparent') +
        annotate(
            'path',
            x=circle$x,
            y=circle$y,
            colour='black',
            size=bordersize) +
        coord_fixed() + 
        theme_nothing()

    return(plot)
}

# function to generate identity and proportion of colors from experiment display name
GenMetaData <- function(names) {
  
outmeta <-  names %>%
    str_remove(., "Null: ") %>% # get names, tidy and split into two parts
    str_remove(., ' [0-9]{1}') %>%
    str_split(., ' v ') %>%
    as.data.frame() %>%
    t() %>% # bodge!
    as.data.frame() %>%
    rename(color1 = V1, color2 = V2) %>%
    mutate(isfraccolor1 = ifelse(str_detect(color1, "\\^"), 1, 0)) %>% # detect if colors are simple of mixtures
    mutate(isfraccolor2 = ifelse(str_detect(color2, "\\^"), 1, 0)) %>%
    mutate(color1A = ifelse(isfraccolor1==1, str_extract(color1, "[A-z]*(?=\\^)"), as.character(color1))) %>% # get individual colors for first color
    mutate(color1B = ifelse(isfraccolor1==1, str_extract(color1, "(?<=\\+).*(?=\\^)"), as.character(color1))) %>%
    mutate(color2A = ifelse(isfraccolor2==1, str_extract(color2, "[A-z]*(?=\\^)"), as.character(color2))) %>% # get individual colors for second color
    mutate(color2B = ifelse(isfraccolor2==1, str_extract(color2, "(?<=\\+).*(?=\\^)"), as.character(color2))) %>%
    mutate(frac1A = ifelse(isfraccolor1 == 1, as.numeric(str_extract(color1, "[0-9].(?=\\+)")), 0)) %>% # get actual proportions for color 1
    mutate(frac1B = ifelse(isfraccolor1 == 1, as.numeric(str_extract(color1, "[0-9].$")), 0)) %>% 
    mutate(frac2A = ifelse(isfraccolor2 == 1, as.numeric(str_extract(color2, "[0-9].(?=\\+)")), 0)) %>% # get actual proportions for color 2
    mutate(frac2B = ifelse(isfraccolor2 == 1, as.numeric(str_extract(color2, "[0-9].$")), 0)) %>%
    mutate(frac2B = ifelse(color2 == 'UV^22+Yellow', 100-22, frac2B)) %>% # fix differently formatted yellow experiment
    mutate(color2B = ifelse(color2 == 'UV^22+Yellow', 'Yellow', color2B)) %>%
    mutate(sum1 = frac1A + frac1B) %>%
    mutate(sum2 =  frac2A + frac2B) %>%
    mutate(frac1A = ifelse(is.nan(frac1A/sum1), 0, frac1A/sum1)) %>% # get fractions from proportions, where nan (steady colors) have zeros
    mutate(frac1B = ifelse(is.nan(frac1B/sum1), 0, frac1B/sum1)) %>%
    mutate(frac2A = ifelse(is.nan(frac2A/sum2), 0, frac2A/sum2)) %>%
    mutate(frac2B = ifelse(is.nan(frac2B/sum2), 0, frac2B/sum2)) %>%
    mutate(exp_unique_name = names) # recover the original experiment names too

return(outmeta)

}


########################################
## data read-in and cleaning
########################################

args <- commandArgs(trailingOnly=TRUE)
mcmc_fit_path <-args[1]
data_path <- args[2]
outpath <- args[3]

cat('\n\nReading data...\n')
dat <- read.csv(data_path, stringsAsFactors=FALSE)
metadata <- GenMetaData(unique(dat$ExperimentDisplayName))
fit <- readRDS(mcmc_fit_path)

## handle regression model on uv/green and uv/red from
## same script as others
rg_uv_only <- grepl("regression", mcmc_fit_path)

if (rg_uv_only){
    cat("\nOnly considering red-uv and green-uv axis data...\n")
    dat <- dat[dat$axis %in% c(1, 2),]

    ## renumber experiments...
    dat$ExperimentNumber <- with(dat,
                                 match(ExperimentNumber,
                                       unique(ExperimentNumber)))
}


## order factors in dat correctly
dat$experimentClass <-
  factor(dat$experimentClass,
         levels=c(
           "spectral",
           "nonspectral",
           "null"))

## get one row per experiment
exp_dat <- dat %>% group_by(ExperimentNumber) %>%
    slice(1) %>% ungroup() %>%
    select(
        # Experiment=Experiment,
        exp_unique_name=ExperimentDisplayName,
        exp_id=ExperimentNumber,
        euclid=Euclid,
        experimentClass=experimentClass)

exp_dat <- exp_dat %>% inner_join(metadata,
                                  by="exp_unique_name")                                  

order_by_euclid <- TRUE 

if(order_by_euclid){

    exp_name_levels <- exp_dat$exp_unique_name[order(exp_dat$euclid,
                                                     decreasing=T)]
    exp_dat$exp_unique_name <-
        factor(exp_dat$exp_unique_name,
               levels=exp_name_levels,
               ordered=T)
} else { # else use their order in the data
  
    exp_dat$exp_unique_name <- factor(exp_dat$exp_unique_name,
                                      levels=exp_dat$exp_unique_name,
                                      ordered=T)
}

exp_dat <- exp_dat[order(exp_dat$exp_unique_name,
                         decreasing=F),]



tidyfit <- fit %>%
    spread_draws(
        unbiased_p_correct_experienced[exp_id]) %>%
    ungroup()

tidyfit <- tidyfit %>%
    left_join(exp_dat, by="exp_id")

cat("plotting ridgelines...\n")

ridgeplot <- tidyfit %>%
    ungroup() %>%
    ggplot(aes(x=unbiased_p_correct_experienced,
               y=fct_rev(exp_unique_name),
               colour=euclid,
               fill=euclid)) +
    geom_density_ridges(
        rel_min_height=rel_min_height,
        scale=den_scale,
        alpha=1,
        quantile_lines = TRUE,
        quantiles = c(0.025, 0.5, 0.975)) +
    scale_fill_distiller(direction=1,
                         name="Euclidean distance\nbetween colours") +
    scale_colour_gradient(low = lower_ridgecolor,
                          high = upper_ridgecolor,
                          guide='none') + 
    scale_y_discrete(labels=label_fmt) +
    xlab("Probability of choosing rewarded color") +
    ylab("experiment") +
    xlim(0.5, 1) + 
    coord_cartesian(clip='off') +
    theme_ridges() +
    theme_classic() +
    theme(axis.text.y = element_text(margin = margin(r = 75))) + 
    theme(axis.title.y = element_blank()) + 
    theme(legend.position = "none",
          axis.text = element_text(color='black',
                                   size = 12))

n_exp <- length(exp_dat$exp_unique_name)

cat("plotting circles...\n")
for(k in 1:n_exp){
    verbose <- FALSE
    if (verbose) {
        cat(sprintf("plotting circle for experiment %s",
                    exp_dat$exp_unique_name[k]))
    }

    circle_sugar <- ggplotGrob(ggGradientCircle(
        palette[[exp_dat$color1B[k]]],
        palette[[exp_dat$color1A[k]]],
        fraction1=exp_dat$frac1A[k],
        fineness=50,
        bordersize=0.25))
    
    circle_water <- ggplotGrob(ggGradientCircle(
        palette[[exp_dat$color2B[k]]],
        palette[[exp_dat$color2A[k]]],
        fraction1=exp_dat$frac2A[k],
        fineness=50,
        bordersize=0.25))
    
    ridgeplot <- ridgeplot +
        annotation_custom(circle_sugar,
                          xmin=x_light,
                          xmax=x_light + light_wid,
                          ymin=n_exp + 1 - k - 0.5,
                          ymax=n_exp + 1 - k + 0.5) +
        annotation_custom(circle_water,
                          xmin=x_light + light_pad,
                          xmax=x_light + light_pad + light_wid,
                          ymin=n_exp + 1 - k - 0.5,
                          ymax=n_exp + 1 - k + 0.5)
}

save_plot(outpath,
          ridgeplot,
          base_height=8,
          base_aspect_ratio=.8)
