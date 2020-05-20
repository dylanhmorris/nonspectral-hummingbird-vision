#!/usr/bin/env Rscript

#####################################
## name: install_needed_packages.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## installs needed packages for
## reproducing hummingbird non-spectral
## color perception /learning
##
##
####################################

install_if_absent <- function(package_name){
    if (!suppressPackageStartupMessages(
             require(package_name, character.only=TRUE))){
      install.packages(pkgs=package_name,
                       repos="http://cloud.r-project.org")
  }
  else
      cat(sprintf("Package %s already installed\n", package_name))
}

needed_packages <- c(
    "rstan",
    "dplyr",
    "readr",
    "magrittr",
    "stringr",
    "bayesplot",
    "tidybayes",
    "ggplot2",
    "ggthemes",
    "ggridges",
    "cowplot",
    "latex2exp")

for (package in needed_packages)
    install_if_absent(package)
