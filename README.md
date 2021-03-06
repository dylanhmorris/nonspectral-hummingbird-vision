# Wild hummingbirds discriminate nonspectral colors
[Mary Caswell Stoddard](https://www.marycstoddard.com/)(1, 2, \*), [Harold N. Eyster](https://orcid.org/0000-0002-5571-3126)(2, 3, †), [Benedict G. Hogan](https://orcid.org/0000-0001-6762-7738 )(1, 2, †), [Dylan H. Morris](https://dylanhmorris.com)(1), [Edward R. Soucy](http://orcid.org/0000-0002-1187-5596)(4), [David W. Inouye](http://biology.umd.edu/david-inouye.html)(2, 5)

\* Corresponding author
† These authors contributed equally to this work. 

1. Department of Ecology and Evolutionary Biology, Princeton University, Princeton, NJ 08544, USA. 
2. Rocky Mountain Biological Laboratory, Crested Butte, CO 81224, USA. 
3. Institute for Resources, Environment and Sustainability, University of British Columbia, Vancouver, British Columbia, Canada V6T 1Z4
4. Center for Brain Science, Harvard University, Cambridge, MA 02138, USA
5. Department of Biology, University of Maryland, College Park, MD 20742, USA


## Repository information
This repository accompanies the paper "Wild hummingbirds discriminate nonspectral colors" (M.C. Stoddard et al). It provides data and code for reproducing statistical analysis and recreating display figures from the paper.

## License and citation information
If you use the code or data provided here, please make sure to do so in light of the project [license](LICENSE) and please cite our work. We have provided [citation guidelines](CITATION.md) for your reference.

## Article abstract 
Many animals have the potential to discriminate nonspectral colors. For humans, purple is the clearest example of a nonspectral color. It is perceived when two color cone types in the retina (blue and red) with nonadjacent spectral sensitivity curves are predominantly stimulated. Purple is considered nonspectral because no monochromatic light (such as from a rainbow) can evoke this simultaneous stimulation. Except in primates and bees, few behavioral experiments have directly examined nonspectral color discrimination, and little is known about nonspectral color perception in animals with more than three types of color photoreceptors. Birds have four color cone types (compared to three in humans) and might perceive additional nonspectral colors such as UV+red and UV+green. Can birds discriminate nonspectral colors, and are these colors behaviorally and ecologically relevant? Here, using comprehensive behavioral experiments, we show that wild hummingbirds can discriminate a variety of nonspectral colors. We also show that hummingbirds, relative to humans, likely perceive a greater proportion of natural colors as nonspectral. Our analysis of plumage and plant spectra reveals many colors that would be perceived as nonspectral by birds but not by humans: Birds’ extra cone type allows them not just to see UV light but also to discriminate additional nonspectral colors. Our results support the idea that birds can distinguish colors throughout tetrachromatic color space and indicate that nonspectral color perception is vital for signaling and foraging. Since tetrachromacy appears to have evolved early in vertebrates, this capacity for rich nonspectral color perception is likely widespread.

## Directories
- ``src``: all code, including model specification, model fitting, post-processing of model output, and figure generation
- ``dat``: data files, in ``.csv`` format. ``dat/cleaned`` contains data files that have been generated from raw files for model fitting.
- ``out``: output files from model analysis, including MCMC chains (as serialized ``.Rds`` datasets), autogenerated figures (as ``.pdf`` files), calculated statistical quantities of interest (as ``.csv`` files), and chain convergence diagnostics (as ``.csv`` files).

## Reproducing analysis

A guide to reproducing the analysis from the paper follows. If you encounter issues, see the **Troubleshooting** section at the end of this README.

All seeds for pseudorandom number generation should be fixed, so results should be consistent on any given machine, so long as seeds and parameters are kept constant. Note, however, that due to differences in computer architecture, there may be minor differences in MCMC output from machine to machine. These should not qualitatively effect results in any way.

### Getting the code
First download this repository. The recommended way is to ``git clone`` it from the command line:

    git clone https://github.com/dylanhmorris/nonspectral-hummingbird-vision.git

Downloading it manually via Github's download button or from [OSF](https://doi.org/10.17605/OSF.IO/5MRKS) should also work.

### Dependency installation
The analysis can be auto-run from the project ``Makefile``, but you may need to install some external dependencies first. See the **Dependency installation guide** below for a complete walkthrough. In the first instance, you'll need a working installation of the statistical programming language R and an implementation of GNU-style ``make``. If you already have these, dependencies can be installed from the command line by typing

    make depend

### Running the analysis

The simplest approach is simply to type ``make`` at the command line, which should produce a full set of figures and MCMC output (saved as R Dataset ``.Rds`` files in the ``out/mcmc-chains/`` directory as ``<model_name>_chains.Rds``). These can be loaded in any working R installation, as long as the package ``rstan`` is also installed.

If you want to do things piecewise, typing ``make <filename>`` for any of the files listed in the ``dat/cleaned`` or ``out`` directories below should run the steps needed to produce that file.

Some shortcuts are available:

- ``make figs`` produces all figures
- ``make chains`` produces all MCMC output.
- ``make clean`` removes all generated files, leaving only source code (though it does not uninstall packages)

### Examining Stan code

Examining the raw Stan code is the place to start to understand how models have been specified. But note that hyperparameters for the prior distributions are set at runtime rather than hard-coded into the ``.stan`` files, so that recompilation is not required when parameter choices are changed (this makes it easier to try the models using different priors, for sensitivity analysis).

Prior hyperparameter choices are specified in the model fitting file, ``src/fit_stan_model.R``.

Note that the ``learning_effect`` variables should be understood to refer to what we call "experience effects" in the text, and ``unbiased_correct_experienced`` is our population average color discrimination parameter.

## Project structure when complete

The project structure once the full analysis has been run should look like this:
```
├── CITATION.md
├── LICENSE
├── Makefile
├── README.md
├── dat
│   └── cleaned
│       └── CleanedDataTable.csv
├── out
│   ├── chain_diagnostics.csv
│   ├── figures
│   │   ├── figure_axis_difference.pdf
│   │   ├── figure_axis_difference_multilevel.pdf
│   │   ├── figure_experience.pdf
│   │   ├── figure_main_stan_model_ridgeline.pdf
│   │   ├── figure_multilevel_model_ridgeline.pdf
│   │   ├── figure_no_error_model_ridgeline.pdf
│   │   ├── figure_position_effect.pdf
│   │   ├── figure_pp_check.pdf
│   │   └── figure_pp_check_all.pdf
│   ├── mcmc_chains
│   │   ├── main_stan_model_chains.Rds
│   │   ├── multilevel_model_chains.Rds
│   │   └── no_error_model_chains.Rds
│   └── reported_quantities.csv
└── src
    ├── chain_diagnostics.R
    ├── extract_reported_quantities.R
    ├── figure_axis_difference.R
    ├── figure_experience.R
    ├── figure_position_effect.R
    ├── figure_posterior_ridgeline.R
    ├── figure_pp_check.R
    ├── fit_stan_model.R
    ├── install_needed_packages.R
    ├── main_stan_model.rds
    ├── main_stan_model.stan
    ├── multilevel_model.rds
    ├── multilevel_model.stan
    ├── no_error_model.rds
    └── no_error_model.stan
```

## Dependency installation guide
You will need a working R installation with the command line interpreter ``Rscript`` (macOS and Linux) or ``Rscript.exe`` (Windows). On mac and Linux, you can check that you have an accessible ``Rscript`` by typing ``which Rscript``at the command line and seeing if one is found.

If you do not have an R installation, you can install it from [the R project website](https://www.r-project.org/) or from the command line using a package manager such as [Homebrew](https://brew.sh/) on macOS or ``apt-get`` on Linux. macOS users may also need to install the macOS "command line tools" by typing ``xcode-select --install`` at a command prompt.

Once R is installed, you can automatically install all other dependencies (including the Hamiltonian Monte Carlo software Stan and its R interface rstan) on most systems using ``make``. In the top level project directory, type the following at the command line:

    make depend

Alternatively, you can run the script ``src/install_needed_packages.R`` manually. 

Note that installing Stan and RStan can be time-consuming, Stan is a large program that must be compiled from source. Some of the packages in the very valuable [tidyverse](https://www.tidyverse.org/) may also take some time to install.

## Data column descriptions

The columns of dataset ``CleanedDataTable.csv`` found in dat/cleaned`` are explained below:

- ``ExperimentNumber``: Enumeration of experiments.
- ``Trial``: Consecutive trial number. Between each trial, the positions of reward (sugar) and non-reward (water) were switched along with the colors presented on the TetraColorTubes.
- ``ChoseSugar``: Number of visits to the rewarded (sugar) color during trial.
- ``ChoseWater``: Number of visits to non-rewarded color during trial.
- ``nHits``: Total number of visits recorded in trial.
- ``cumulativeHitsAtTrialStart``: The cumulative number of visits to the experiment preceeding the current trial.
- ``CorrectSide``: Indicator variable which indicates whether the reward was in position 1 (1) or 2 (-1), used to account for position bias in bird visitation.
- ``PositionSet``: Indicates pairs of trials with swapped color/reward positions.
- ``Euclid``: Estimated Euclidean distance in hummingbird tetrahedral colorspace between rewarded and not-rewarded colors for the current experiment(see article for details).
- ``dtStart``: Date and time of trial initiation.
- ``dtEnd``: Date and time of trial end.
- ``trialDuration``: Duration (minutes) of trial.
- ``cumulativeDurationAtTrialStart``: Cumulative duration (minutes) of all trials in current experment preceeding current trial.
- ``stoppingCondition``: Indicates whether trials during the current experiment were ended after a 25 bird visits regardless of duration, of after a 15 minute duration (see article for details).
- ``ExperimentDisplayName``: Display name for the current experiment as reflected in the article. 
- ``experimentClass``: Indicates the type of experiment. Not-nonspectral colors that we expect birds to discriminate ('spectral'), control experiments with identically colored lights that we expect birds to be unable to discriminate ('null'), and nonspectral colors that are the main focus of the discrimination experiments ('nonspectral', see article for details).
- ``axisID``: Indicates whether the current experiment falls on the UV-green axis of color space (1) or the UV-red axis of color space (2), or neither of these (NA, see article for details).

## Troubleshooting
- RStan is extremely memory-hungry, and in testing in lower memory environments, we sometimes encountered the following error after an MCMC chain had been sampled but before it could be saved:

```
Error in FUN(X[[i]], ...) : 
  trying to get slot "mode" from an object of a basic class ("NULL") with no slots
Calls: stan ... sampling -> sampling -> .local -> sapply -> lapply -> FUN
Execution halted
```
This is not due to an error in the source code, and it appears to simply be due to lack of memory. It can usually be avoided by opening a fresh command prompt or if necessary restarting your machine.

