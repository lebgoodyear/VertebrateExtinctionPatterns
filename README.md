# Vertebrate Extinction Patterns

Scripts for the analysis of the effect of human population density on vertebrate exintctions.

Contains the following scripts:

1) conda_env_r_vertex.sh

Bash script that creates a conda environment containing the required R packages for running the rest of the scripts.

2) prep_data.R

Prepares the extinction and human population data by filtering, calculating human population density change, binning the data by custom year blocks and combining the two data sets.

3) explore_data.R

Analyses distribution of numbers of extinctions, produces initial plots and calculates correlations.

4) spatial_analysis.R

Produces heatmaps of human population density overlaid with extinction point data for custom periods in time.

5) modelling.R

Uses GLMs (quasipoisson with log link and also negative binomial) to model relationship between number of extinctions and human population density/change. Also takes UN population prediction scenarios and predicts number of extinctions for each scenario.
