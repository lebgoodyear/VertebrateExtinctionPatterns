# Vertebrate Extinction Patterns

Scripts for the analysis of both temporal effects and human population growth on vertebrate exintctions.

Contains the following scripts:

1) conda_env_r_vertex.sh

Bash script that creates a conda environment containing the required R packages for running the rest of the scripts.

2) prep_data.R

Prepares the extinction data and human population data by filtering, binning the data by custom year blocks and combining the two data sets.

3) explore_data.R

Analyses distribution of numbers of extinctions, produces initial plots and calculates correlations.

4) time_analysis.R

Analyses trends in numbers of extinctions through time with Mann-kendall tests and compares numbers of extinctions pre- and post-industrial revolution with Chi-squared tests.

5) modelling.R

Uses GLMs (quasipoisson with log link and also negative binomial) to model relationship between number of extinctions and human population growth. Also takes UN population prediction scenarios and predicts number of extinctions for each scenario.
