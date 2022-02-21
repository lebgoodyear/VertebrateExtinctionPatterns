#!/bin/bash
# #SBATCH --time=03:00:00
# #SBATCH --partition=k2-hipri
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --mem=10G

# Script: conda_env_r_vertex.sh

# Author: Luke Goodyear lgoodyear01@qub.ac.uk
# Date created: Feb 2022
# Last edited: Feb 2022

# Description:
# 1) cleans conda environment 
# 2) creates new environment and installs packages
# 3) cleans conda environment
# Arguments: None


#echo "Loading anaconda..."
#module load apps/anaconda3/2021.05/bin

# remove any unused packages and caches
conda clean -a -y

echo "Creating new R conda environment..."
conda create -n r_vertex -c conda-forge r-base=4.1.0 -y
# activate new environment
source activate r_vertex

# install required packages
conda install -c conda-forge r-sf=1.0_4 -y
conda install -c conda-forge r-rnaturalearth=0.1.0 -y
conda install -c conda-forge r-rnaturalearthdata=0.1.0 -y
conda install -c conda-forge r-ggplot2=3.3.5 -y
conda install -c conda-forge r-dplyr=1.0.7 -y
conda install -c conda-forge r-rcolorbrewer=1.1_2 -y
conda install -c conda-forge r-viridislite=0.4.0 -y
conda install -c conda-forge r-tidyr=1.1.4 -y
conda install -c conda-forge r-gridextra=2.3 -y
conda install -c conda-forge r-spatialpack=0.3 -y 
conda install -c conda-forge r-mass=7.3_55 -y

# remove any unused packages and caches after installations
conda clean -a -y

# view packages in environment as check
conda list

# deactivate environment
conda deactivate


# end of script
