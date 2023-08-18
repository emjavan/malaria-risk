####################################################################
## Code to run and save all Malaria simulations 2023-07-22s
####################################################################
# Run from inside the code/frontera_parallel_code/ dir 
# When on Frontera submit job with: sbatch launcher_malaria.sh

# Frontera/TACC's version of R has an issue downloading tidycensus, so must use a container
# we will not do any map/figure plotting in this file, only do simulations and generate epi prob per case detected

library(tidyverse)
#library(rtZIKVrisk) # Note no longer pulling from the original package
source("outbreak_analysis_fxns.R") # edited code from rtZIKVrisk to save numerator and denominator of epidemic probability
source("sim-malaria-outbreaks.R")

## Setup all of the parameters that need to be run
args             = commandArgs(TRUE)
base_r_not       = as.double(args[1]) # R0
intro_rate       = as.double(args[2]) # Intro_rate
path             = "../../processed_data/full_run_processed_data/"
date             = "2023-08-15"  #Sys.Date() # run/save with one system date even if the job runs over mulitple days
run_df           = expand_grid(base_r_not, intro_rate, path, date)
num_runs         = 10000

## Run and save simulations across all parameter combinations
if(!dir.exists(path)){
  dir.create(path)
}

run_df %>% # pipe the 2 inputs into save_malaria_runs function, when refresh is FALSE it will not overwrite an existing output
  pmap(.f = save_malaria_runs, num_reps = num_runs, refresh=TRUE) %>% 
  unlist()

run_df %>% 
  pmap(.f = get_save_path, num_reps = num_runs) %>% # simply gets the path of simulation file generated above
  map(get_epi_prob) # changed this function name from the zika/covid code to be more logical, fnc creates the epidemic prob table
  



