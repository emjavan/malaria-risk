####################################################################
## Code to run and save all Malaria simulations 2023-07-07s
####################################################################
# Run from inside the code/frontera_parallel_code/ dir 
# When on Frontera submit job with: sbatch launcher_malaria.sh

# Frontera/TACC's version of R has an issue downloading tidycensus, so must use a container
# we will not do any map/figure plotting in this file, only do simulations and generate epi prob per case detected

library(tidyverse)
library(rtZIKVrisk) # Note generated when loading package: Please note that 'maptools' will be retired during 2023

source("sim-malaria-outbreaks.R")

## Setup all of the parameters that need to be run
args             = commandArgs(TRUE)
base_r_not       = as.double(args[1]) # R0
intro_rate       = as.double(args[2]) # Intro_rate
path             = "../../processed_data/full_run_processed_data/"
date             = Sys.Date() # run with one system date even if the job runs over mulitple days
run_df           = expand_grid(base_r_not, intro_rate, path, date)
num_runs         = 10000

## Run and save simulations across all parameter combinations
if(!dir.exists(path)){
  dir.create(path)
}
start_time <- Sys.time()
run_df %>% # pipe the 2 inputs into save_malaria_runs function, when refresh is FALSE it will not overwrite an existing output
  pmap(.f = save_malaria_runs, num_reps = num_runs, refresh=TRUE) %>% 
  unlist()
end_time <- Sys.time()
total_time = end_time - start_time
print(paste0("total run time for R0=", base_r_not, " intro_rate=", intro_rate, " num_reps=", num_runs,
             " is ", round(total_time, 2), " sec" ))

run_df %>% 
  pmap(.f = get_save_path, num_reps = num_runs) %>% # simply gets the path of simulation file generated above
  map(get_epi_prob) # changed this function name from the zika/covid code to be more logical, fnc creates the epidemic prob table
  



