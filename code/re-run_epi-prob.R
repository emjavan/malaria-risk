# Code to re-run epi_prob calculation so it will output the numerator and denominator as columns
# Only needed if you have simulations but added something new to the output of epi_prob calculation
# This is way too slow to do in serial for all 16K combinations, so use it as a test to see if your code is behaving 

# on frontera always first: module load Rstats
# run from the parent directory malaria-risk
library(tidyverse)
source("code/frontera_parallel_code/outbreak_analysis_fxns.R")
source("code/frontera_parallel_code/sim-malaria-outbreaks.R")

path             = "processed_data/full_run_processed_data/"
run_df = list.files(path = path, pattern = ".rda$") %>%
  as_tibble() %>%
  separate(value, into = c(NA, "base_r_not", "intro_rate", "num_reps", "date"), sep="_") %>%
  mutate(date = substr(date, 1, 10),
         path = path)

run_df %>% 
  pmap(.f = get_save_path) %>% # simply gets the path of simulation file generated above
  map(get_epi_prob) # changed this function name from the zika/covid code to be more logical, fnc creates the epidemic prob table
