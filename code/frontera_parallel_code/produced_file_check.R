# Get the expected and produced count of all files to ensure all files produced
# If something wasn't run, re-run it
# Use inside the code/frontera_parallel_code folder

library(tidyverse)

# Get expected file count
commands_file = read_table("commands_malaria.txt", col_names = F) %>%
  select(X4, X5) %>%
  rename(rnot=X4, import=X5) %>%
  mutate(expect_exists = 1,
         rnot = as.character(rnot),
         import = as.character(import) )
n_expect = nrow(commands_file) # 16315 is expected total epi_prob and sim files

# Get actual file counts
folder_path = "../../processed_data/full_run_processed_data/" #full_run_processed_data
epi_files = list.files(folder_path, pattern = ".csv$") %>%
  enframe(value = "filename") %>%
  select(-name) %>%
  separate(filename, into = c(NA, NA, NA, "rnot", "import", "sims", "date_csv"), 
           sep = "_" ) %>%
  mutate(date_csv = substring(date_csv, 1, 10),
         date_csv = as.Date(date_csv),
         csv_exists = 1,
         rnot = as.character(rnot),
         import = as.character(import) )
sim_files = list.files(folder_path, pattern = ".rda$") %>%
  enframe(value = "filename") %>%
  select(-name) %>%
  separate(filename, into = c(NA, "rnot", "import", "sims", "date_rda"), 
           sep = "_" ) %>%
  mutate(date_rda = substring(date_rda, 1, 10),
         date_rda = as.Date(date_rda),
         rda_exists = 1,
         rnot = as.character(rnot),
         import = as.character(import) )
all_expect_sim_epi = commands_file %>%
  left_join(sim_files, by=c("rnot", "import")) %>%
  distinct() %>%
  left_join(epi_files, by=c("rnot", "import", "sims")) %>%
  select(-sims, -starts_with("date"))
all_expect_sim_epi$rda_exists[is.na(all_expect_sim_epi$rda_exists)]=0
all_expect_sim_epi$csv_exists[is.na(all_expect_sim_epi$csv_exists)]=0

# Check files counts and figure out which files is any need to be processed
n_epi = nrow(epi_files); n_sim = nrow(sim_files)
source("sim-malaria-outbreaks.R")
if( !(n_expect==n_epi) ){
  print(paste0("Epi prob files missing. Only ", n_epi, " when ", n_expect, " expected")) 
  epi_missing = all_expect_sim_epi %>%
    filter(csv_exists==0 & rda_exists==1)
  write.csv(epi_missing, "../../processed_data/missing_epiprobs.csv", row.names = F)
  base_r_not = unique(epi_missing$rnot)
  intro_rate = unique(epi_missing$import)
  path       = folder_path
  date       = "2023-07-22"
  run_df     = expand_grid(base_r_not, intro_rate, path, date)
  num_runs   = 10000
  
  run_df %>% 
    pmap(.f = get_save_path, num_reps = num_runs) %>% # simply gets the path of simulation file generated above
    map(get_epi_prob) # changed this function name from the zika/covid code to be more logical, fnc creates the epidemic prob table
}else{
  print("All epi_probs expected available")
} # end if epi csv does not meet expected count

epi_files_new = list.files(folder_path, pattern = ".csv$")
print(paste0("now total epi risk files are", length(epi_files_new)))

# Confirm if all sim rda exist and if not write what needs to be run to file
if( !(n_expect==n_sim) ){
  print(paste0("Sim files missing. Only ", n_sim, " when ", n_expect, " expected")) 
  sims_missing = all_expect_sim_epi %>%
    filter(rda_exists==0) %>%
    select(rnot, import)
  write.csv(sims_missing, "../../processed_data/missing_simsrda.csv", row.names = F)
}else{
  print("All sims expected available")
} # end if sims rda does not meet expected count


if( !(length(epi_files_new)==n_sim) ){
  print("Epi prob files and total rda still don't match")  
}else{
  print("Sims and epi prob have same num of files")
  
}














