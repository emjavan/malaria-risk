# Summarize data and plot malaria figures
# Assumes you run on Frontera from inside code folder

# Open libraries and source functions to run below
library(tidyverse)
source("outbreak_analysis_fxns.R")
source("malaria-plotting-code.R")

# File run from parent directory malaria-risk (assumes you opened the .rproj to perform this task on local machine)
#  if this is the case update rl=""
rl = "../"
fig_dir = paste0(rl, "figures/")
if(!dir.exists(fig_dir)){
  dir.create(fig_dir)
}

############################################################################
# Estimate distribution of cumulative annual cases from import=0 simulations
# Returns the data frame of cumulative cases detected
# Writes to "rl/processed_data/full_run_processed_data/clean_epi_files/expect_detect_per_year.csv"
# Only depends on sims with 0 import so if all those have finished successfully you can run

yr_detect_file=paste0(rl, "processed_data/full_run_processed_data/clean_epi_files/expect_detect_per_year.csv")
if(!file.exists(yr_detect_file)){
  print("Creating file")
  est_detect = import_0_cum_cases(run_location=rl, est_detect=100, expect_import=2044)
}else{
  print("File detected and loaded")
  est_detect = read_csv(yr_detect_file)
}

# Plot histogram of the cases
ggplot(est_detect, aes(x=sum_cum_detect))+
  geom_histogram()+
  theme_bw()





#####################################################################
# Organize data for heat map of epi prob by R0, import, and 0-7 cases
# Interpolates the epi prob for needed but missing combinations of R0 and import
# produces all_epi_by_r0_import.csv
epi_prob_by_r0_import(run_location=rl) # , date="2023-07-30"

# Plot heat map of epi prob by R0 and importation across 0-7 cases
p2=epi_r0_import_plot(run_location=rl)

# Save the plot returned by function. Could go in the fuction, but this gives you some ability to look at it and adjust
ggsave(
  paste0(fig_dir, "epi_r0_import_heatmap.png"),
  p2,
  width = 5.25,
  height = 3.25,
  dpi = 1200
)


################################################################################
# Serial operation, looping over counties to get their simulations and summarize
# produces *_all_epi_data_by_case.csv, *_all_epi_prob_by_case.csv, and county_real_epi_prob_and_trigger.csv
process_all_epi_prob_to_mean_df(run_location=rl, date="2023-07-30")

# Plot the summary stats on maps
# If .rda of county shapes not on Frontera, then use tidycensus container. In Emily's $WORK if needed, too big for github
plot_continental_US_epi_and_trigger(run_location=rl)

########################################################################################
# Plot Sarasota County with updated R0 
# include vertical lines of the date Sarasota hit those cumulative case counts
# *12049_all_epi_prob_by_case.csv has the mean epi prob for all cases, min/max have to pull from og files












