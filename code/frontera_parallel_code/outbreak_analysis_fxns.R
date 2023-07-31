################################################
# Functions to analyze the epidemic simulations
################################################


#' Function to join the epidemic risk for all R0 and import combinations to county case data, then summarize by county
#'
#' @param run_location where to look for folders in relation to current working dir run_location="" when opening Rproj in parent dir
#' @param date the date of the files were created, should match date of code/frontera_parallel_code/run-malaria-sims-only_parallel.R
#'
#'
process_all_epi_prob_to_mean_df = function(run_location="../", date="2023-07-22"){
  # Open files
  # Case data for each county: This file created by Shraddah 
  county_case_data = read_csv(paste0(run_location, "input_data/uscounties_malria_case_2023-07-19.csv")) %>%
    rename(fips=county_fips) %>%
    mutate(fips=as.character(fips),
           fips=ifelse(nchar(fips)==5, fips, paste0(0, fips)) ) # padding all fips to be 5 char
  county_case_data_w_epi = county_case_data %>% # initalize data frame where the values for making maps will be stored
    mutate(sim_run_date          = NA,
           epi_prob              = NA,
           case_when_epi_over_50 = NA)
  # R0 and import data: This file created by code/clean_r0_import_write_commands.R
  r0_import_data = read_csv(paste0(run_location, "processed_data/all_rnot_import_per_county_", date, ".csv")) %>%
    mutate(fips=ifelse(nchar(fips)==5, fips, paste0(0, fips)) ) %>%
    group_by(fips, rnot_round2, daily_import_round3) %>%
    summarize(total_sims = n()) %>% # only keep the count of unique combinations 
    ungroup() %>%
    mutate(across(everything(), as.character)) %>%
    filter(!grepl('^02|^15', fips)) # remove Alaska and Hawaii as they didn't have importation risk estimated
  
  # Largest group of unique R0 and import combinations is 120, so doing the join and calucation by county is fine
  # Need all 101 case epi numer/denom per county to get the correct epi prob and cases when epi_prob>0.5
  # grp_sz = r0_import_data %>%
  #   group_by(fips) %>%
  #   summarise(size_of_group = n())
  
  # Looks like FIPS 51678 is missing from the R0 data - Emily J. 2023-07-31
  #sort(county_case_data$fips[!(county_case_data$fips %in% r0_import_data$fips)])
  
  folder_path = paste0(run_location, "processed_data/full_run_processed_data/")
  output_folder_path = paste0(folder_path, "clean_epi_files/")
  files = dir(folder_path, pattern = "*.csv$") # get file names
  if(!dir.exists(output_folder_path)){dir.create(output_folder_path)}
  all_county_fips = unique(r0_import_data$fips)
  
  for(i in 1:length(all_county_fips)){ # 3106, expecting 3107
    temp_county = r0_import_data %>% # filter to only one county
      filter(fips == all_county_fips[i])
    
    # Join all the epi data by R0 and import for the county filtered above
    # I don't love this implementation, bc we have to open the files so many times
    # Can't think of something more robust, trying to avoid nested for loops!
    epi_data = tibble(filename = files) %>% # create a data frame
      separate(filename, sep="_", remove = F,
               into=c(NA, NA, NA, "rnot_round2", "daily_import_round3", "sims", "date")) %>%
      mutate(fips = all_county_fips[i],
             across(everything(), as.character)) %>%
      filter( rnot_round2 %in% temp_county$rnot_round2 ) %>%
      filter( daily_import_round3 %in% temp_county$daily_import_round3 ) %>%
      mutate(date = substr(date, 1, 10),
             file_contents = map(filename, ~read_csv(file.path(folder_path, .))) ) %>%
      unnest(cols = c(file_contents)) %>% # opens up the joined file contents into new columns in the tibble
      rename(cases= detected) %>%
      select(fips, everything(), -prob_epidemic) %>% # the true prob_epidemic will be calculated below
      left_join(temp_county, by=c("fips", "rnot_round2", "daily_import_round3")) %>%
      mutate(total_sims = as.numeric(total_sims))
    # write the full data for epi prob calculation to a file
    write.csv(epi_data, paste0(output_folder_path, all_county_fips[i], "_all_epi_data_by_case.csv"), row.names = F)
    
    # Get epi_prob by all cases for the county
    epi_sum = epi_data %>%
      mutate(total_sims = as.numeric(total_sims)) %>%
      group_by(fips, date, cases) %>%
      summarise(epi_prob = sum(epi_numer*total_sims)/sum(epi_denom*total_sims)) %>%
      ungroup() %>%
      rename(sim_run_date = date)
    write.csv(epi_data, paste0(output_folder_path, all_county_fips[i], "_all_epi_prob_by_case.csv"), row.names = F)
    
    # Trigger threshold of cases when the epi prob is over 0.5
    epi_over_50 = epi_sum %>%
      filter(epi_prob > 0.5) %>%
      arrange(epi_prob) %>%
      slice(1) %>%
      rename(case_when_epi_over_50 = cases) %>%
      select(-epi_prob)
    
    county_case_data_w_epi = county_case_data_w_epi %>%
      left_join(epi_sum, by=c("fips", "cases") ) %>%
      mutate(epi_prob = ifelse(is.na(epi_prob.y), epi_prob.x, epi_prob.y),
             sim_run_date = ifelse(is.na(sim_run_date.y), sim_run_date.x, sim_run_date.y)) %>%
      select(-ends_with(".x"), -ends_with(".y")) %>%
      left_join(epi_over_50, by=c("fips", "sim_run_date")) %>%
      mutate(case_when_epi_over_50 = ifelse(is.na(case_when_epi_over_50.y), 
                                            case_when_epi_over_50.x, case_when_epi_over_50.y)) %>%
      select(-ends_with(".x"), -ends_with(".y"))
    
    # write the updated file every iteration in case it crashes or something, we can keep whatever finished
    write.csv(county_case_data_w_epi, paste0(output_folder_path, "county_real_epi_prob_and_trigger.csv"), row.names = F)
    
  } # end for i loop over counties in the continental US
} # end function process_all_epi_prob_to_mean_df



#######################################################################################
# Code from rtZIKVrisk package outbreak_analysis_fxns.R
#######################################################################################

#' Get single final cumulative infected numbers from sims
#'
#' @param x a single trial run
#' @return the max cumulative infected individuals for each simulation
#' @examples
#'
last_cuminfect_local <- function(x) {
  x[nrow(x), "Cumulative_Infections"] - x[nrow(x), "Cumulative_Intro_Infections"]
}

#' Get all final cumulative infected numbers from sims
#'
#' @param trials A list with simulated zika outbreaks. Usually utput of run_n_zika_sims()
#' @return A vector of all the max cumulatively infected individuals for each simulation.
#' Calls \code{\link{last_cuminfect_local}}
#' @examples
#'
all_max_cum_infect <- function(trials) {
  
  unlist(plyr::laply(trials, last_cuminfect_local))
}

#' Get single maximum of local prevalence
#'
#' @param x a single trial run
#' @return the maximum local prevalence from run
#' @examples
#'
max_local_prev <- function(x){
  return(max(x[,"Total_Infections"]-x[,"Total_Intro_Infections"]))
}


#' Returns cumulative vector of locally detected cases
#'
#' @param x a single trial run
#' @return Cumulative detectections in vector that are local
#' @examples
#'
cum_detect_local <- function(x){
  ## Returns column of cumulative local detections
  x[, "Cum_Detections"] - x[, "Cum_Intro_Detections"]
}



#' Cumulative cases for detects for single simulation run
#'
#' @param df A dataframe that is the result of a single zika simulation run
#' @param max_detect the maximum number of detections desired to be analyzed
#' @return the maximum cumulative cases and max prevalence for each unique detection
cumcases_by_detects <- function(df, max_detect){
  all_detects <- cum_detect_local(df)
  unique_detects <- unique(all_detects)
  unique_detects <- unique_detects[unique_detects<=max_detect]
  
  data.frame(detected = unique_detects, cum_infections = last_cuminfect_local(df), max_prevalence = max_local_prev(df))
}

#' Get all cumcases and prevalence by detects
#'
#' Vectorized version of \code{\link{cumcases_by_detects}}, to get all trial information.
#'
#' @param trials a list of zika trial simulations
#' @param max_detect the maximum amount of detections interested in for analysis
#'
#' @return a single probability of epidemic for single detected value
get_cumcases_by_detects_all <- function(trials, max_detect){
  ## Returns data frame of all prevalence by detections for all trials
  plyr::ldply(trials, cumcases_by_detects, max_detect)
}


#' Get frequency cases above cum threshold and prev threshold for unique detected
#'
#' Returns the frequency that cases were found to be above a threshold.
#'
#' @param df dataframe formatted output from \code{\link{get_cumcases_by_detects_all}}. Has one columns for detected, cum_infections, and max_prevalence.
#' @param detected single integer number of interest
#' @param cum_threshold cumulative threshold for epidemic classification
#' @param prev_threshold prevalence threshold for epidemic classification
#' @param num_necessary number of trials necessary for epi classification usually 1\% of trials.
#'
#' @return a single probability of epidemic for single detected value.
freq_above_thresh <- function(df, detected, cum_threshold, prev_threshold, num_necessary){
  ## Takes in dataframe of all prevalence by detect
  ## Returns a single frequency of times that
  ## prevalence for a specific detection criteria is below a threshold
  rows <- which(df[,"detected"] == detected)
  if(length(rows)<=num_necessary){
    return( c(NA, NA, NA))
  }else{
    ## Return number of rows that excede both thresholds divided by the total rows
    # Edited by Emily Javan on 07-28-23 to get mean epi prob over a vector of R0 for a county (see malaria-risk)
    epi_numer = sum(df[rows, "cum_infections"] >= cum_threshold & df[rows,"max_prevalence"] >= prev_threshold)
    epi_denom = length(rows)
    epi_prob  = epi_numer/epi_denom
    df = data.frame(epi_numer=epi_numer, epi_denom=epi_denom, epi_prob=epi_prob)
    return(df)
  }
}

#' Vectorizes the freq_above_thresh for the number of detected
#'
#' Called from the get_epidemic_prob_by_d
#'
#' @param detected vector of detected of interest
#' @inheritParams freq_above_thresh
#' @return A vector of epidemic probs for the given detected values.
#'
freq_above_thresh_vec <- Vectorize(freq_above_thresh, vectorize.args = "detected")


#' Get the epidemic probability as a function of number of reported cases
#'
#' @param trials A list with simulated zika outbreaks. Usually utput of \code{\link{run_n_zika_sims}}
#' @param prev_threshold The maximum autochthonous prevalence necessary to be classifed as an epidemic - depends on individual scenario run values
#' @param cum_threshold The cumulative autochthonous infections necessary to be classified as an epidemic - usually the e_thresh value of the runs
#' @param max_detect The maximum number of detections to go to.
#' @param num_necessary Number of instances of trials to be necessary before it gets a probability. Usually want to be ~1\% of total runs
#'
#' @return A dataframe of rows max_detect+1 that has column for detected, and column for prob_epidemic given that detected number
#' @export
#' @examples
#' \dontrun{
#' get_epidemic_prob_by_d(trials, 5, 100, 15, 1)
#' }
get_epidemic_prob_by_d <- function(trials, prev_threshold, cum_threshold, max_detect=50, num_necessary=1){
  detected <- seq(0, max_detect)
  ## Gets the number of cumulative cases for each run for each unique detected.
  data <- get_cumcases_by_detects_all(trials, max_detect = max_detect)

  ## Calculates the frequency that the cum_threshold and prev_thresholds are met in for each detected number
  probs <- freq_above_thresh_vec(data, detected, cum_threshold, prev_threshold, num_necessary) %>%
    as_tibble() %>%
    t() %>%
    as.data.frame()
  names(probs) = c("epi_numer", "epi_denom", "prob_epidemic")
  probs = probs %>%
    mutate(detected = detected) %>%
    select(detected, everything()) %>%
    mutate(across(everything(), as.character)) # convert all columns to characters since one is being saved as list
  return(probs)
}



