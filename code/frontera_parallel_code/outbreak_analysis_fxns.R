

#' Function to join the epidemic risk for all R0 and import combinations to county case data, then summarize by county
#'
#' @param run_location where to look for folders in relation to current working dir run_location="" when opening Rproj in parent dir
#' @param date the date of the files were created, should match date of code/frontera_parallel_code/run-malaria-sims-only_parallel.R
#'
process_all_epi_prob_to_mean_df = function(run_location="../", date="2023-07-22"){
  # Open files
  county_case_data = read_csv(paste0(run_location, "input_data/uscounties_malria_case_2023-07-19.csv")) %>%
    rename(fips=county_fips) %>%
    mutate(fips=as.character(fips),
           fips=ifelse(nchar(fips)==5, fips, paste0(0, fips)) ) # padding all fips to be 5 char
  r0_import_data = read_csv(paste0(run_location, "processed_data/all_rnot_import_per_county_", date, ".csv")) %>%
    left_join(county_case_data, by="fips") %>%
    mutate(fips=ifelse(nchar(fips)==5, fips, paste0(0, fips)) )
  
  # Join files and initialize variables for loop
  all_data = r0_import_data %>%
    mutate(rnot_round2         = as.character(rnot_round2), # have to be characters to match on join
           daily_import_round3 = as.character(daily_import_round3),
           cases    = as.character(cases),
           epi_numer= NA, # initialize the columns for the join
           epi_denom= NA,
           epi_prob = NA) 
  # Getting the expected R0 and import from files, rather than running on what is in the folder
  all_r0 = unique(all_data$rnot_round2)
  all_import = unique(all_data$daily_import_round3)
  folder_path = paste0(run_location, "processed_data/full_run_processed_data/")

  # Loop over all files to join epi_prob by R0, import, cases
  for(i in 1:length(all_r0)){ # 1:length(all_r0)
    for(j in 1:length(all_import)){ # 1:length(all_import)
      open_file = # open the file of epidemic probabilities and clean-up for join
        read_csv(
          paste0(folder_path, "epi_prob_data_", all_r0[i], "_", all_import[j], "_10000_", date, ".csv")) %>%
        rename(cases    = detected,
               epi_prob = prob_epidemic) %>%
        mutate(rnot_round2         = all_r0[i], # doubles have to be characters to join exactly
               daily_import_round3 = all_import[j] ) %>%
        mutate(across(everything(), as.character)) # make everything a character
      
      #
      epi_over_50_temp = open_file %>%
        filter(epi_prob > 0.5) %>%
        slice(1) %>%
        select(detected)
      
      all_data = all_data %>% # join the epi_probs and clean-up extra columns
        left_join(open_file, by = c("rnot_round2", "daily_import_round3", "cases")) %>%
        # the epi_prob we want is the .y joined, so keep it or leave as epi_prob.x
        mutate(epi_numer = ifelse(is.na(epi_numer.y), epi_numer.x, epi_numer.y),
               epi_denom = ifelse(is.na(epi_denom.y), epi_denom.x, epi_denom.y),
               epi_prob = ifelse(is.na(epi_prob.y), epi_prob.x, epi_prob.y) ) %>%
        select(-epi_prob.x, -epi_prob.y, -epi_numer.x, -epi_numer.y, -epi_denom.x, -epi_denom.y)
      
      
      print(paste0("finished R0=", all_r0[i], " import=", all_import[j]))
      print(paste0("R0 iter=", i, " import iter=", j))
      
      # write the full file of 1K R0 per county to test 
      write.csv(all_data, paste0(run_location, "processed_data/full_run_all_data_", date, ".csv"), row.names=F)
      
      
      
    } # end for j over import
  } # end for i over R0
  
  # write the full file of 1K R0 per county and the summarized one 
  write.csv(all_data, paste0(run_location, "processed_data/full_run_all_data_", date, ".csv"), row.names=F)
  
  # Get mean epi prob per county from all 1000 R0 and mean/min/max R0?
  mean_all_data = all_data %>%
    mutate(rnot_round2 = as.double(rnot_round2)) %>%
    group_by(fips, county, county_full, state_id, state_name,
             daily_import_round3, population, cases) %>%
    summarise(rnot_mean = mean(rnot_round2), # , na.rm=T
              rnot_min = min(rnot_round2),
              rnot_max = max(rnot_round2),
              epi_prob_mean = sum(epi_numer)/sum(epi_denom), #mean(epi_prob),
              epi_prob_min  = min(epi_prob),
              epi_prob_max  = max(epi_prob))
  write.csv(mean_all_data, paste0(run_location, "processed_data/full_run_mean_all_data_", date, ".csv"), row.names=F)
  
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



