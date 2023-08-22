###################################################################
# Functions to analyze the epidemic simulations/probability files
###################################################################

#' Function to estimate cumulative detected cases in one year
#' 0. Process sims from all R0 and import=0 to get 365th or last row of df
#' 1. Distribute X=2044 importations to counties based on daily importation prob
#' 2. For each county, draw n total importation simulations from its set of 1K R0s and prob next import = 0
#' 2b. This requires choosing an R0, then opening that set of 10K sims and sampling 1
#' 3. Sum together the cumulative detections of sims for each importation
#' 
#' @param run_location ="../" if running inside the code folder
#' @param est_detect number of times to estimate cumulative detections in a year
#' @param expect_import total infectious malaria importations to the US, 2044 is for 2023. 
#'   Use 2044*4=8176 if you think we only detect 25% of importations and true value is 4x observed
#'
import_0_cum_cases = function(run_location="../", est_detect=100, expect_import=2044){
  require(tidyverse)
  # Define folder paths, get total files, make needed folders
  folder_path = paste0(run_location, "processed_data/full_run_processed_data/")
  files = dir(folder_path, pattern = "*_0_10000_2023-08-15.rda$") # get file names
  file_df = as_tibble(files) %>%
    rename(rda_name = value) %>%
    mutate(rda_exists = 1)
  output_folder_path = paste0(folder_path, "sim_last_or_year_end/")
  if(!dir.exists(output_folder_path)){dir.create(output_folder_path)}
  if(!dir.exists( paste0(folder_path, "clean_epi_files/") )){
    dir.create( paste0(folder_path, "clean_epi_files/") )}
  output_files = dir(output_folder_path, pattern = "*_0_10000_2023-08-15.csv$")  %>% # get file names
    as_tibble() %>%
    rename(csv_name = value) %>%
    mutate(rda_name = gsub("csv", "rda", csv_name),
           rda_name = gsub("cum_detect_all_", "", rda_name)) %>%
    right_join(file_df, by="rda_name") %>%
    mutate(csv_missing = ifelse(is.na(csv_name), 1, 0) )
  
  n_missing_csv = sum(output_files$csv_missing) 
  if(n_missing_csv > 0){ # check num of missing csv files
    # Function to choose row to slice from each simulation's data frame. Each row is day
    slice_last_or_year_end = function(df){
      total_row = nrow(df)
      if(total_row>365){ # Only want expected cases detected in a year so no sim day > 365
        end_df = df %>%
          slice(365)
      }else{
        end_df = df %>%
          tail(n=1) # get the last row of data frame
      } # end if
    } # end function
    
    # Open simulation files to get the last row or the 365th day of data frame
    missing_files = output_files %>%
      filter(csv_missing==1)
    total_files = nrow(missing_files)
    for(i in 1:total_files){
      load(paste0(folder_path, missing_files$rda_name[i])) # opens list of data frames of all 10K sims
      sim_tail_detect = map(sims, ~slice_last_or_year_end(.) ) %>%
        map_dfr(~as_tibble(.)) #%>%
      # select(Cum_Detections) # going to write entire last row of all dataframes in 10K
      csv_name = gsub("rda", "csv", files[i])
      write.csv(sim_tail_detect, paste0(output_folder_path, "cum_detect_all_", csv_name), row.names=F )
      print(paste0("Finished ", i, " of ", total_files))
    } # end for i
  } # end if any csv missing for analysis below
  
  # All the cumulative detection files once missing files replaced
  cum_detect_files = dir(output_folder_path, pattern = "*.csv$")  %>% # new file set to look through for cumulative detections
    as_tibble() %>%
    rename(filename=value) %>%
    separate(filename, into = c(NA, NA, NA, NA, "rnot_round1", "import", "sims", NA), sep="_", remove=F) %>%
    mutate(rnot_round1 = as.character(rnot_round1))
  # All the importation probabilities
  impProbByCounty_new = read_csv("../input_data/malaria_probability_next_import_comes_from_county_2023-08-14.csv") %>%
    mutate(fips=as.character(fips),
           fips=ifelse(nchar(fips)==5, fips, paste0(0, fips)) )
  # All county R0 samples
  load(paste0(run_location, "input_data/county-single-rnot-estimates_2023-08-09.rda"))
  county_single_rnot_samples_new = county_single_rnot_samples %>% 
    mutate(rnot_round1 = round(rnot, 1),
           fips=as.character(fips),
           fips=ifelse(nchar(fips)==5, fips, paste0(0, fips)),
           rnot_round1 = as.character(rnot_round1))
  
  # Loop for the number of times we want to estimate the cumulative detected cases in a year
  expect_detect_per_year = data.frame()
  for(j in 1:est_detect){ # 100 is just a test, probably do 1K?
    # Randomly sample N=2044 expected importations for 2023 to counties based on their probability of next import
    n_county_imports = sample(impProbByCounty_new$fips, size = expect_import, 
                              prob=impProbByCounty_new$ProbNextOcc, replace=T) %>%
      as_tibble() %>%
      rename(fips=value) %>%
      group_by(fips) %>%
      summarise(total_import = n()) %>%
      ungroup()
    
    # Draw total_import simulations from county i’s set based on the 1K R0s and prob next import = 0
    sub_county = county_single_rnot_samples_new %>% # all R0 samples for all counties
      filter(fips %in% n_county_imports$fips) %>% # get only fips in the sample
      left_join(n_county_imports, by="fips") %>% # need to randomly sample n imports per county in sample
      group_by(fips) %>%
      mutate(samp = sample(n())) %>% # new column of random order of integers 1:nrow(group)
      filter(samp <= total_import) %>% # pick just the rows where samp is less than or equal to total_imports
      ungroup() %>%
      left_join(cum_detect_files, by="rnot_round1") %>% # Open the file randomly chosen for each county to get last row/day of sim
      mutate(file_contents = map(filename, ~(read_csv(file.path(output_folder_path, .)) %>% slice_sample(n=1)) )) %>%
      unnest(cols = c(file_contents))
    
    # Sum the total autochthonous detected infections in the US after one year (365d) assuming all simulations start on day 1
    sum_sub_county = sub_county %>%
      summarise(sum_cum_detect = sum(Cum_Detections - Cum_Intro_Detections), # value we're interested
                sum_cum_inf = sum(Cumulative_Infections - Cumulative_Intro_Infections), # should be ~4x larger if we only detect 25% of all infections
                sum_cum_intro_inf = sum(Cumulative_Intro_Infections)) %>%# just a check, should be 1 for each simulation
      mutate(true_prop_auto_inf_detect = sum_cum_detect/sum_cum_inf)
    expect_detect_per_year = rbind(expect_detect_per_year, sum_sub_county)
    
    # Redundant to write in loop, but don't want to lose progress in case job times out
    write.csv(expect_detect_per_year, paste0(folder_path, "clean_epi_files/expect_detect_per_year.csv"), row.names=F)
  } # end for j
  
  return(expect_detect_per_year)
} # end function


#' Function to join the epidemic risk for all R0 and import combinations to county case data, then summarize by county
#'
#' @param run_location where to look for folders in relation to current working dir run_location="" when opening Rproj in parent dir
#' @param date the date of the files were created, should match date of code/frontera_parallel_code/run-malaria-sims-only_parallel.R
#'
process_all_epi_prob_to_mean_df = function(run_location="../", date="2023-08-15"){
  require(tidyverse)
  # Open files
  # Case data for each county: This file created by Shraddah 
  county_case_data = read_csv(paste0(run_location, "input_data/uscounties_malria_case_2023-08-18.csv")) %>%
    rename(fips=county_fips) %>%
    mutate(fips=as.character(fips),
           fips=ifelse(nchar(fips)==5, fips, paste0(0, fips)) ) %>% # padding all fips to be 5 char
    filter(!grepl('^02|^15', fips)) %>% # remove Alaska and Hawaii as they didn't have importation risk estimated
    filter(!(fips=="11001")) # remove Washington DC
  county_case_data_w_epi = county_case_data %>% # initalize data frame where the values for making maps will be stored
    mutate(sim_run_date = NA, min_R0 = NA, max_R0   = NA,  mean_R0 = NA,
           epi_prob     = NA, case_when_epi_over_50 = NA)
  # R0 and import data: This file created by code/clean_r0_import_write_commands.R
  r0_import_data = read_csv(paste0(run_location, "processed_data/all_rnot_import_per_county_2023-08-09.csv")) %>%
    mutate(fips=ifelse(nchar(fips)==5, fips, paste0(0, fips)) ) %>%
    group_by(fips, rnot_round1, daily_import_round3) %>%
    summarize(total_sims = n()) %>% # only keep the count of unique combinations 
    ungroup() %>%
    mutate(across(everything(), as.character)) %>%
    filter(!grepl('^02|^15', fips)) # remove Alaska and Hawaii as they didn't have importation risk estimated

  write.csv(r0_import_data, paste0(run_location, "processed_data/summarized_rnot_import_per_county_", date, ".csv"),
            row.names = F)
  
  folder_path = paste0(run_location, "processed_data/full_run_processed_data/")
  output_folder_path = paste0(folder_path, "clean_epi_files/")
  files = dir(folder_path, pattern = "*.csv$") # get file names
  if(!dir.exists(output_folder_path)){dir.create(output_folder_path)}
  all_county_fips = unique(r0_import_data$fips)
  
  # County FIPS and location in vector
  # which(all_county_fips=="12115") # position 348 # Sarasota County
  # which(all_county_fips=="12081") # position 330 # Manatee County - neighbor
  # which(all_county_fips=="12027") # position 303 # DeSoto County - neighbor
  # which(all_county_fips=="12015") # position 298 # Charlotte County - neighbor
  # which(all_county_fips=="12049") # position 314 # Hardee County - closest county with max R0>1

  for(i in 1:length(all_county_fips)){ # expecting 3107
    print(paste0("Starting ", i))
    temp_county = r0_import_data %>% # filter to only one county
      filter(fips == all_county_fips[i])
    
    # Join all the epi data by R0 and import for the county filtered above
    # I don't love this implementation, bc we have to open the files so many times
    # Can't think of something more robust, trying to avoid nested for loops!
    epi_data = tibble(filename = files) %>% # create a data frame
      separate(filename, sep="_", remove = F,
               into=c(NA, NA, NA, "rnot_round1", "daily_import_round3", "sims", "date")) %>%
      mutate(fips = all_county_fips[i],
             across(everything(), as.character)) %>%
      filter( rnot_round1 %in% temp_county$rnot_round1 ) %>%
      filter( daily_import_round3 %in% temp_county$daily_import_round3 ) %>%
      mutate(date = substr(date, 1, 10),
             file_contents = map(filename, ~read_csv(file.path(folder_path, .))) ) %>%
      unnest(cols = c(file_contents)) # opens up the joined file contents into new columns in the tibble
      
    if(nrow(epi_data)>0){
      epi_data = epi_data %>%
        mutate(epi_denom = ifelse(is.na(epi_denom), 0, epi_denom),
               epi_numer = ifelse(is.na(epi_numer), 0, epi_numer)) %>%
        rename(cases= detected) %>%
        select(fips, everything(), -prob_epidemic) %>% # the true prob_epidemic will be calculated below
        left_join(temp_county, by=c("fips", "rnot_round1", "daily_import_round3")) %>%
        mutate(total_sims = as.numeric(total_sims))
      # write the full data for epi prob calculation to a file
      write.csv(epi_data, paste0(output_folder_path, all_county_fips[i], "_all_epi_data_by_case.csv"), row.names = F)
      
      # Get epi_prob by all cases for the county
      epi_sum = epi_data %>%
        mutate(total_sims          = as.numeric(total_sims),
               rnot_round1         = as.numeric(rnot_round1),
               daily_import_round3 = as.numeric(daily_import_round3)) %>%
        group_by(fips, date, cases) %>%
        # removing NAs for now to just get best guess with the sims we have completed, still waiting on many to finish
        summarise(mean_R0  = sum(rnot_round1*total_sims)/sum(total_sims),
               max_R0   = max(rnot_round1),
               min_R0   = min(rnot_round1),
               #import   = mean(daily_import_round3),
               epi_prob = sum(epi_numer*total_sims, na.rm = T)/sum(epi_denom*total_sims, na.rm = T)) %>%
        ungroup() %>%
        rename(sim_run_date = date) %>%
        mutate(epi_prob = ifelse(is.nan(epi_prob), 0, epi_prob))
      
      write.csv(epi_sum, paste0(output_folder_path, all_county_fips[i], "_all_epi_prob_by_case.csv"), row.names = F)
      
      # Trigger threshold of cases when the epi prob is over 0.5
      epi_over_50 = epi_sum %>%
        filter(epi_prob > 0.5) %>%
        arrange(epi_prob) %>%
        slice(1) %>%
        rename(case_when_epi_over_50 = cases) %>%
        select(fips, sim_run_date, case_when_epi_over_50)
      
      county_case_data_w_epi = county_case_data_w_epi %>%
        mutate(cases = as.character(cases)) %>%
        left_join(epi_sum %>% mutate(cases = as.character(cases)), by=c("fips", "cases") ) %>%
        mutate(sim_run_date = ifelse(is.na(sim_run_date.y), sim_run_date.x, sim_run_date.y),
               min_R0 = ifelse(is.na(min_R0.y), min_R0.x, min_R0.y),
               max_R0 = ifelse(is.na(max_R0.y), max_R0.x, max_R0.y),
               mean_R0 = ifelse(is.na(mean_R0.y), mean_R0.x, mean_R0.y),
               epi_prob = ifelse(is.na(epi_prob.y), epi_prob.x, epi_prob.y)) %>%
        select(-ends_with(".x"), -ends_with(".y")) %>%
        left_join(epi_over_50, by=c("fips", "sim_run_date")) %>%
        mutate(case_when_epi_over_50 = ifelse(is.na(case_when_epi_over_50.y), 
                                              case_when_epi_over_50.x, case_when_epi_over_50.y)) %>%
        select(-ends_with(".x"), -ends_with(".y")) %>%
        distinct()
        
      # write the updated file every iteration in case it crashes or something, we can keep whatever finished
      write.csv(county_case_data_w_epi, paste0(output_folder_path, "county_real_epi_prob_and_trigger.csv"), row.names = F)
    } # end if to check there is data to analyze

    print(paste0("Finished ", i))
  } # end for i loop over counties in the continental US
  return("All county epi files summarized successfully")
} # end function process_all_epi_prob_to_mean_df

# Function to organize data for heat map of epi prob by R0, import, and cases
#'
#' @param run_location where to look for folders in relation to current working dir run_location="" when opening Rproj in parent dir
#'
epi_prob_by_r0_import = function(run_location="../"){ # , date="2023-07-30"
  
  folder_path = paste0(run_location, "processed_data/full_run_processed_data/")
  output_folder_path = paste0(folder_path, "clean_epi_files/")
  files = dir(folder_path, pattern = paste0("*.csv$") ) # get file names
  if(!dir.exists(output_folder_path)){dir.create(output_folder_path)}
  
  # Open/create file of 0-7 cases for all simulations available
  file_name = paste0(output_folder_path, "all_epi_by_r0_import_og.csv")
  if(!file.exists(file_name)){
    all_epi_data = tibble(filename = files) %>% # create a data frame
      separate(filename, sep="_", remove = F,
               into=c(NA, NA, NA, "rnot_round1", "daily_import_round3", NA, NA)) %>%
      mutate(file_contents = map(filename, ~(read_csv(file.path(folder_path, .)) %>%
                                               slice(1:8)) )) %>% # get cases 0-7 (Sarasota cty max)
      unnest(cols = c(file_contents))
    
    write.csv(all_epi_data, file_name, row.names=F)
  }
  # else{
  #   all_epi_data = read_csv(file_name)
  # } # end if else
  # 
  # # From all the CSV files we have epi data, expand the grid to get all possible combinations
  # rnot_round1 = seq(min(all_epi_data$rnot_round1), max(all_epi_data$rnot_round1), 0.1)
  # daily_import_round3 = seq(min(all_epi_data$daily_import_round3), max(all_epi_data$daily_import_round3), 0.001)
  # detected = seq(0, 7, 1) # get cases 0-7 (Sarasota cty max)
  # every_combo = expand_grid(rnot_round1, daily_import_round3, detected) %>%
  #   mutate(across(everything(), as.character))
  # 
  # 
  # # below is failing because the R0=0 with detect=0.273 is missing
  # 
  # # Since epi prob is NA when a simulation does not reach that case count we will replace those NA with 0
  # # The NAs generated because we didn't run a simulation, we'll interpolate the missing epi prob
  # all_epi_data_interpolate =  all_epi_data %>%
  #   mutate(rnot_round1 = as.character(rnot_round1), # values have to be characters for join
  #          daily_import_round3 = as.character(daily_import_round3),
  #          detected = as.character(detected),
  #          prob_epidemic_interp = prob_epidemic,
  #          epi_origin = ifelse(!is.na(prob_epidemic_interp), "simulation", "fill")) %>%
  #   group_by(rnot_round1, daily_import_round3) %>%
  #   fill(prob_epidemic_interp, .direction = "down") %>% # fill values with the number above, so from top down of group
  #   ungroup() %>% 
  #   
  #   # interpolate between R0 values 
  #   full_join(every_combo, by=c("rnot_round1", "daily_import_round3", "detected")) %>%
  #   mutate(epi_origin = ifelse(is.na(prob_epidemic_interp), "interpolate_r0", epi_origin)) %>%
  #   group_by(daily_import_round3, detected) %>%
  #   arrange(rnot_round1) %>%
  #   mutate(total_na = sum(!is.na(prob_epidemic_interp))) %>%
  #   filter(total_na > 0) %>% # can't interpolate between values where there is no data
  #   #mutate(prob_epidemic_interp = zoo::na.approx(prob_epidemic_interp) ) %>% # Interpolate values
  #   ungroup() %>%
  # 
  #   # interpolate between importation values
  #   full_join(every_combo, by=c("rnot_round1", "daily_import_round3", "detected")) %>%
  #   mutate(epi_origin = ifelse(is.na(prob_epidemic_interp), "interpolate_import", epi_origin)) %>%
  #   group_by(rnot_round1, detected) %>%
  #   arrange(daily_import_round3) %>%
  #   mutate(total_na = sum(!is.na(prob_epidemic_interp))) %>%
  #   filter(total_na > 0) %>% # can't interpolate between values where there is no data
  #   mutate(prob_epidemic_interp = zoo::na.approx(prob_epidemic_interp) ) %>% # Interpolate values
  #   ungroup()
  #   
  # write.csv(all_epi_data_interpolate, paste0(output_folder_path, "all_epi_by_r0_import_interp.csv"), row.names=F)
  
} # end function epi_prob_by_r0_import


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



