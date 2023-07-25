###############################
## Code for simulations
###############################

# Overdispersion parameter estimation
# Based on this paper: https://onlinelibrary.wiley.com/doi/10.1111/tbed.14655#pane-pcw-figures
fit_dispersion <- function(dispersion=0.16, rnot = 1.1){
  P <- 0.2 # For malaria from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3128496/
  Q <- 0.8 # and https://www.nature.com/articles/438293a
  xs <- 0:100 # X = 100
  cum_sum_vals <- cumsum(xs * dnbinom(xs, mu = rnot, size = dispersion) / rnot) # 1 - Q
  x_val <- xs[ which(cum_sum_vals<(1-Q))[length( which(cum_sum_vals<(1-Q)) )] ] # get upper limit
  
  return( abs(1-P - pnbinom(q = x_val, mu = rnot, size = dispersion)) )
} # end fit_dispersion

get_optimal_dispersion <- function(rnot){
  max_search = 10
  if(rnot==0){ # When 0 neg binom could return anything, so set as the max of interval below
    min_disp = max_search
  }else{
    fit_disp_param <- optimise(fit_dispersion, interval = c(0.0,max_search), rnot = rnot)
    min_disp = fit_disp_param$minimum
  } # end if
  return(min_disp)
} # end get_optimal_dispersion

###################################################
#' Default malaria parameters for sims
#' 
#' A function that gives a list of parameters for running a malaria simulation.
#' Defaults to parameters for malaria with no introductions.
#' @param base_r_not Effective reproduction number: value specifying outbreak R0, which automatically modulates transmission rate
#' @param base_det_prob Daily detecting rate: The daily probability of an infectious individual being detected
#' @param intro_rate Daily rate of new introductions according to poisson distribution
#' @param infBoxes Integer specifying number of infectious class boxes
#' @param incBoxes Integer specifying number of incubation class boxes
#' @param gen_time Generation time (days): Average length of time between consecutive exposures T_G = T_E + (T_I/2)
#' @param inf_period Infectious period (days)
#' @param exp_period Latent period (days)
#' @param incub_rate Daily probability of transitioning for each incubation box
#' @param recov_p Daily probability of recovery for each infectious box
#' @param base_dispersion Total dispersion parameter of negative binomial distribution
#' @param dispersion Daily dispersion parameter of negative binomial distribution
#' @param prop_p Daily Poisson rate of transmission for each infectious individual
#' @param e_thresh Maximum number of cumulative infections before sim ends
#' @param prob_symp Probability that cases are symptomatic (haven't tested anything <1)
#' @param dis_prob_symp Probability of discovery for symptomatic individuals
#' @param dis_prob_asymp Probability of discovery for asymptomatic individuals
#' @return A list of length num_reps, where each component is a single call to run_malaria_sim.
malaria_params_fn = function(base_r_not           = 1.5,                       ## county specific 
                             base_det_prob        = 0.25,                      ## symptomatic_rate*self_report_rate = 0.25*0.154
                             intro_rate           = 0.01,                      ## county specific (prob of import * 2044)/365
                             
                             gen_time             = 51,                        ## https://www.mdpi.com/2076-2607/8/7/984
                             inf_period           = 34,                        ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4323116/
                             infBoxes             = 22,                        ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4323116/
                             exp_period           = gen_time - 0.5*inf_period, ## Roberts & Heesterbeek 2007 
                             incBoxes             = 28,                        ## (7.9+10+10) very rough minimum
                             incub_rate           = incBoxes/exp_period, 	     ## nu 	  
                             recov_p              = infBoxes/inf_period,       ## delta

                             base_dispersion      = get_optimal_dispersion(base_r_not), ## see above
                             dispersion           = base_dispersion/inf_period,  
                             prop_p               = base_dispersion/(base_r_not+base_dispersion), # beta
                             e_thresh             = 1000,                       # was 2000 for covid and 500 for zika
                             prob_symp            = 1,
                             dis_prob_symp        = base_det_prob/inf_period,
                             dis_prob_asymp       = 0.0
){ return(as.list(environment())) }

#' Run multiple malaria simulations
#' 
#' A function to run multiple malaria simulations with the same parameters
#' @return A list of length num_reps, where each component is a single call to run_malaria_sim.
run_n_malaria_sims <- function(num_reps, parms) {
  #plyr::rlply(.n = num_reps, .expr = run_malaria_sim(...) )
  map(seq_len(num_reps), ~run_malaria_sim(parms))
}

#' Update state movement
#'
#' A function that calculates how many individuals transition out of a state
#' @param n An integer number of draws.
#' @param prob A probability of transitioning.
#' @return Number of individuals leaving state.
#' @examples
#' update_state(100, 0.3)
update_state <- function(n, prob){
  ## Draws n uniform random numbers
  ## Returns a sum of all numbers less than the probability
  sum(runif(n = n) < prob)
}

#' Update multiple state's movement
#'
#' A vectorized function of update_state()
#' @param n An integer vector for the number of draws for each state.
#' @param prob A probability of transitioning.
#' @return Number of individuals leaving each state.
#' @examples
#' update_states(c(50,100,20,30), 0.3)
update_states <- Vectorize(update_state, vectorize.args = c("n"))


#' Track movement through infectious classes
#'
#' A function that tracks movement between undetected and detected classes of individuals.
#' Progresses individuals according to one step (day), and accounts for the fact that individuals
#' can recover and be reported in the same time period
#'
#' @param UI Integer vector length=number of boxes
#' @param recov_prob Probability of recovery
#' @param disc_prob Probability of being reported
#' @return Outputs a list UI recovering, UI reporting, and amount to add to each DI box
#' @examples
#' track_infectious(c(10,50,20), 0.3, 0.5)
track_infectious <- function(UI, recov_prob, disc_prob){
  UI_progressing = vector("numeric", length(UI))  ## The number of UI moving to next UI box
  UI_detecting = vector("numeric", length(UI))    ## The number of UI that were detected
  DI_adding = vector("numeric", length(UI)+1)     ## The number of UI that are detected by stage of recovery
  
  for(box in 1:length(UI)){
    UI_recovery_draws = runif(UI[box])
    UI_detected_draws = runif(UI[box])
    
    ## Get the total numbers recovering, detected, and (detected and recovered)
    UI_recovering = sum(UI_recovery_draws < recov_prob)
    UI_detected = sum(UI_detected_draws < disc_prob)
    UI_detected_and_recovered = sum( (UI_detected_draws < disc_prob) & (UI_recovery_draws < recov_prob) )
    
    ## Keep track of total leaving the current UI box into another UI box
    UI_progressing[box] = UI_recovering - UI_detected_and_recovered
    
    ## Neeed to keep track of those leaving due to detection
    UI_detecting[box] = UI_detected
    
    ## Now update those entering the DI class
    ## Those that didn't recover, but were detected, are added to current DI box
    DI_adding[box] = DI_adding[box] + UI_detected - UI_detected_and_recovered
    
    ## Those that recovered and detected are added to next DI class
    DI_adding[box+1] = DI_adding[box+1] + UI_detected_and_recovered
  }
  
  return(list(UI_progressing=UI_progressing, UI_detecting=UI_detecting, DI_adding=DI_adding))
}


#' Runs a single malaria simulation
#' @return Returns a data frame that contains the daily evolution of the model compartments

run_malaria_sim <- function(params) {
  with(params,{
    ## The 4 infected classes all have that number of boxes
    UI_Symp = vector("numeric", length=infBoxes)
    UI_Asymp = vector("numeric", length=infBoxes)
    DI_Symp = vector("numeric", length=infBoxes)
    DI_Asymp = vector("numeric", length=infBoxes)
    
    ## Classes for storing the introduced infectious
    UI_Intro_Symp = vector("numeric", length=infBoxes)
    UI_Intro_Asymp = vector("numeric", length=infBoxes)
    DI_Intro_Symp = vector("numeric", length=infBoxes)
    DI_Intro_Asymp = vector("numeric", length=infBoxes)
    
    ## Initialize the epidemic with one undetected infectious
    UI_Intro_Symp[1] = 1
    
    ## UI and DI keep track of the total number of undetected and detected
    ## cases at any given time
    UI = sum(UI_Symp, UI_Asymp, UI_Intro_Symp, UI_Intro_Asymp)
    DI = sum(DI_Symp, DI_Asymp, DI_Intro_Symp, DI_Intro_Asymp)
    
    ## Only those that are introduced
    UI_Intro = sum(UI_Intro_Symp, UI_Intro_Asymp)
    DI_Intro = sum(DI_Intro_Symp, DI_Intro_Asymp)
    
    ## I keeps track of total infectious
    I = UI + DI
    
    ## Only introduced ones
    I_Intro = UI_Intro + DI_Intro
    
    cumD = 0  ## Cum detected
    cumD_Intro = 0  ## Cum detected of introduced
    
    ## Keeps track of all exposed classes
    incubationInfecteds = vector("numeric", length = incBoxes)
    
    ## Cumulative incidence
    cumI = I
    cumI_Intro = I_Intro
    
    ## Data.frame for holding time series
    time_record = data.frame(matrix(data = 0, ncol = 11))
    colnames(time_record) = c("New_Exposed", "New_Infection", "Total_Infections", "Cumulative_Infections",
                              "New_Intro", "Total_Intro_Infections", "Cumulative_Intro_Infections",
                              "New_Detections", "Cum_Detections", "New_Intro_Detections", "Cum_Intro_Detections")
    
    ## Initialize appropriately
    time_record$New_Infection = time_record$Total_Infections = time_record$Cumulative_Infections = I
    
    while  (( (cumI-cumI_Intro) < e_thresh) & ((I+sum(incubationInfecteds)) > 0)) {
      #while Number of infected is below epidemic threshold and more than 0 infected
      # or while number of infecteds is above 0 and the number of detected is below threshold
      
      ########################### First Introduction Rate
      ## intro_rate is rate of introductions, drawn from poisson distribution daily
      intro_count = rpois(n = 1, lambda = intro_rate)
      
      ## get total introductions that are symptomatic
      intro_symp = update_state(n=intro_count, prob = prob_symp)
      intro_asymp = intro_count-intro_symp
      
      ############################## INCUBATION
      ## Keep track of transitions out of each incubation box
      leavingInc = update_states(n = incubationInfecteds, prob = incub_rate)
      
      ## DETERMING ASYMP/SYMP of those leaving last incubation box (becoming infectious)
      newUI_Symp = update_state(n = leavingInc[incBoxes], prob = prob_symp)
      newUI_Asymp = leavingInc[incBoxes] - newUI_Symp
      
      ############################## INFECTION
      #Calculate the number of new exposed from all infectious
      ## Add in second line if decide that symp/asymp should have different transmission rates
      ## potentially detected/undetected as well.
      newExposed = sum(rnbinom(n = I, prob = prop_p, size = dispersion))
      
      ############################## INFECTIOUS CLASSES MOVEMENT
      # Infectious classes movement Undetected infectious
      ## Can progress in UI boxes, or get transferred
      ## to DI boxes
      infectious_update_Symp = track_infectious(UI = UI_Symp, recov_prob = recov_p, disc_prob = dis_prob_symp)
      infectious_update_Asymp = track_infectious(UI = UI_Asymp, recov_prob = recov_p, disc_prob = dis_prob_asymp)
      
      ## Update introduction infectious classes
      infectious_intro_update_Symp = track_infectious(UI = UI_Intro_Symp, recov_prob = recov_p, disc_prob = dis_prob_symp)
      infectious_intro_update_Asymp = track_infectious(UI = UI_Intro_Asymp, recov_prob = recov_p, disc_prob = dis_prob_asymp)
      
      ## DI movement (can only progress towards recovery)
      DI_symp_leaving <- update_states(n=DI_Symp, prob=recov_p)
      DI_asymp_leaving <- update_states(n=DI_Asymp, prob=recov_p)
      
      DI_intro_symp_leaving <- update_states(n=DI_Intro_Symp, prob=recov_p)
      DI_intro_asymp_leaving <- update_states(n=DI_Intro_Asymp, prob=recov_p)
      
      #########################################################
      # UPDATING
      # Incubation classes
      for(i in 1:incBoxes){
        if(i==1){
          incubationInfecteds[i] = incubationInfecteds[i] + newExposed - leavingInc[i]
        } else {
          incubationInfecteds[i] = incubationInfecteds[i] + leavingInc[i-1] - leavingInc[i]
        }
      }
      
      # UI_leaving=UI_leaving, DI_adding=DI_adding
      # Infectious classes updating
      for(i in 1:infBoxes){
        if(i==1){
          ## Update the UI first box
          UI_Symp[i] = UI_Symp[i] +  newUI_Symp - infectious_update_Symp$UI_progressing[i] -  infectious_update_Symp$UI_detecting[i]
          UI_Asymp[i] = UI_Asymp[i] + newUI_Asymp - infectious_update_Asymp$UI_progressing[i] -  infectious_update_Asymp$UI_detecting[i]
          
        } else {
          ## Update the other UI boxes
          UI_Symp[i] = UI_Symp[i] + infectious_update_Symp$UI_progressing[i-1] - infectious_update_Symp$UI_progressing[i] -  infectious_update_Symp$UI_detecting[i]
          UI_Asymp[i] = UI_Asymp[i] + infectious_update_Asymp$UI_progressing[i-1] - infectious_update_Asymp$UI_progressing[i]-  infectious_update_Asymp$UI_detecting[i]
          
        }
        
        ## Update the DI boxes (all the same process, so not in conditionals)
        DI_Symp[i] = DI_Symp[i] + infectious_update_Symp$DI_adding[i] - DI_symp_leaving[i]
        DI_Asymp[i] = DI_Asymp[i] + infectious_update_Asymp$DI_adding[i] - DI_asymp_leaving[i]
      }
      
      for(i in 1:infBoxes){
        if(i==1){
          ## Update the UI first box
          UI_Intro_Symp[i] = UI_Intro_Symp[i] + intro_symp - infectious_intro_update_Symp$UI_progressing[i] -  infectious_intro_update_Symp$UI_detecting[i]
          UI_Intro_Asymp[i] = UI_Intro_Asymp[i] + intro_asymp +  infectious_intro_update_Asymp$UI_progressing[i] -  infectious_intro_update_Asymp$UI_detecting[i]
          
        } else {
          ## Update the other UI boxes
          UI_Intro_Symp[i] = UI_Intro_Symp[i] + infectious_intro_update_Symp$UI_progressing[i-1] - infectious_intro_update_Symp$UI_progressing[i] -  infectious_intro_update_Symp$UI_detecting[i]
          UI_Intro_Asymp[i] = UI_Intro_Asymp[i] + infectious_intro_update_Asymp$UI_progressing[i-1] - infectious_intro_update_Asymp$UI_progressing[i]-  infectious_intro_update_Asymp$UI_detecting[i]
          
        }
        
        ## Update the DI boxes (all the same process, so not in conditionals)
        DI_Intro_Symp[i] = DI_Intro_Symp[i] + infectious_intro_update_Symp$DI_adding[i] - DI_intro_symp_leaving[i]
        DI_Intro_Asymp[i] = DI_Intro_Asymp[i] + infectious_intro_update_Asymp$DI_adding[i] - DI_intro_asymp_leaving[i]
      }
      
      
      # Newly Infectious
      newInf = leavingInc[incBoxes] + intro_symp + intro_asymp
      newInf_Intro = intro_symp + intro_asymp
      
      # Update cumulative infectious
      cumI = cumI + newInf
      cumI_Intro = cumI_Intro + newInf_Intro
      
      # Newly Detected = symptomatic and asymptomatic cases discovered
      newlyDisc = sum(infectious_update_Symp$DI_adding, infectious_update_Asymp$DI_adding, infectious_intro_update_Symp$DI_adding, infectious_intro_update_Asymp$DI_adding)
      newlyDisc_Intro = sum(infectious_intro_update_Symp$DI_adding, infectious_intro_update_Asymp$DI_adding)
      
      #Cumulatve Detected = those already detected + newly detected
      cumD = cumD + newlyDisc
      cumD_Intro = cumD_Intro + newlyDisc_Intro
      
      UI = sum(UI_Symp, UI_Asymp, UI_Intro_Symp, UI_Intro_Asymp)
      DI = sum(DI_Symp, DI_Asymp, DI_Intro_Symp, DI_Intro_Asymp)
      
      ## Only those that are introduced
      UI_Intro = sum(UI_Intro_Symp, UI_Intro_Asymp)
      DI_Intro = sum(DI_Intro_Symp, DI_Intro_Asymp)
      
      #Detected = Detected + Recovered Individuals + Newly Detected - Previously UI that have been detected and Recovered
      I = UI + DI #Undiscovered and Discovered Infected
      I_Intro = UI_Intro + DI_Intro
      
      #adding time step data
      time_record <- rbind(time_record, c(newExposed, newInf, I, cumI, newInf_Intro, I_Intro, cumI_Intro, newlyDisc, cumD, newlyDisc_Intro, cumD_Intro))
    }
    time <- data.frame(time=seq(1:nrow(time_record)))
    time_record <- cbind(time, time_record)
    
    return(time_record)
  })
} # end function run_malaria_sim

# path for desired .rda
get_save_path <- function(base_r_not, 
                          intro_rate,
                          num_reps,
                          path,
                          summary = FALSE){
  if(summary){
    paste0(path, "county-summary_", base_r_not, "_", intro_rate, "_", num_reps, "_", Sys.Date(), ".rda")
  }else{
    paste0(path, "sim_", base_r_not, "_", intro_rate, "_", num_reps, "_", Sys.Date(), ".rda")  
  } # end if else
} # end function get_save_path

save_malaria_runs <- function(num_reps,
                              refresh=FALSE, # when TRUE will re-write files
                              ...) {
  parms <- malaria_params_fn() # get the parameters defined above or those passed
  
  saved_file_path <- get_save_path(base_r_not, intro_rate, num_reps, path) # pass params that vary to change file names
  if(!file.exists(saved_file_path) | refresh){
    print("Running the simulation")
    print(parms)
    sims <- run_n_malaria_sims(num_reps, parms)  
    save(sims, file = saved_file_path)
    print("Simulation completed")
    return("Refreshed")
  } else{
    return("Already run before")
  } # end if else
} # end function save_malaria_runs

###########################################################
# Convert simulation data to epidemic probability table
get_epi_prob = function(data_path){
  
  load(data_path)
  csv_path=str_replace(string = data_path, pattern = "rda", replacement = "csv")
  epi_path=str_replace(string = csv_path, pattern = "sim", replacement = "epi_prob_data")
  
  if(file.exists(epi_path)){
    epi_probs = read_csv(epi_path)
  }else{
    epi_probs <- get_epidemic_prob_by_d(trials = sims, 
                                        prev_threshold = 50,
                                        cum_threshold = 1000, # should match e_thresh
                                        max_detect = 100) # the max number of cases to get epi_prob for
    epi_probs$prob_epidemic[is.na(epi_probs$prob_epidemic)]=1
    
    ## write to processed_data folder to use for plotting
    write.csv(epi_probs, epi_path, row.names = FALSE)
  }
  return("Epi Probs generated")
} # end function get_epi_data
















