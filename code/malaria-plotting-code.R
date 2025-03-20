###############################
## Code for plotting
###############################



plot_continental_US_epi_and_trigger=function(run_location="../"){
 
  
  # Read in county shapefile data if it doesn't exist, otherwise create rda
  rda_file = paste0(run_location, "processed_data/us_county_geometries_2021.rda")
  if(file.exists(rda_file) ){
    load(rda_file) # us_county_geo
    print("US county data loaded")
  } else{
    us_county_geo = tidycensus::get_acs(geography="county", variables=c("B01001_001"), 
                                        geometry=TRUE, year=2021, shift_geo = TRUE) %>%
      separate(NAME, into = c("county", "state"), sep = ", ") %>%
      rename(fips=GEOID)
    
    save(us_county_geo, file = "processed_data/us_county_geometries.rda")
    print("us county data saved")
  }
  
  # Summarized epidemic probability across all 1K R0 values for a county
  # Some counties have missing data bc sims didn't run => will change to be look-up later
  epi_file_path = 
    paste0(run_location, "processed_data/full_run_processed_data/clean_epi_files/county_real_epi_prob_and_trigger.csv")
  
  county_data = read_csv(epi_file_path) %>%
    left_join(us_county_geo %>% select(fips, geometry), by="fips") %>% 
    mutate(case_when_epi_over_50_cut_off = ifelse(case_when_epi_over_50>20, 20, case_when_epi_over_50))
  
  require(sf)
  p1 = ggplot(county_data)+
    geom_sf(mapping=aes(geometry=geometry, fill=epi_prob), size = 0.05, color="black")+
    #geom_polygon(aes(group = group, fill = epi_prob*100), color = "black", size = 0.1) +
    scale_fill_gradient(low = "gainsboro", high = "darkred", name = "Epi Prob", na.value = "white",
                        #breaks=c(0.0, 0.25, 0.50, 0.75, 1.0), 
                        #labels=c(0.0, 0.25, 0.50, 0.75, 1.0), limits=c(0,1.0) 
                        )+
    scale_x_continuous("", breaks = NULL) + scale_y_continuous("", breaks = NULL)+
    theme_void()+
    theme(
      legend.position = "right",
      legend.title=element_text(size=10),
      legend.text=element_text(size=8),
      legend.text.align = 1,
      #legend.key.size = unit(0.5, "lines"),
      #plot.margin=unit(c(0, 0, 0, 0),"cm"),
      panel.background  = element_rect(color = "white", fill = "white"),
      legend.background = element_rect(color = "white", fill = "white"))+
    guides(fill = guide_colourbar(ticks.colour='black', ticks.linewidth = 0.25, # barheight = unit(2, "in"), 
                                  border=element_line(color='black')))
  
  p2 = ggplot(county_data )+
    geom_sf(mapping=aes(geometry=geometry, fill=case_when_epi_over_50_cut_off), size = 0.05, color="black")+
    #geom_polygon(aes(group = group, fill = epi_prob*100), color = "black", size = 0.1) +
    scale_fill_gradient(high = "gainsboro", low = "#2c7fb8", name = "Case Prob>0.5", na.value = "white",
                        #breaks=c(0.0, 0.25, 0.50, 0.75, 1.0), 
                        #labels=c(0.0, 0.25, 0.50, 0.75, 1.0), limits=c(0,1.0) 
    )+
    scale_x_continuous("", breaks = NULL) + scale_y_continuous("", breaks = NULL)+
    theme_void()+
    theme(
      legend.position = "right",
      legend.title=element_text(size=10),
      legend.text=element_text(size=8),
      legend.text.align = 1,
      #legend.key.size = unit(0.5, "lines"),
      #plot.margin=unit(c(0, 0, 0, 0),"cm"),
      panel.background  = element_rect(color = "white", fill = "white"),
      legend.background = element_rect(color = "white", fill = "white"))+
    guides(fill = guide_colourbar(ticks.colour='black', ticks.linewidth = 0.25, # barheight = unit(2, "in"), 
                                  border=element_line(color='black')))
  
  
  plot_list = list(p1, p2)
  
  return(plot_list)
}

#' Plot heat map of epi prob by R0 and importation across 0-7 cases
#' 
#' @param run_location ="../" assumes you're running from inside code, change to "" if in malaria-risk
#' 
epi_r0_import_plot = function(run_location="../"){
  # file created by epi_prob_by_r0_import in outbreak_analysis_fxns.R
  folder_path = paste0(run_location, "processed_data/full_run_processed_data/clean_epi_files/")
  epi_data = read_csv(paste0(folder_path, "all_epi_by_r0_import_og.csv"))
  
  heat_map_og = ggplot(epi_data, aes(x=rnot_round1, y=daily_import_round3))+
    geom_tile(aes(fill=prob_epidemic))+
    #geom_raster(aes(fill=prob_epidemic), interpolate=TRUE)+ # continuous interpolation, use geom_tile for discrete
    scale_fill_gradient2(low="grey40", mid="grey20", high="black", na.value = "orange",
                         midpoint=0.5)+ # , limits=range(epi_data$prob_epidemic)
    facet_wrap(~detected, nrow = 2)+
    labs(x="R0", y="Daily Importation Risk", fill="Epi Prob")+
    theme_classic(base_size = 10)
  
  # heat_map_interp = ggplot(epi_data, aes(x=rnot_round1, y=daily_import_round3))+
  #   geom_tile(aes(fill=prob_epidemic_interp))+
  #   #geom_raster(aes(fill=prob_epidemic), interpolate=TRUE)+ # continuous interpolation, use geom_tile for discrete
  #   scale_fill_gradient2(low="grey40", mid="grey20", high="black", na.value = "orange",
  #                        midpoint=0.5)+ # , limits=range(epi_data$prob_epidemic)
  #   facet_wrap(~detected, nrow = 2)+
  #   labs(x="R0", y="Daily Importation Risk", fill="Epi Prob")+
  #   theme_classic(base_size = 10)
  
  return(heat_map_og)
} # end epi_r0_import_plot



plot_case_by_prob = function(run_location="../"){
  
  # Open needed files for plot
  # Case detect date by the week of detection in Florida reports
  sara_detect_dates = read_csv(paste0(run_location, "input_data/sarasota_county_cases_through_time_2023-07-29.csv"))
  # County mean epi_prob by case
  sarasota = 
    read_csv(paste0(run_location, "processed_data/full_run_processed_data/clean_epi_files/12115_all_epi_prob_by_case.csv")) %>%
    select(fips, cases, epi_prob) %>%
    filter(cases <= 20)
  prince_g = 
    read_csv(paste0(run_location, "processed_data/full_run_processed_data/clean_epi_files/24033_all_epi_prob_by_case.csv")) %>%
    select(fips, cases, epi_prob) %>%
    filter(cases <= 20)
  
  # # Summarized epi_prob just to check what case epi>50 for plot
  # epi_prob = read_csv(paste0(run_location,
  #                            "processed_data/full_run_processed_data/clean_epi_files/county_real_epi_prob_and_trigger.csv")) %>%
  #   filter(fips %in% c("12115", "24033")) %>%
  #   select(fips, epi_prob, case_when_epi_over_50)
  # Sarasota = 8 and Prince George's is NA, so we'll plot epi risk for up to 20 cases
  
  # R0 and Importation data for counties of interest
  county_of_interest = read_csv(paste0(run_location, 
                                        "processed_data/all_rnot_import_per_county_2023-08-09.csv")) %>%
    filter(fips %in% c("12115", "24033")) %>%
    mutate(county_name = ifelse(fips=="12115", "Sarasota", "Prince George's") )
  
  p1=ggplot(county_of_interest, aes(x=rnot_round1, color=county_name, fill=county_name))+
    #geom_density(alpha=0.3)+
    geom_histogram(alpha=0.3)+
    #facet_grid(~fips)+
    labs(x="R0", color="", fill="")+
    theme_bw()+
    theme(legend.position = c(0.9, 0.5), legend.justification="right")
  
  boundary_path = paste0(run_location, "processed_data/full_run_processed_data/epi_files_for_county_plot/")
  r0_from_boundaries = list.files(boundary_path) %>%
    as_tibble() %>%
    rename(filename = value) %>%
    separate(filename, into=c(NA, NA, NA, "rnot_round1", "import", NA, NA), remove=F, sep="_" ) %>%
    mutate(across(everything(), as.character ))
  
  sum_county_of_interest = county_of_interest %>%
    # county_single_rnot_summary %>%
    # filter(fips %in% c("12115", "24033"))
    mutate(rnot_round1 = round(rnot, 1)) %>%
    group_by(fips) %>%
    summarise(import = mean(daily_import_round3),
              #mean_rnot = round(mean(rnot_round1), 1) ,
              lwr_qnt = quantile(rnot_round1, probs=c(0.025)),
              upr_qnt = quantile(rnot_round1, probs=c(0.975)),
              #min_rnot = min(rnot_round1),
              #max_rnot = max(rnot_round1) 
              ) %>%
    ungroup() %>%
    gather(key=rnot_type, value=rnot_round1, -fips, -import) %>%
    mutate(across(everything(), as.character)) %>%
    left_join(r0_from_boundaries, by=c("rnot_round1", "import")) %>%
    mutate(file_contents = map(filename, ~read_csv(file.path(boundary_path, .)) %>%
                                 select(detected, prob_epidemic) %>%
                                 filter(detected <= 20) ) ) %>%
    unnest(cols = c(file_contents))
  
  
  

}







