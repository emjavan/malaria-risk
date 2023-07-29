# Get all the unique combinations of R0 and importation prob
# R0 rounded to 2 decimal and import to 3

# Importation risk from maximum entropy model
impProbByCounty_new = read_csv("input_data/malaria_probability_next_import_comes_from_county_2023-07-22.csv") %>%
  mutate(daily_import = (ProbNextOcc*2044 - 1)/365, # shifting all values down by 1 bc the sim starts with 1 import
         daily_import_round3 = round(daily_import, 3),
         fips=as.character(fips),
         fips=ifelse(nchar(fips)==5, fips, paste0(0, fips)) )
impProbByCounty_new$daily_import_round3[impProbByCounty_new$daily_import_round3<0]=0

uniq_imports = unique(impProbByCounty_new$daily_import_round3)
length(uniq_imports) # 65 unique daily import probs
range(uniq_imports) # 0 to 0.273

# County specific R0
load("input_data/county-single-rnot-estimates_2023-07-13.rda")
county_single_rnot_samples_new = county_single_rnot_samples %>%
  mutate(rnot_round2 = round(rnot, 2),
         fips=as.character(fips),
         fips=ifelse(nchar(fips)==5, fips, paste0(0, fips)) ) %>%
  left_join(impProbByCounty_new, by="fips")
county_single_rnot_samples_new[is.na(county_single_rnot_samples_new)]=0
uniq_r0 = unique(county_single_rnot_samples_new$rnot_round2)
length(uniq_r0) # 251 unique R0
range(uniq_r0) # 0.00 to 3.65
write.csv(county_single_rnot_samples_new, 
          "processed_data/all_rnot_import_per_county_2023-07-22.csv", row.names = F)

# sorting is just for ease of looking at final commands txt file
all_r0_import = expand.grid(sort(uniq_r0), sort(uniq_imports) ) # 251 * 65 = 16315

distinct_r0_import = county_single_rnot_samples_new %>%
  select(fips, rnot_round2, daily_import_round3) %>%
  group_by(fips) %>%
  distinct() %>%
  ungroup() # 124007 >> 16315, so it's best to run all unique combination and look up from file to plot map

# Selecting 4 for a test run of the extremes
# Frontera may need the jobs batched out into sets of 1K commands
all_r0_import_command = all_r0_import %>%
  rename(R0 = Var1, Import=Var2) %>%
  #filter(R0 %in% c(min(uniq_r0), max(uniq_r0)) ) %>% # to test just the extremes of R0
  #filter(Import %in% c(min(uniq_imports), max(uniq_imports)) ) %>% # to test just the extremes of importation risk
  mutate(command = paste0("Rscript --no-save run-malaria-sims-only_parallel.R ", R0, " ", Import ) ) %>%
  select(command)
write.table(all_r0_import_command, "code/frontera_parallel_code/commands_malaria.txt", 
            row.names=F, col.names=F, quote=F)









