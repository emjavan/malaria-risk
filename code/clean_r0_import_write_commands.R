

# Importation risk from maximum entropy model 
load("input_data/impProbByCounty_2023-07-14.RData")

impProbByCounty_new = impProbByCounty %>%
  mutate(daily_import = (ProbNextOcc*2044)/365,
         daily_import_round3 = round(daily_import, 3))

uniq_imports = unique(impProbByCounty_new$daily_import_round3)
length(uniq_imports) # 52 unique daily import probs
range(uniq_imports) # 0 to 0.063

# County specific R0
load("input_data/county-single-rnot-estimates_2023-07-13.rda")

county_single_rnot_samples_new = county_single_rnot_samples %>%
  mutate(rnot_round2 = round(rnot, 2)  ) %>%
  left_join(impProbByCounty_new %>% rename(fips=FIPS), by="fips")
uniq_r0 = unique(county_single_rnot_samples_new$rnot_round2)
length(uniq_r0) # 251 unique daily import probs
range(uniq_r0) # 0.00 to 3.65

# sorting is just for ease of looking at final commands txt file
all_r0_import = expand.grid(sort(uniq_r0), sort(uniq_imports) ) # 251 * 52 = 13052

distinct_r0_import = county_single_rnot_samples_new %>%
  select(fips, rnot_round2, daily_import_round3) %>%
  group_by(fips) %>%
  distinct() %>%
  ungroup() # 124007 > 13052, so it's best to run all unique combination and look up from file to plot map

all_r0_import_command = all_r0_import %>%
  rename(R0 = Var1, Import=Var2) %>%
  mutate(command = paste0("Rscript --no-save run-malaria-sims-only_parallel.R ", R0, " ", Import ) ) %>%
  select(command)
write.table(all_r0_import_command, "code/frontera_parallel_code/commands_malaria.txt", 
            row.names=F, col.names=F, quote=F)









