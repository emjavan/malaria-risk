###################################################################################################
# Get all the unique combinations of R0 and importation prob
# R0 rounded to 1 decimal and import to 3
# Import kept at 3 decimals bc 1 import in a year is 1/365~0.003 daily probability of importation
# Run from inside code folder where script is kept
###################################################################################################

# Importation risk from random forest - newest version
impProbByCounty_new = read_csv("../input_data/malaria_probability_next_import_comes_from_county_2023-08-14.csv") %>%
  mutate(daily_import = (ProbNextOcc*2044 - 1)/365, # shifting all values down by 1 bc the sim starts with 1 import
         daily_import_round3 = round(daily_import, 3),
         fips=as.character(fips),
         fips=ifelse(nchar(fips)==5, fips, paste0(0, fips)) )
impProbByCounty_new$daily_import_round3[impProbByCounty_new$daily_import_round3<0]=0

uniq_imports = sort(unique(impProbByCounty_new$daily_import_round3))
length(uniq_imports) # 50 unique daily import probs
range(uniq_imports) # 0 to 0.178 # these are the new values from 2023-07-31

# County specific R0
load("../input_data/county-single-rnot-estimates_2023-08-09.rda")
county_single_rnot_samples_new = county_single_rnot_samples %>%
  mutate(rnot_round1 = round(rnot, 1),
         fips=as.character(fips),
         fips=ifelse(nchar(fips)==5, fips, paste0(0, fips)) ) %>%
  left_join(impProbByCounty_new, by="fips")
county_single_rnot_samples_new[is.na(county_single_rnot_samples_new)]=0
uniq_r0 = sort(unique(county_single_rnot_samples_new$rnot_round1))
length(uniq_r0) # 52
range(uniq_r0) # 0.0 to 5.4
write.csv(county_single_rnot_samples_new,
          "../processed_data/all_rnot_import_per_county_2023-08-09.csv", row.names = F)

rnot_round1=uniq_r0; daily_import_round3=uniq_imports
all_r0_import = expand_grid(rnot_round1, daily_import_round3) # 52 * 50 = 2600

# Generate all possible commands needed for Frontera
all_r0_import_command = all_r0_import %>%
  mutate(command = 
           paste0("Rscript --no-save run-malaria-sims-only_parallel.R ", rnot_round1, " ", daily_import_round3)) %>%
  select(command)
write.table(all_r0_import_command, "../code/frontera_parallel_code/commands_malaria.txt", 
            row.names=F, col.names=F, quote=F)







