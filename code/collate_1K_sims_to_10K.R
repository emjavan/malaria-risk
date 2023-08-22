require(tidyverse)

rl="../"
sim_0.2_0.178 = sim_0.3_0.161 = sim_0.3_0.178 = sim_0.4_0.178 = list()
epi_prob_data_0.2_0.178 = epi_prob_data_0.3_0.161 = epi_prob_data_0.3_0.178 = epi_prob_data_0.4_0.178 = data.frame()
for(i in 1:10){
  path=paste0(rl, "processed_data/full_run_processed_data/set", i, "/") # path to open 1K files
  
  # Concatenate lists from 10 sim files
  load(paste0(path, "sim_0.2_0.178_1000_2023-08-15.rda")); sim_0.2_0.178_temp = sims;
  sim_0.2_0.178 = c(sim_0.2_0.178_temp, sim_0.2_0.178)
  load(paste0(path, "sim_0.3_0.161_1000_2023-08-15.rda")); sim_0.3_0.161_temp = sims;
  sim_0.3_0.161 = c(sim_0.3_0.161_temp, sim_0.3_0.161)
  load(paste0(path, "sim_0.3_0.178_1000_2023-08-15.rda")); sim_0.3_0.178_temp = sims;
  sim_0.3_0.178 = c(sim_0.3_0.178_temp, sim_0.3_0.178)
  load(paste0(path, "sim_0.4_0.178_1000_2023-08-15.rda")); sim_0.4_0.178_temp = sims;
  sim_0.4_0.178 = c(sim_0.4_0.178_temp, sim_0.4_0.178)
  
  # Row bind the 10 epi prob files together
  t1 = read_csv(paste0(path, "epi_prob_data_0.2_0.178_1000_2023-08-15.csv"))
  epi_prob_data_0.2_0.178 = rbind(epi_prob_data_0.2_0.178, t1)
  t2 = read_csv(paste0(path, "epi_prob_data_0.3_0.161_1000_2023-08-15.csv"))
  epi_prob_data_0.3_0.161 = rbind(epi_prob_data_0.3_0.161, t2)
  t3 = read_csv(paste0(path, "epi_prob_data_0.3_0.178_1000_2023-08-15.csv"))
  epi_prob_data_0.3_0.178 = rbind(epi_prob_data_0.3_0.178, t3)
  t4 = read_csv(paste0(path, "epi_prob_data_0.4_0.178_1000_2023-08-15.csv"))
  epi_prob_data_0.4_0.178 = rbind(epi_prob_data_0.4_0.178, t4)
} # end for i

out_path=paste0(rl, "processed_data/full_run_processed_data/")

save(sim_0.2_0.178, file=paste0(out_path, "sim_0.2_0.178_10000_2023-08-15.rda") )
save(sim_0.3_0.161, file=paste0(out_path, "sim_0.3_0.161_10000_2023-08-15.rda") )
save(sim_0.3_0.178, file=paste0(out_path, "sim_0.3_0.178_10000_2023-08-15.rda") )
save(sim_0.4_0.178, file=paste0(out_path, "sim_0.4_0.178_10000_2023-08-15.rda") )

epi_prob_data_0.2_0.178_10000 = epi_prob_data_0.2_0.178 %>%
  group_by(detected) %>%
  summarise(epi_numer = sum(epi_numer),
            epi_denom = sum(epi_denom),
            prob_epidemic = epi_numer/epi_denom) %>%
  ungroup()
write.csv(epi_prob_data_0.2_0.178_10000, paste0(out_path, "epi_prob_data_0.2_0.178_10000_2023-08-15.csv"), row.names=F)
  
epi_prob_data_0.3_0.161_10000 = epi_prob_data_0.3_0.161 %>%
  group_by(detected) %>%
  summarise(epi_numer = sum(epi_numer),
            epi_denom = sum(epi_denom),
            prob_epidemic = epi_numer/epi_denom) %>%
  ungroup()
write.csv(epi_prob_data_0.3_0.161_10000, paste0(out_path, "epi_prob_data_0.3_0.161_10000_2023-08-15.csv"), row.names=F)

epi_prob_data_0.3_0.178_10000 = epi_prob_data_0.3_0.178 %>%
  group_by(detected) %>%
  summarise(epi_numer = sum(epi_numer),
            epi_denom = sum(epi_denom),
            prob_epidemic = epi_numer/epi_denom) %>%
  ungroup()
write.csv(epi_prob_data_0.3_0.178_10000, paste0(out_path, "epi_prob_data_0.3_0.178_10000_2023-08-15.csv"), row.names=F)

epi_prob_data_0.4_0.178_10000 = epi_prob_data_0.4_0.178 %>%
  group_by(detected) %>%
  summarise(epi_numer = sum(epi_numer),
            epi_denom = sum(epi_denom),
            prob_epidemic = epi_numer/epi_denom) %>%
  ungroup()
write.csv(epi_prob_data_0.4_0.178_10000, paste0(out_path, "epi_prob_data_0.4_0.178_10000_2023-08-15.csv"), row.names=F)

