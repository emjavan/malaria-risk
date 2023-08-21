rm commands_malaria_missing.txt

for i in {1..10}; do
   echo "Rscript --no-save run-malaria-sims-only_parallel.R 0.2 0.178 ../../processed_data/full_run_processed_data/set${i}/
Rscript --no-save run-malaria-sims-only_parallel.R 0.3 0.178 ../../processed_data/full_run_processed_data/set${i}/
Rscript --no-save run-malaria-sims-only_parallel.R 0.3 0.161 ../../processed_data/full_run_processed_data/set${i}/
Rscript --no-save run-malaria-sims-only_parallel.R 0.4 0.178 ../../processed_data/full_run_processed_data/set${i}/" \
   >> commands_malaria_missing.txt
done
