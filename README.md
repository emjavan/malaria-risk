# Estimate the risk/probability of an undetected Malaria epidemic in each county of continental US

The stochastic SEIR model begins each simulation with one imported infectious individual at day 1 of their infectious period. The simulation terminates when there are no more infectious individuals or the total autochthonous infections reach 1,000. The model assumes an unlimited pool of susceptible individuals, as is often the case in a naive population, thus will never self-limit. Each day of the simulation the number of new infections caused by an infectious individual is chosen by a random draw from a negative binomial distribution with an overdispersion parameter (lambda) dependent on the R<sub>0</sub>. As R<sub>0</sub> increases, lambda decays  Infectious individuals have an equal daily probability of being detected and removed. We assume 25% of all infectious individuals in a simulation will be detected in accordance with the proportion of *Plasmodium vivax* infections that become symptomatic. We do not explicitly model asymptomatic infections, although the mechanics to do so are in the model. 

Model code and parameters can be found in `sim-malaria-outbreaks.R`, but parameters are listed below in Table 1 as well.



Table 1. 
| Parameter   | Value  | Source   |
| :---        | ----:  |     ---: |
| Generation time (days) T<sub>G</sub>    | 51          | https://www.mdpi.com/2076-2607/8/7/984   |
| Latent period (days) T<sub>E</sub>      | 34          | T<sub>G</sub> - 0.5\*T<sub>I</sub>, Roberts & Heesterbeek 2007      |
| Infectious period (days) T<sub>I</sub>  | 34          | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4323116/   |
| Min day exposed *e*                     | 28          | (7.9 + 10 + 10)      |
| Min day infectious *n*                  | 22          | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4323116/   |
| R<sub>0</sub>                           | [0, 5.4]    | Estimated     |
| Daily importation                       | [0, 0.178]  | (ProbNextOcc\*2044 - 1)/365   |

## Ingest R0 and probability of next importation from `input_data`, then write commands file
1. Run `clean_r0_import_write_commands.R` to get the commands file to run all parameter combinations in parallel on Frontera. `generate_commands.sh` is simply to generate an example for testing if you do not want to queue entire job. 

## Run simulations on TACC Frontera
2. `sbatch launch_malaria.sh` from `code/frontera_parallel_code/`. You can check on the number of finished jobs by running `grep "sec" *.o | wc -l` this will count the number of times "sec" appears in the standard output file. You could also `ls sim_* | wc -l`inside `full_run_processed_data`. Job will run full 48hr.

As a brief hierarchy, `launch_malaria.sh` < `commands_malaria.txt` < `run-malaria-sims-only_parallel.R` < `sim-malaria-outbreaks.R`

## Ensure all expected files produced
3. Check all sims finished with `produced_file_check.R`. If for some reason file date aren't consistent, change with `change_file_name.sh`.

## Re-run files that cannot finish in time limit
4. Some sims with a low R0 but high importation will not finish in the 48hr max time limit, so they are run as 1K x10 sets
run `batch_sims_wont_finish.sh` to re-make `commands_malaria_missing.txt`. After updating `run-malaria-sims-only_parallel.R` to 1K sims and to take a third command line input, run `launch_malaria_missing.sh`, ~8.5hr to 40 tasks finish. Then post process the 10 1K sets to make one 10K set with `XXX`. Future iterations of code should put a limit on number of simulation days each runs for.

## Generate summary stats and figures
5. All functions run from `run-malaria-plotting-code.R` this file calls `frontera_parallel_code/outbreak_analysis_fxns.R` and `malaria-plotting-code.R`.

As the first job is running to collate final results (mean epi prob across 1K R0 per county) you can check on it with `frontera_parallel_code/check_on_running_code.sh`. 

Making Figure 4 doesn't require any epi prob data, just the raw sims, so this could be run while waiting for the missing sims.

















