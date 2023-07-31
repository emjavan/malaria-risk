#!/bin/bash

# Generate the commands to run a test in parallel on Frontera
# First input (outer loop) is base_r_not, second (inner) is the intro_rate
# Not sure why I get some weird rounding issue, but including the 0.01 after the last value in sequence to ensure inclusion
for i in $(seq 0.0 0.1 1.51); do
   for j in $(seq 0.0 0.05 0.31); do
      echo "Rscript --no-save run-malaria-sims-only_parallel.R" $i $j >> commands_malaria.txt
   done
done

# get the number of simulation to run
# good to confirm this number matches your expectation as >> just add to file and will not over-write existing one
# add to launcher script as the number of tasks, then estimate total computers and how much time you need
wc -l commands_malaria.txt
