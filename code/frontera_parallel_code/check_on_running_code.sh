#!/bin/bash

# Code to check while running R code on Frontera if values you expect are being output
# substitute the names of files and expected columns as needed
# example below run inside processed_data/full_run_processed_data/clean_epi_files

# 3108 = header + 3107 counties in continental US
wc -l county_real_epi_prob_and_trigger.csv

# Find the column number with name case_when_epi_over_50 or epi_prob
col_num=$(awk -v RS=',' '/case_when_epi_over_50/{print NR; exit}' county_real_epi_prob_and_trigger.csv)

# fixed previous join issue and now stable 
awk -v RS=',' '/epi_prob/{print NR; exit}' county_real_epi_prob_and_trigger.csv
col2_num=$(head -1 county_real_epi_prob_and_trigger.csv | tr -s ',' '\n' | nl -nln |  grep "epi_prob" | cut -f1)

# Print the first n rows of that column
cut -d"," -f$col_num2 county_real_epi_prob_and_trigger.csv | \
sort -k1 -nr | awk '!/NA/' | head -20

# Gives total rows in column that are not NA
cut -d"," -f$col_num county_real_epi_prob_and_trigger.csv | \
sort -k1 -nr | awk '!/NA/' | wc -l
