
# if the file creation date changes during the job, make the names consistent

for file in *31.csv; do mv "$file" "${file/31.csv/30.csv}"; done
for file in *31.rda; do mv "$file" "${file/31.rda/30.rda}"; done

for file in *08-01.csv; do mv "$file" "${file/08-01.csv/07-30.csv}"; done
for file in *08-01.rda; do mv "$file" "${file/08-01.rda/07-30.rda}"; done

ls *30.csv | wc -l  # 15721 # on 08-01-2023
ls *30.rda | wc -l  # 15721 

# new expected is 16315 (og imports) + 6275 (new imports) = 22590
