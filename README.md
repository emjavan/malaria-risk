# Estimate the risk/probability of an undetected Malaria epidemic in each county of continental US
The stochastic SEIR model begins each simulation with one imported infectious individual at day 1 of their infectious period. The simulation terminates when there are no more infectious individuals or the total autochthonous infections reaches 1,000. The model assumes an unlimited pool of susceptible individuals, often the case early in an emerging epidemic, thus will never self-limit. Model parameters and their source are provided below.

Table 1. 
| Parameter   | Assumption  | Source        |
| :---        |    ----:    |          ---: |
| Generation time     | Title       | Here's this   |
| Exposure period    | Text        | And more      |
| Infectious period      | Title       | Here's this   |
| Min day infectious   | Text        | And more      |
| Min day     | Title       | Here's this   |
| Paragraph   | Text        | And more      |
| Header      | Title       | Here's this   |
| Paragraph   | Text        | And more      |
| Header      | Title       | Here's this   |
| Paragraph   | Text        | And more      |
| Header      | Title       | Here's this   |
| Paragraph   | Text        | And more      |





# Ingest R0 and probability of next importation from input_data, then write commands file
`clean_r0_import_write_commands.R`   

# Run simulations on TACC Frontera

