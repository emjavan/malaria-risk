#!/bin/bash

#SBATCH -J new_r0_import                 # Job name
#SBATCH -o new_r0_import.%j.o            # Name of stdout output file (%j expands to jobId)
#SBATCH -e new_r0_import.%j.e            # Name of stderr output file (%j expands to jobId)
#SBATCH -p normal                        # Queue name, small is for <=2 nodes, normal 3+
#SBATCH -N 48                  	         # Total number of nodes requested
#SBATCH -n 2600                          # Total number of tasks to run 56 cores/node (28 per socket)
#SBATCH -t 48:00:00            	         # Run time (hh:mm:ss)
#SBATCH -A A-ib1                         # Allocation name
#SBATCH --mail-user=emjavan@utexas.edu   # Email for notifications
#SBATCH --mail-type=all                  # Type of notifications, begin, end, fail, all

# Load newest R module on Frontera
module load Rstats/4.0.3

# Load launcher
module load launcher

# Configure launcher
EXECUTABLE=$TACC_LAUNCHER_DIR/init_launcher
PRUN=$TACC_LAUNCHER_DIR/paramrun
CONTROL_FILE=commands_malaria.txt
export LAUNCHER_JOB_FILE=commands_malaria.txt
export LAUNCHER_WORKDIR=`pwd`
export LAUNCHER_SCHED=interleaved

# Start launcher
$PRUN $EXECUTABLE $CONTROL_FILE
