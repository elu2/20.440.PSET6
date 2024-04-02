#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --mincpus=8
#SBATCH --job-name=SCVelo   # Job name
#SBATCH --output=./slurmOut/SCV.%j.out   # Output file
#SBATCH --error=./slurmOut/SCV.%j.err    # Error file
#SBATCH --partition=sched_any

cd /home/ericjlu/20.440Project/

python3 RNAVelo.py %%SRAID%%