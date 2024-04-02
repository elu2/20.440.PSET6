#!/bin/bash
#SBATCH -t 18:00:00         
#SBATCH --mincpus=4
#SBATCH --job-name=SRA   # Job name
#SBATCH --output=SRA.%j.out   # Output file
#SBATCH --error=SRA.%j.err    # Error file
# Go to working directory
cd /om/user/azding
# Run SRA Toolkit to extract fastq files
sratoolkit.3.1.0-ubuntu64/bin/prefetch -X 9999999999999 %%SRAID%% 
sratoolkit.3.1.0-ubuntu64/bin/fasterq-dump %%SRAID%% --include-technical -S
