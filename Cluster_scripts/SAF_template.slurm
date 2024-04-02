#!/bin/bash
#SBATCH -t 24:00:00                 # walltime = 24 hrs
#SBATCH --mincpus=10                # Request 10 CPUs
#SBATCH --job-name=alevinfry        # Job name
#SBATCH --output=alevin.%j.out      # Output file
#SBATCH --error=alevin.%j.err       # Error file
# Go to working directory
cd /om/user/azding
# Activate conda environment
source /home/azding/anaconda3/bin/activate /home/azding/anaconda3/envs/af
# Set directory paths
export AF_SAMPLE_DIR=$PWD/af_test_workdir
cd $AF_SAMPLE_DIR
FASTQ_DIR="$AF_SAMPLE_DIR/data/eightsras"
export ALEVIN_FRY_HOME="$AF_SAMPLE_DIR/af_home"
simpleaf set-paths
ulimit -n 2048
cd $AF_SAMPLE_DIR
REF_DIR="$AF_SAMPLE_DIR/data/refdata-gex-GRCh38-2020-A"
IDX_DIR="$AF_SAMPLE_DIR/%%INDEX%%/human-2020-A_splici3"
# Identify reads
reads1_pat="%%SRAID%%_2.fastq"
reads2_pat="%%SRAID%%_3.fastq"
reads1="$(find -L ${FASTQ_DIR} -name "*$reads1_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"
reads2="$(find -L ${FASTQ_DIR} -name "*$reads2_pat*" -type f | sort | awk -v OFS=, '{$1=$1;print}' | paste -sd, -)"
# Run simpleaf quant
simpleaf quant \
--reads1 $reads1 \
--reads2 $reads2 \
--threads 16 \
--index $IDX_DIR/index \
--chemistry 10xv3 --resolution cr-like \
--expected-ori fw --unfiltered-pl \
--t2g-map $IDX_DIR/index/t2g_3col.tsv \
--output $AF_SAMPLE_DIR/%%SRAID%%_quant_v3

