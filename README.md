# 20.440.PSET6

## Overview
This repo contains code to produce a UMAP plot with velocity trajectories. Much of the code leans on scVelo, a popular python package used to analyze AnnData-type objects for single-cell RNA (scRNA) velocity. 

Then, there is a very hand-wavy, tentative method for assigning cell types. Cells are first clustered with Louvain clustering and t-tests are performed between each cluster and the rest of the cells (standard). $\log_2$ fold-change of genes at the cluster-level are also reported to give a quick measure for differences in expression of cells between clusters. 

My tentative method for assigning cell-types relies on a list of marker genes in each cell-type defined by [1]. The cell type maximizing the sum of $\log_2$ fold-changes in the cell-type's associated list of marker genes for a given cluster is assigned to the cluster. If the sum of $\log_2$ fold-changes is less than 0, then the cluster is assigned as "Misc." This method is intended to be a sufficient temporary solution and suffers from a lack of statistical rigor for the sake of expediency.


## Data
The data being used are single-cell RNA sequencing samples from [1], studying the single-cell transcriptome in endometriosis. We are currently focusing on 8 samples split between 2 individuals -- each individual with 2 endometrial samples and 2 control samples.

FASTQ files were retrieved using SRAToolkit and converted into AnnData-type objects with Alevin-fry [2]. Both required the use of cluster computing and FASTQ file sizes were ~50GB and AnnData objects were ~5-15GB. As such, we will not be uploading our data to GitHub. Our initial 8 samples were retrieved with the following SRA accession IDs: SRR2638397[2-5] and SRR2638398[2-5]. After the data retrieval process, a large directory (`af_quant`) suitable for reading by scVelo is returned. This repo assumes that this hypothetical `af_quant` exists in the home directory.


## Folder structure
The `Cluster_scripts` directory contains slurm scripts (`*.slurm`) and shell scripts (`*.sh`) for retrieving data and running analyses. RNAVelo.py is used for downstream data analysis and generates figures (including the one for this assignment).

`Figures` contains (one of) the figures that RNAVelo generates.

The hypothetical `af_quant` has a complicated structure that is not relevant to know during analysis -- scVelo asks only for the directory's path as an argument to read in. This directory would be stored in the home directory.


## Installation
First, SRAToolkit and Alevin-fry should be installed. One slurm script uses SRAToolkit to prefetch, then dump FASTQ files. Another slurm script uses alevin-fry to preprocess the FASTQ into analyzable AnnData structure. Each has an accompanying shell script to automatically generate slurm scripts with different SRA IDs. Paths in each file will not make sense or be consistent as they are unchanged from how they were run on a cluster.

All requisite python packages for downstream analyses are enumerated in the `requirements.txt` file under Python version 3.11.3. Generating the figure only requires running the following command `python3 RNAVelo.py [SRA########]`. The figure in `Figures` uses SRR26383984.


## References
1. Fonseca, M.A.S., Haro, M., Wright, K.N. et al. Single-cell transcriptomic analysis of endometriosis. Nat Genet 55, 255–267 (2023). https://doi.org/10.1038/s41588-022-01254-1
2. He, D., Zakeri, M., Sarkar, H., Soneson, C., Srivastava, A., and Patro, R. Alevin-fry unlocks rapid, accurate and memory-frugal quantification of single-cell RNA-seq data. Nat Methods 19, 316–322 (2022). https://doi.org/10.1038/s41592-022-01408-3
