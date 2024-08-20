#!/bin/bash

# Load R module, if necessary
ml  R/4.2.0-foss-2021b

cd /nemo/lab/gandhis/home/users/wagena/SNCA/

# Run the R script with the eQTL number as an argument
Rscript r_scripts/02i_run_coloc_gwas_immune_aquino_rizig.R $1