#!/bin/bash

# Define the number of eQTL datasets
NUM_EQTLS=2258  # Update this with the actual number of eQTL datasets

# Loop over the eQTL datasets and submit each as a separate job
for ((i=2257; i<=NUM_EQTLS; i++)); do
    sbatch --job-name="qtl_$i" \
           --output="aquino_shigemizu_logs/eqtl_${i}.out" \
           --error="aquino_shigemizu_logs/eqtl_${i}.err" \
           --time=4:00:00 \
           --mem=16G \
           --partition=ncpu \
          02l_submit_aquino_shigemizu_coloc_job.sh $i
done