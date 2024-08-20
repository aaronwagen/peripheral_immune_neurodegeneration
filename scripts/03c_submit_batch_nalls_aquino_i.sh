#!/bin/bash

# Define the number of eQTL datasets
NUM_EQTLS=1500  # Update this with the actual number of eQTL datasets

# Loop over the eQTL datasets and submit each as a separate job
for ((i=1; i<=NUM_EQTLS; i++)); do
    sbatch --job-name="eqtl_$i" \
           --output="aquino_nalls_logs/eqtl_${i}.out" \
           --error="aquino_nalls_logs/eqtl_${i}.err" \
           --time=1:30:00 \
           --mem=16G \
           --partition=ncpu \
          02f_submit_aquino_nalls_coloc_job.sh $i
done