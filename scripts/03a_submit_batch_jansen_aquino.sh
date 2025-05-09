#!/bin/bash

# Define the number of eQTL datasets
NUM_EQTLS=2904  # Update this with the actual number of eQTL datasets

# Loop over the eQTL datasets and submit each as a separate job
for ((i=1; i<=NUM_EQTLS; i++)); do
    sbatch --job-name="qtl_$i" \
           --output="aquino_jansen_logs/eqtl_${i}.out" \
           --error="aquino_jansen_logs/eqtl_${i}.err" \
           --time=3:00:00 \
           --mem=20G \
           --partition=ncpu \
          02k_submit_aquino_jansen_coloc_job.sh $i
done