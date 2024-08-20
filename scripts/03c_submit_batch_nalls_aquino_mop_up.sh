#!/bin/bash

# Define the specific eQTL dataset indices you want to process
specific_indices=(201 2278 2290 2310 399 448 482)


# Loop over the eQTL datasets and submit each as a separate job
for i in "${specific_indices[@]}"; do
    sbatch --job-name="eqtl_$i" \
           --output="aquino_nalls_logs/eqtl_${i}.out" \
           --error="aquino_nalls_logs/eqtl_${i}.err" \
           --time=4:00:00 \
           --mem=16G \
           --partition=ncpu \
          02f_submit_aquino_nalls_coloc_job.sh $i
done
