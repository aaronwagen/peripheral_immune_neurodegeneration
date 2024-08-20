#!/bin/bash

# Define the specific eQTL dataset indices you want to process
specific_indices=(474 554 1465 2166 2195)


# Loop over the eQTL datasets and submit each as a separate job
for i in "${specific_indices[@]}"; do
    sbatch --job-name="eqtl_$i" \
           --output="aquino_rizig_logs/qtl_${i}.out" \
           --error="aquino_rizig_logs/qtl_${i}.err" \
           --time=4:00:00 \
           --mem=16G \
           --partition=ncpu \
          02i_submit_aquino_rizig_coloc_job.sh $i
done
