#! /bin/bash 

# To run,
# sbatch --array 1-20 cluster_scripts/stability_kmeans.sh
source activate goldwrap

SEED=${SLURM_ARRAY_TASK_ID}
for N_CLUSTERS in {5..20};
do
    Rscript run_clustering.R ${N_CLUSTERS} ${SEED}
done;
