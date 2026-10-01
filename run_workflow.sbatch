#!/usr/bin/env bash

#SBATCH --job-name=CutandTag_ReplicatePeak
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2


# Modules needed by Snakemake
module load python/3.10
module load pandas/2.2.3
module load numpy/1.22.3
module load matplotlib/3.7.1

# Submit workflow jobs using config/cluster_config.yml.
snakemake \
    -j 50 \
    --use-envmodules \
    --latency-wait 300 \
    --cluster-config config/cluster_config.yml \
    --cluster "sbatch \
        -A {cluster.account} \
        -p {cluster.partition} \
        --cpus-per-task {threads} \
        -t {cluster.time} \
        --mem {cluster.mem} \
        --job-name {cluster.name} \
        --output {cluster.output} \
        --error {cluster.error}" \
    "$@"
