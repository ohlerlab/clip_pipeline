#!/bin/bash

# script that will run snakemake on the max cluster

#$ -V
#$ -cwd
#$ -j yes
#$ -l longrun
#$ -l m_mem_free=4G
#$ -l h_rt=168:0:0

test -d logs/cluster || { >&2 echo "logs/cluster does not exist"; exit 1; }

eval "$(conda shell.bash hook)"

conda activate snakemake

if [ -f mounts.txt ]; then
    >&2 echo "mounts.txt file is present..."
    while read s t; do
        MOUNT="${MOUNT}-B ${s}:${t} "
    done < mounts.txt
    >&2 echo "Adding the following mounts to singularity : \"${MOUNT}\""
fi

export SGE_ROOT="/opt/uge"

# Start snakemake
snakemake --snakefile workflow/Snakefile \
          --use-singularity \
          --singularity-args "--nv ${MOUNT}" \
          --cluster "qsub -V -cwd -pe smp {threads} -l m_mem_free={resources.mem} -l h_rt {resources.runtime} {resources.misc} -j yes " \
          --default-resources mem="4G" runtime="2:0:0" \
          --directory "${PWD}" \
          --jobs 100 \
          --latency-wait 30 \
          --keep-going \
          "$@"
          
