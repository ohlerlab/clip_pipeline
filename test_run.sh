#!/bin/bash

# script that will run snakemake on the max cluster

#$ -V
#$ -m ea
#$ -cwd
#$ -j yes
#$ -l m_mem_free=4G
#$ -l h_rt=4:0:0
#$ -o logs/

#test -d logs/cluster || { >&2 echo "logs/cluster does not exist"; exit 1; }
mkdir -p logs/

eval "$(/${HOME}/miniconda3/bin/conda shell.bash hook)"
conda activate snakemake

if [ -f mounts.txt ]; then
    >&2 echo "mounts.txt file is present..."
    while read s t; do
        MOUNT="${MOUNT}-B ${s}:${t} "
    done < mounts.txt
    >&2 echo "Adding the following mounts to singularity : \"${MOUNT}\""
fi

#export SGE_ROOT="/opt/uge"

## get test genome and annotation

if [ ! -f chr1.fa ]; then
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz &&
    gunzip chr1.fa.gz 
fi

if [ ! -f gencode.v19.annotation.gtf ]; then
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz &&
    gunzip gencode.v19.annotation.gtf.gz
fi
if [ ! -f gencode.v19.annotation.gff3 ]; then
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz &&
    gunzip gencode.v19.annotation.gff3.gz 
fi


# Start snakemake locally
snakemake -rp --snakefile workflow/Snakefile \
          --use-singularity \
	      --singularity-args "--nv ${MOUNT}" \
          --directory "${PWD}" \
          --jobs 100 \
          --rerun-incomplete \
          --latency-wait 30 \
          --keep-going \
          --show-failed-logs \
          "$@"
          
