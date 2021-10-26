# Snakemake workflow: CLIP

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.2.1-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/CLIP.svg?branch=master)](https://travis-ci.org/snakemake-workflows/CLIP)

This is the template for a new Snakemake workflow. Replace this text with a comprehensive description covering the purpose and domain.
Insert your code into the respective folders, i.e. `scripts`, `rules`, and `envs`. Define the entry point of the workflow in the `Snakefile` and the main configuration in the `config.yaml` file.

## Authors

* Ohler (@lebedeva)

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.tsv` to specify your sample setup.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake mamba
    conda activate snakemake
    mamba install snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores. 

To submit a job that runs snakemake, you can use `run.sh`, which contains some sensible default parameters for the MDC max-cluster.

    qsub run.sh 

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.


## Testing

Test cases are in the subfolder `.test`. They are automatically executed via continuous integration with [Github Actions](https://github.com/features/actions).


# Pipeline steps

### Step 1. Adapter trimming

### Step 2. Collapsing identical reads

Modifiable depending on whether your data has UMIs or not (older datasets will often not have this).  

#### Option a. You have UMIs: use rule umi-tools.

#### Option b. You do not have UMIs: use rule collapse-reads.

This will collapse all identical reads into one.

### Step 3. Align reads.

#### Option a. We use Segemehl aligner because we found it to be suitable for short CLIP reads. Find details of aligner comparison for simulatted data here: PMID: 26776207.

#### Option b. Use a combination of STAR and segemehl.

### Step 4. Call peaks.

Ohler lab has 2 peak callers: 

#### PARALYZER:

##### Input: a bam file of mapped CLIP reads.

##### Output: a bed file of called peaks.

#### omniCLIP:

##### Input:

- bam files of at least 2 replicates of mapped CLIP reads
- at least one bam file of the background sample (it could be an input sample matched to your CLIP, or RNA-seq for the respective cell line/organism)

##### Output: a bed file of called peaks.

The use of peak caller depends on the type of CLIP data and whether you have replicates.

#### Option a. You have PAR-CLIP data and no replicates. Use PARALYZER.

#### Option b. You have PAR-CLIP data and replicates. Use either PARALYZER or omniCLIP.

#### Option c. You have other CLIP (eCLIP, iCLIP, ...) that does not rely on 4-thiouridine incorporation and no replicates. Ohler lab cannot help you.

#### Option d. You have non-PAR CLIP data and replicates. Use omniCLIP.


