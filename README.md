# Snakemake workflow: CLIP

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.2.1-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/CLIP.svg?branch=master)](https://travis-ci.org/snakemake-workflows/CLIP)


## Authors

* Ohler (@lebedeva)

# CLIP Pipeline description

CLIP is a series of methods of determining RNA-binding protein (RBP) binding sites on the transcriptome. It uses UV crosslinking to covalently attach RBPs to RNA, which can be enhanced using nucleotide analogs such is 4-thiouridine (called PAR-CLIP). The RNA-protein complexes are digested with RNase to create a short footprint, enirched by immunoprecipitation of the RBP and sequenced. This pipeline provides tools for mapping of Illumina-sequenced CLIP libraries and calling the RBP binding sites.


### Step 1. Adapter trimming

rule: cutadapt

Different CLIP protocols use custom sequencing adapters, as well as standard Illumina small RNA sequencing adapters. We tried to include most common adapter sequences into the pipeline, so this step does not have to be modified if multiple different CLIP experiments need to be processed together. It is always advisable to check the output of the adapter trimming. If most reads did not contain one of the adapters, the adapter sequence is possibly missing from the pipeline and needs to be added in the config file.

### Step 2. Collapsing identical reads

rule: collapse_reads

Identical reads which are an artifact of PCR duplication need to be removed before further processing. Modern CLIP libraries usually have unique molecular identifiers (4-7 N bases at the end of the adapters) which help to identify if reads with identical sequence were independent molecules before the PCR. UMI-tools can then be used to deduplicate the reads. Older CLIP libraries do not contain UMIs, and the most conservative way in this case is to collapse all idenitcal reads. The latter step is implemented in this pipeline to make it compatbile with the older libraries. 

To Do: add umi-tools rules  

#### Option a. You have UMIs: use rule umi-tools.

#### Option b. You do not have UMIs or you do not know: use rule collapse_reads.


### Step 3. Align reads.

Crosslinking of the RNA-binding protein to the RNA affects reverese transcription step during library preparation, because the crosslinked base will induce base miscorporation or termination of reverse transcriptase (RT) with increased frequency compared to non-crosslinked bases. There are two main variants of CLIP library preparation:

1. Two adapters are ligated to the RNA on both ends. Reverse transcriptase reads through and only full length cDNA is sequenced. RBP binding site is in the middle of the read and is determined by base mutations and deletions. (PAR-CLIP, HITS-CLIP)

2. One adapter is ligated at the 3'end of the RNA and reverse transcription is performed. The cDNA is cirularized and re-cut or the 5'adapter is ligated to the cDNA. RBP binding site is at the end of the read (RT stop). (iCLIP, iCLIP2, eCLIP)

Our pipeline is focused on the first, readthrough type of CLIP libraries. Because the reads are short and often have mutations, the optimal alignment is not trivial. Several aligners can be used:

#### Option a. Segemehl.

rule: segemehl

Segemehl aligner is designed to align short reads with deletions and mismatches so we use it by default for the CLIP reads. Find details of aligner comparison for simulatted data here: PMID: 26776207.

#### Option b. STAR.

STAR is designed for longer, especially paired reads. This would be a more appropriate aligner for the RT stop type of CLIP libraries. 

#### Option c. bowtie.

Bowtie aligner was used previously for the readthrough kinds of CLIP libraries. It will output less alignments than segemehl.

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
    snakemake --use-singularity -n

Execute the workflow locally via

    snakemake --use-conda --cores $N
    snakemake --use-singularity --cores $N

using `$N` cores. 

To submit a job that runs snakemake, you can use `run.sh`, which contains some sensible default parameters for an SGE queueing system.

    qsub run.sh 

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.


## Testing

Test cases are in the subfolder `.test`. They are automatically executed via continuous integration with [Github Actions](https://github.com/features/actions).



