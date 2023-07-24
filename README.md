# Pipeline overview

This pipeline is to map reads and call peaks for CLIP data. Our lab developed two peak callers that can be used in this workflow: [PARALYZER](https://ohlerlab.mdc-berlin.de/software/PARalyzer_85/) and [omniCLIP](https://github.com/ohlerlab/omniCLIP).

We are grateful to the people who created containers for those tools (Biocontainers for everything and Zavolan lab for omniCLIP).

## Authors

* Svetlana Lebedeva (@lebedeva), Senior postdoc, Ohler group
* René Schmiedler (@reschmi), intern, Ohler group

# Prerequisites

This pipeline is configured to run on a cluster with Sun Grid Engine queuing system. It requires snakemake and singularity to run.

# Description of the input

1. CLIP raw reads: a single end fastq.gz from Illumina (file for testing provided in `test_data/CLIP.fastq.gz`).
2. Reference genome: 
    1. A single .fa file (specify path in `config.yaml` REFERENCE_GENOME; use [GRCh37.p13.genome.fa](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz) for the test data)
    2. [omniCLIP only] A directory with a fasta.gz file for each chromosome. (specify path in `config.yaml` GENOME_DIR). You can use [UCSC faSplit](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSplit) to split your reference fasta by chromosome  like this `faSplit byname GRCm38.p6.genome.fa [GENOME_DIR]`.
    3. An annotation: gtf and [omniCLIP only] gff (specify path in `config.yaml` GTF and GFF; use [gencode.v19.annotation.gtf](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz) and [gencode.v19.annotation.gff3](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gff3.gz) for the test data)
3. [omniCLIP only] A background bam file (provided in `test_data/test_backgr.bam`). We use total RNA-seq; make sure that it is indexed. Also make sure that it is the same strandedness as the CLIP bam (standard RNA-seq is sequencing the reverse strand and small RNA-seq the forward strand, in this case you need to invert the strand of the background bam file for peak calling to work correctly).

## Usage

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and, if available, its DOI (see above).

### Step 1: Obtain a copy of this workflow

1. Create a new github repository using this workflow [as a template](https://help.github.com/en/articles/creating-a-repository-from-a-template).
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the newly created repository to your local system, into the place where you want to perform the data analysis.

```
git clone https://github.com/ohlerlab/clip_pipeline
cd clip_pipeline/
```

### Step 2: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake mamba
    conda activate snakemake
    mamba install snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 3: Install additional packages

Install singularity, if it is not available in your system (make sure that it is version 3.x):

    mamba install -c conda-forge singularity
 
### Step 4: Execute workflow on test data


Test your configuration by performing a dry-run via

   bash test_run.sh --dry-run 

Execute the test workflow locally via

   bash test_run.sh 

It will download sequence for human chr1 and hg19 gtf/gff3 annotation. 

### Step 5: Execute workflow on your data

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.tsv` to specify your sample setup.

Change the paths in `config/config.yaml` to the actual path for your genome and annotation:

```
## Reference genome and annotation:
REFERENCE_GENOME:
  "GRCh37.p13.genome.fa"
GFF:
  "gencode.v19.annotation.gff3"
GTF:
  "gencode.v19.annotation.gtf"
```

Specify your input files in "config/samples.tsv"

Provide a background RNA-seq bam file if you would like to run omniCLIP. Make sure it has the same strandedness as the CLIP input file:

```
## RNA-seq background for omniCLIP:
BACKGROUND:
  "test_data/test_backgr.bam"
```

Provide directories to mount to containers in 'mounts.txt', if necessary.

You can run the pipeline locally:

    snakemake --use-singularity

Note that you will need at least 30G RAM ([recommended 64G](https://www.bioinf.uni-leipzig.de/Software/segemehl/)) for generating the segemehl index of human genome. 

To submit a job that runs snakemake to SGE cluster, you can use `run.sh`, which contains some sensible default parameters for an SGE queueing system.

    qsub run.sh


See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

---

# Expected output

The output is in the folder `results`. The `[sample]` is the sample name you gave in `samples.tsv`.

## Aligned reads

`results/prepare_aligned/[sample]_sorted_deduplicated.bam` contains the aligned reads from the CLIP sample.

## Called peaks

### PARALYZER

`results/call_peaks` contains output from PARALYZER. 

You are mostly interested in `[sample].clusters` file which has the following structure:

```
Chromosome,Strand,ClusterStart,ClusterEnd,ClusterID,ClusterSequence,ReadCount,ModeLocation,ModeScore,ConversionLocationCount,ConversionEventCount,NonConversionEventCount
chr1,+,86048856,86048921,G7.1,CTTGAGGAGCATTAAGGTATTTCGAAACTGCCAAGGGTGCTGGTGCGGATGGACACTAATGCAGCC,63,86048905,0.9407689813401,12,16,305
chr1,+,86048970,86048992,G7.2,TTCATTTTGGAGCTTGTGGAGTT,21,86048982,0.8242463078856529,6,9,150
chr1,+,86048998,86049001,G7.3,CTTT,22,86049000,0.5344784636094164,2,2,60
```

See [readme](https://ohlerlab.mdc-berlin.de/files/duke/PARalyzer/README.txt) for the description of the complete output.

Note that `config/config.yaml` contains settings for PARALYZER peak caller parameters which you can modify to achieve desired peak calling results.


### omniCLIP

`results/omniclip/output_[sample]` contains output of omniCLIP:

1. `pred.txt` - this is the raw peaks output

- Gene : gene
- ChrName   : chr
- Start : peak start
- Stop : peak end
- Strand : peak strand
- SiteScore : the log-likelihood ratio of the peak state versus the other states
- Coverage : peak coverage
- A.A A.C A.G A.T A.D C.A C.C C.G C.T C.D G.A G.C G.G G.T G.D T.A T.C T.G T.T T.D : counts of mutations and deletions 
- mean_mat_fg, var_mat_fg, mean_mat_bg, var_mat_bg : the means and variances of the Negative binomial models for modeling the read coverage 
- counts_fg, counts_bg : read coverage in foreground (CLIP) and background 
- pv : the log10 of the p-value
- max_pos : position of maximum coverage
- dir_score: the likelihood of the Dirichlet model for modeling the diagnostic events


2. `pred.bed` - these are the final selected peaks

Score : SiteScore normalized to be maximum 1000:
Score = Score/ ( max(SiteScore) * 0.001 )

ThickStart and ThickEnd: position of the maximum coverage

### Visualize the results

An igv screenshot after successful run on the test data:
<img src="https://github.com/ohlerlab/clip_pipeline/blob/main/test_data/expected_output.png">
[https://github.com/ohlerlab/clip_pipeline/blob/omniclip/test_data/expected_output.png](https://github.com/ohlerlab/clip_pipeline/blob/main/test_data/expected_output.png)

Navigate to chr1:86,047,775-86,050,811 (or any of the regions in `test_data/regions.bed` if you used whole genome reference) to check the results. 

Load the following files into igv:

- `results/omniclip/pred.bed` is the peaks called by omniCLIP
- `results/prepare_aligned/test_sorted_deduplicated.bam` is the CLIP bam file used for peak calling
- you can create a bed file from PARALYZER output as follows: `awk '{FS=","}{OFS="\t"}{print $1,$3,$4,".",".",$2}' results/call_peaks/test.clusters | tail -n+2 > results/call_peaks/test_clusters.bed`


# Detailed description of steps

Part 1. Map reads (implemented)

- Input: fastq
- Output: bam
- Parameters:
    - Adapter sequences
    - config: has_UMI (T/F)
    - aligner: segemehl (default) (STAR, bowtie - not implemented)
    
Part 2. Call peaks (PARalyzer and omniCLIP)

- Input: bam, omniCLIP requires background = total RNA sequencing bam; omniCLIP requires genome direstory with each chromosome as a separate fasta.gz.
- Output: bed file with peak coordinates and scores
- Parameters: peak caller (PARALYZER, omniCLIP)



# CLIP description

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

Crosslinking of the RNA-binding protein to the RNA affects reverse transcription step during library preparation, because the crosslinked base will induce base miscorporation or termination of reverse transcriptase (RT) with increased frequency compared to non-crosslinked bases. There are two main variants of CLIP library preparation:

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

#### Option a. You have PAR-CLIP data but no backgound RNA-seq. Use PARALYZER.

#### Option b. You have PAR-CLIP data and backgound RNA-seq. Use either PARALYZER or omniCLIP.

#### Option c. You have non-PAR CLIP data and background RNA-seq. Use omniCLIP.

#### Option d. You have non PAR-CLIP data and no background RNA-seq. Ohler lab cannot help you.





# Known bugs

## omniCLIP

- Caution: CLIP bam file HAS TO have reads mapping to EVERY chromosome of the reference genome. If the data is very sparse and there is a chromosome with no reads mapping to it, the pipeline will fail. This is one of the reasons we do not currently use scaffolds but only main chromosomes.

