def get_mem_mb(wildcards, attempt):
  return str(attempt**2 * 16) + 'G'

"""Pre-align reads to ncRNA using bowtie"""

from pathlib import Path

# Bowtie indexes
IDX_NUMS = [1, 2, 3, 4, "rev.1", "rev.2"]
BT_IDX_PATH = "results/align_reads/bowtie_idx_ncRNA"
BT_IDXS = [f"{BT_IDX_PATH}.{x}.ebwt" for x in IDX_NUMS]

rule bowtie_idx:
    """Build bowtie index from ncRNA sequences.
    """
    input:
      fasta = config["NCRNA"],
    output:
      BT_IDXS,
    conda:
      "../envs/bowtie.yaml"
    container:
      "docker://quay.io/biocontainers/bowtie:1.3.1--py312hf8dbd9f_10"
    params:
        outdir = "results/align_reads/bowtie_idx_ncRNA",
    log:
        "results/logs/bowtie_idx_ncRNA.log",
    message:
        "Building bowtie index for ncRNA"
    shell:
        """
        bowtie-build {input.fasta} {params.outdir} 2> {log}
        """

rule bowtie_align_ncRNA:
    """ Align fastqs to ncRNA fasta.
        """ 
    input:
      fastq=lambda wildcards: "results/process_reads/{sample}_trim_umi-extr.fastq.gz".format(sample=wildcards.sample) if config['UMI-BARCODE'] is not False else "results/process_reads/{sample}_trim_collapsed.fastq.gz".format(sample=wildcards.sample),
      idx = BT_IDXS
    output:
      bam="results/align_reads/{sample}_ncRNA.bam",
      unmapped="results/align_reads/{sample}_umi_trim_noncRNA.fastq.gz"
    params:
      idx_path = BT_IDX_PATH,
    conda:
      "../envs/bowtie.yaml"
    container:
      # Contains bowtie 1 and samtools 1
      "docker://quay.io/biocontainers/mulled-v2-3ffca7af725971ca3232b6e98526cfb5ed81828c:76f24b7c13470ff9abbe70f3dcfabf8961dd2a60-0"
    threads: config["THR"]
    log:
        "logs/bowtie_ercc_rRNA/{sample}.log",
    message: "Aligning {wildcards.sample} to ncRNA"
    shell:
      """
        bowtie -p {threads} -q \
        --best --sam \
        --un >(gzip -c > {output.unmapped}) \
        {params.idx_path} \
        <(gunzip -c {input.fastq}) \
        | samtools view -bhSF4 - \
        | samtools sort - -o {output.bam} 
        
        samtools index {output.bam} 
      """

# TODO call peaks on those too

rule segemehl_idx:
    input:
      fasta=config["REFERENCE_GENOME"],
    output:
      sege_idx=config["REFERENCE_GENOME"] + ".idx",
    conda:
      "../envs/segemehl.yaml"
    container:
      "docker://quay.io/biocontainers/segemehl:0.3.1--h39379e4_3"
    log:
      "results/logs/segemehl_idx.log"
    resources:
      mem_mb_per_cpu=get_mem_mb
    threads: config["THR"]
    retries: 2
    shell:
      """
      segemehl.x -t {threads} -x {output.sege_idx} -d {input.fasta} 2> {log}
      """


rule segemehl:
    input:
      fastq=rules.bowtie_align_ncRNA.output.unmapped,
      sege_idx=rules.segemehl_idx.output.sege_idx,
      ref=config["REFERENCE_GENOME"]
    output:
      sam="results/align_reads/{sample}_aligned.sam",
      unmapped="results/align_reads/{sample}_segemehl_unmapped.fastq"
    conda:
      "../envs/segemehl.yaml"     
    container:
      "docker://quay.io/biocontainers/segemehl:0.3.1--h39379e4_3"
    log:
      "results/logs/segemehl_{sample}.log"
    threads: config["THR"]
    resources:
      mem_mb_per_cpu='8G'
    params:
      max_occurences=1,
      differences=2
    shell:
      """
      segemehl.x -S -D {params.differences} -M {params.max_occurences} --briefcigar -t {threads} -i {input.sege_idx} -d {input.ref} -q {input.fastq} -u {output.unmapped} > {output.sam} 2> {log}
      """


