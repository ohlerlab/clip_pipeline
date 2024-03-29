def get_mem_mb(wildcards, attempt):
  return str(attempt**2 * 16) + 'G'

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
      fastq=lambda wildcards: "results/process_reads/{sample}_trim_umi-extr.fastq.gz".format(sample=wildcards.sample) if config['UMI-BARCODE'] is not False else "results/process_reads/{sample}_trim_collapsed.fastq.gz".format(sample=wildcards.sample),
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
