rule cutadapt:
    input:
      fastq=lambda wc : samples.loc[wc.sample, "path"],
    output:
      fastq="results/process_reads/{sample}_trim.fastq.gz",
    conda:
      "../envs/cutadapt.yaml"
    singularity:
      "docker://quay.io/biocontainers/cutadapt:3.5--py39h38f01e4_0"
    log: 
      "results/logs/cutadapt_{sample}.log" 
    threads: 2
    params:
      repeat_n_times=2,
      match_read_wildcards=18
    shell:
      """
     cutadapt -a {config[THREE_PRIME_ADAPTER_SEQUENCE]} -g {config[FIVE_PRIME_ADAPTER_SEQUENCE]}  \
                -b {config[NINETEEN_NUCLEOTIDE_MARKER_SEQUENCE]} -b {config[TWENTY_FOUR_NUCLEOTIDE_MARKER_SEQUENCE]}  \
                -b {config[KIT_5ADAPTER]} -b {config[ICLIP_ADAPTER]} -b {config[HITSCLIP_3ADAPTER]} -b {config[HITSCLIP_5ADAPTER]}  \
                -b {config[SPECIAL_AD_1]} -b {config[SPECIAL_AD_2]} -b {config[SPECIAL_AD_3]}  \
                -j {threads} -n {params.repeat_n_times} -m {params.match_read_wildcards} -o {output.fastq} {input.fastq} 2> {log}
      """

rule extract_UMIs:
    input:
      fastq="results/process_reads/{sample}_trim.fastq.gz",
    output:
      processed="results/process_reads/{sample}_trim_umi-extr.fastq.gz"
    singularity:
      "docker://quay.io/biocontainers/umi_tools:1.1.2--py310h1425a21_1"
    log:
      "results/logs/extract_UMIs_{sample}.log"
    resources:
      mem_mb_per_cpu='8G'
    params:
      pattern=lambda wc: config['UMI-BARCODE']
    shell:
      """
      umi_tools extract --stdin={input.fastq} --extract-method=regex --bc-pattern={params.pattern:q} --log={log} --stdout {output.processed}
      """


rule collapse_all_reads:
    input:
      fastq="results/process_reads/{sample}_trim.fastq.gz",
    output:
      collapsed="results/process_reads/{sample}_trim_collapsed.fastq.gz",
    conda:
      "../envs/fastx_toolkit.yaml"
    singularity:
      "docker://biocontainers/fastxtools:v0.0.14_cv2"
    log:
      "results/logs/collapse_reads_{sample}.log"
    shell:
      """
      zcat {input.fastq} | fastx_collapser | gzip > {output.collapsed} 2> {log}
      """