rule fix_header:
    input:
      fastq=lambda wc : samples.loc[wc.sample, "path"],
    output:
      fastq="results/process_reads/{sample}_fix.fastq.gz",
    shell:
      """
      zcat {input.fastq} | sed s'/ /:/g' | gzip > {output.fastq}
      """

rule cutadapt:
    input:
      fastq="results/process_reads/{sample}_fix.fastq.gz",
    output:
      fastq="results/process_reads/{sample}_trim.fastq.gz",
    conda:
      "../envs/cutadapt.yaml"
    container:
      "docker://quay.io/biocontainers/cutadapt:3.5--py39h38f01e4_0"
    log: 
      "results/logs/cutadapt_{sample}.log" 
    threads: 2
    params:
      three_prime_adapter_sequence=config["THREE_PRIME_ADAPTER_SEQUENCE"],
      five_prime_adapter_sequence=config["FIVE_PRIME_ADAPTER_SEQUENCE"],
      nineteen_nucleotide_marker_sequence=config["NINETEEN_NUCLEOTIDE_MARKER_SEQUENCE"],
      twenty_four_nucleotide_marker_sequence=config["TWENTY_FOUR_NUCLEOTIDE_MARKER_SEQUENCE"],
      kit_5adapter=config["KIT_5ADAPTER"],
      iclip_adapter=config["ICLIP_ADAPTER"],
      hitsclip_3adapter=config["HITSCLIP_3ADAPTER"],
      hitsclip_5adapter=config["HITSCLIP_5ADAPTER"],
      special_ad_1=config["SPECIAL_AD_1"],
      special_ad_2=config["SPECIAL_AD_2"],
      special_ad_3=config["SPECIAL_AD_3"],
      repeat_n_times=2,
      match_read_wildcards=18
    shell:
      """
     cutadapt -a {params.three_prime_adapter_sequence} -g {params.five_prime_adapter_sequence}  \
                -b {params.nineteen_nucleotide_marker_sequence} -b {params.twenty_four_nucleotide_marker_sequence}  \
                -b {params.kit_5adapter} -b {params.iclip_adapter} -b {params.hitsclip_3adapter} -b {params.hitsclip_5adapter}  \
                -b {params.special_ad_1} -b {params.special_ad_2} -b {params.special_ad_3}  \
                -j {threads} -n {params.repeat_n_times} -m {params.match_read_wildcards} -o {output.fastq} {input.fastq} 2> {log}
      """

rule extract_UMIs:
    input:
      fastq="results/process_reads/{sample}_trim.fastq.gz",
    output:
      processed="results/process_reads/{sample}_trim_umi-extr.fastq.gz"
    container:
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
    container:
      "docker://biocontainers/fastxtools:v0.0.14_cv2"
    log:
      "results/logs/collapse_reads_{sample}.log"
    shell:
      """
      zcat {input.fastq} | fastx_collapser | gzip > {output.collapsed} 2> {log}
      """
