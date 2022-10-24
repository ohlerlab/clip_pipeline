rule cutadapt:
    input:
      fastq=lambda wc : samples.loc[wc.sample, "path"],
    output:
      fastq="test_data/{sample}_trim.fastq.gz",
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

