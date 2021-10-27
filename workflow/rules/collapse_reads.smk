rule collapse_reads:
    input:
      fastq="test_data/{sample}_trim.fq.gz",
    output:
      fastq="test_data/{sample}_trim_collapsed.fq.gz",
    conda:
      "../envs/fastx_toolkit.yaml"
    shell:
      """
      zcat {input.fastq} | fastx_collapser | gzip > {output.fastq}
      """
