rule collapse_reads:
    input:
      fastq=rules.cutadapt.output.fastq,
    output:
      fastq="test_data/{sample}_trim_collapsed.fq.gz",
    conda:
      "../envs/fastx_toolkit.yaml"
    shell:
      """
      zcat {input.fastq} | fastx_collapser | gzip > {output.fastq}
      """
