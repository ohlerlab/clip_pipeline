rule collapse_reads:
    input:
      fastq=rules.cutadapt.output.fastq,
    output:
      fastq="test_data/{sample}_trim_collapsed.fq.gz",
    conda:
      "../envs/fastx_toolkit.yaml"
    singularity:
      "docker://biocontainers/fastxtools:v0.0.14_cv2"
    log:
      "results/logs/collapse_reads_{sample}.log"
    shell:
      """
      zcat {input.fastq} | fastx_collapser | gzip > {output.fastq} 2> {log}
      """
