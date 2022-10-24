rule convert_sam2bam:
    input:
      sam="test_data/{sample}_aligned.sam",
    output:
      bam="test_data/{sample}.sorted.bam",
      bai="test_data/{sample}.sorted.bam.bai",
    conda:
      "../envs/samtools.yaml"
    singularity:
      "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    log:
      "results/logs/samtools_{sample}"
    shell:
      """
      samtools view -b {input.sam} | samtools sort -o {output.bam} - &&
      samtools index {output.bam} 2> {log}
      """