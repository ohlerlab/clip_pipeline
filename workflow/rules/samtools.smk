rule samtools:
    input:
      sam="test_data/{sample}.sam",
    output:
      bam="test_data/{sample}.sorted.bam",
      bai = "test_data/{sample}.sorted.bam.bai",
    conda:
      "../envs/samtools.yaml"
    shell:
      """
      samtools view -b {input.sam} | samtools sort -o {output.bam} - &&
      samtools index {output.bam}
      """
