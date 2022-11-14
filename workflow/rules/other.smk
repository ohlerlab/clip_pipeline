rule convert_sam2bam:
    input:
      sam="test_data/{sample}_aligned_filtered.sam",
    output:
      bam="test_data/{sample}_sorted.bam",
      bai="test_data/{sample}_sorted.bam.bai",
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


rule filter_unmapped:
    input: 
      sam="test_data/{sample}_aligned.sam"
    output:
      sam_filtered="test_data/{sample}_aligned_filtered.sam"
    conda:
      "../envs/samtools.yaml"
    singularity:
      "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    shell:
      """
      samtools view -h -F 4 {input.sam} > {output.sam_filtered}
      """


rule UMI_deduplicate:
      input:
        bam="test_data/{sample}_sorted.bam"
      output:
        bam_dedup="test_data/{sample}_sorted_deduplicated.bam"
      singularity:
        "docker://quay.io/biocontainers/umi_tools:1.1.2--py310h1425a21_1"
      shell:
        """
        umi_tools dedup -I {input.bam} --output-stats=deduplicated -S {output.bam_dedup}
        """
