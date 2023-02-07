rule convert_sam2bam:
    input:
      sam="results/prepare_aligned/{sample}_aligned_filtered.sam",
    output:
      bam="results/prepare_aligned/{sample}_sorted.bam",
      bai="results/prepare_aligned/{sample}_sorted.bam.bai",
    conda:
      "../envs/samtools.yaml"
    singularity:
      "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    log:
      "results/logs/samtools_s2b_{sample}"
    shell:
      """
      samtools view -S -b {input.sam} | samtools sort -o {output.bam} - &&
      samtools index {output.bam} 2> {log}
      """


rule convert_bam2sam:
    input:
      bam_dedup="results/prepare_aligned/{sample}_sorted_deduplicated.bam"
    output:
      sam_dedup="results/prepare_aligned/{sample}_sorted_deduplicated.sam"
    conda:
      "../envs/samtools.yaml"
    singularity:
      "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    log:
      "results/logs/samtools_b2s_{sample}"
    shell:
      """
      samtools view -h -o {output.sam_dedup} {input.bam_dedup} 2> {log}
      """


rule filter_unmapped:
    input: 
      sam="results/align_reads/{sample}_aligned.sam"
    output:
      sam_filtered="results/prepare_aligned/{sample}_aligned_filtered.sam"
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
        bam="results/prepare_aligned/{sample}_sorted.bam"
      output:
        bam_dedup="results/prepare_aligned/{sample}_sorted_deduplicated.bam"
      singularity:
        "docker://quay.io/biocontainers/umi_tools:1.1.2--py310h1425a21_1"
      shell:
        """
        umi_tools dedup -I {input.bam} --output-stats=deduplicated -S {output.bam_dedup}
        """
