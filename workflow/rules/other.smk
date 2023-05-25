rule convert_sam2bam:
    input:
      sam="results/prepare_aligned/{sample}_aligned_filtered.sam",
    output:
      bam="results/prepare_aligned/{sample}_sorted.bam",
      bai="results/prepare_aligned/{sample}_sorted.bam.bai",
    conda:
      "../envs/samtools.yaml"
    container:
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
    container:
      "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    log:
      "results/logs/samtools_b2s_{sample}"
    shell:
      """
      samtools index {input.bam_dedup} 2> {log} &&
      samtools view -h -o {output.sam_dedup} {input.bam_dedup} 2> {log}
      """


rule filter_unmapped:
    input: 
      sam="results/align_reads/{sample}_aligned.sam",
      main_chr_bed=config["REFERENCE_GENOME"] + ".main_chr.bed"
    output:
      sam_filtered="results/prepare_aligned/{sample}_aligned_filtered.sam"
    conda:
      "../envs/samtools.yaml"
    container:
      "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    log:
      "results/logs/samtools_filter_unmapped_{sample}"
    shell:
      """
      samtools view -h -F 4 -L {input.main_chr_bed} {input.sam} > {output.sam_filtered}
      """

rule main_chromosomes:
    input:
      fa=config["REFERENCE_GENOME"]
    output:
      fai=config["REFERENCE_GENOME"] + ".fai",
      main_chr_bed=config["REFERENCE_GENOME"] + ".main_chr.bed"
    conda:
      "../envs/samtools.yaml"
    container:
      "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    log:
      "results/logs/main_chr_bed.log"
    shell:
      """
      samtools faidx {input.fa} &&
      cat {output.fai} | grep chr | grep -v chrM | cut -f 1,2 | awk '{{print $1,"1",$2}}' OFS="\t" > {output.main_chr_bed}
      """
 

rule UMI_deduplicate:
    input:
      bam="results/prepare_aligned/{sample}_sorted.bam"
    output:
      bam_dedup="results/prepare_aligned/{sample}_sorted_deduplicated.bam"
    container:
      "docker://quay.io/biocontainers/umi_tools:1.1.2--py310h1425a21_1"
    log:
      "results/logs/umi_tools_deduplicate_{sample}"
    shell:
      """
      umi_tools dedup -I {input.bam} --output-stats=deduplicated -S {output.bam_dedup} 2> {log}
      """

rule index_bam:
    input:
      bam_dedup="results/prepare_aligned/{sample}_sorted_deduplicated.bam"
    output:
      idx="results/prepare_aligned/{sample}_sorted_deduplicated.bam.bai"
    container:
      "docker://quay.io/biocontainers/samtools:1.9--h91753b0_8"
    shell:
      """
      samtools index {input.bam_dedup} 
      """



rule faToTwoBit_fa:
    input:
      fa=config["REFERENCE_GENOME"]
    output:
      tbit=config["REFERENCE_GENOME"] + ".2bit"
    log:
      "results/logs/fa_to_2bit.log"
    container:
      "docker://quay.io/biocontainers/ucsc-fatotwobit:377--ha8a8165_5"
    shell:
      """
      faToTwoBit {input.fa} {output.tbit} 2> {log}
      """
