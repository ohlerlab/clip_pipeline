rule segemehl_idx:
    input:
      fasta=FA,
    output:
      sege_idx=SEGE_IDX,
    conda:
      "../envs/segemehl.yaml"
    singularity:
      "docker://quay.io/biocontainers/segemehl:0.3.1--h39379e4_3"
    shell:
      """
      segemehl.x -x {output.sege_idx} -d {input.fasta}
      """


rule segemehl:
    input:
      fastq="test_data/{sample}_trim_collapsed.fq.gz",
      sege_idx=SEGE_IDX,
    output:
      sam=temp("test_data/{sample}.sam"),
    conda:
      "../envs/segemehl.yaml"     
    singularity:
      "docker://quay.io/biocontainers/segemehl:0.3.1--h39379e4_3"
    shell:
      """
      segemehl.x -S -D 2 -M 1 --briefcigar -t 4 -i {input.sege_idx} -d {FA} -q {input.fastq}  > {output.sam} 
      """
