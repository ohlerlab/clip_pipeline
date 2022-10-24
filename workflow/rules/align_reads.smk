rule segemehl_idx:
    input:
      fasta=config["REFERENCE_GENOME"],
    output:
      sege_idx=SEGE_IDX,
    conda:
      "../envs/segemehl.yaml"
    singularity:
      "docker://quay.io/biocontainers/segemehl:0.3.1--h39379e4_3"
    log:
      "results/logs/segemehl_idx.log"
    shell:
      """
      segemehl.x -x {output.sege_idx} -d {input.fasta} 2> {log}
      """


rule segemehl:
    input:
      fastq="test_data/{sample}_trim_collapsed.fq.gz",
      sege_idx=SEGE_IDX,
      ref=config["REFERENCE_GENOME"]
    output:
      sam="test_data/{sample}_aligned.sam",
    conda:
      "../envs/segemehl.yaml"     
    singularity:
      "docker://quay.io/biocontainers/segemehl:0.3.1--h39379e4_3"
    log:
      "results/logs/segemehl_{sample}.log"
    threads: 8
    params:
      max_occurences=1,
      differences=2
    shell:
      """
      segemehl.x -S -D {params.differences} -M {params.max_occurences} --briefcigar -t {threads} -i {input.sege_idx} -d {input.ref} -q {input.fastq}  > {output.sam} 2> {log}
      """
