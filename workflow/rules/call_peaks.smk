def input_for_paralyzer_ini():
    if config['UMI-BARCODE'] is not False:
      return "results/prepare_aligned/{sample}_sorted_deduplicated.sam"
    else:
      return "results/prepare_aligned/{sample}_aligned_filtered.sam=COLLAPSED"


def input_for_paralyzer():
    if config['UMI-BARCODE'] is not False:
      return "results/prepare_aligned/{sample}_sorted_deduplicated.sam"
    else:
      return "results/prepare_aligned/{sample}_aligned_filtered.sam"


rule paralyzer_ini:
    input:
        tbit=config["REFERENCE_GENOME"] + ".2bit",
        sam=input_for_paralyzer_ini(),
        edit_script="workflow/scripts/editPARalyzerINIfile.pl",
        ini="config/Default_PARalyzer_Parameters.ini"
    output:
        ini=temp("results/paralyzer_params_{sample}.ini")
    singularity:
        "docker://perl:5.36.0-threaded-bullseye"
    log:
        "results/logs/paralyzer_ini_{sample}.log"
    params:
        path="results/call_peaks/{sample}"
    shell:
        "perl {input.edit_script} {input.ini} {params.path} {input.sam} {input.tbit} > {output.ini} 2> {log}"


rule paralyzer:
    input:
        sam=input_for_paralyzer(),
        ini="results/paralyzer_params_{sample}.ini"
    output: 
        multiext("results/call_peaks/{sample}", config["call_peaks"]["OUTPUT_DISTRIBUTIONS_FILE"], 
                                                config["call_peaks"]["OUTPUT_GROUPS_FILE"], 
                                                config["call_peaks"]["OUTPUT_CLUSTERS_FILE"],
                                                config["call_peaks"]["OUTPUT_READS_FILE"])
    singularity:
        "docker://quay.io/biocontainers/paralyzer:1.5--hdfd78af_3"
    log: 
        "results/logs/paralyzer_{sample}.log"
    params:
        path="results/call_peaks/{sample}",
    resources:
        mem_mb_per_cpu='12G'
    shell:
        "PARalyzer {resources.mem_mb_per_cpu} {input.ini} 2> {log}"

