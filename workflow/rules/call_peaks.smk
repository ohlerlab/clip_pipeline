rule paralyzer_ini:
    input:
        tbit=config["REFERENCE_GENOME"] + ".2bit",
        sam=lambda wildcards: "results/prepare_aligned/{sample}_sorted_deduplicated.sam".format(sample=wildcards.sample) if config['UMI-BARCODE'] is not False else "results/prepare_aligned/{sample}_aligned_filtered.sam=COLLAPSED".format(sample=wildcards.sample),
        edit_script="workflow/scripts/editPARalyzerINIfile.pl",
        ini="config/Default_PARalyzer_Parameters.ini"
    output:
        ini=temp("results/paralyzer_params_{sample}.ini")
    container:
        "docker://perl:5.36.0-threaded-bullseye"
    log:
        "results/logs/paralyzer_ini_{sample}.log"
    params:
        path="results/call_peaks/{sample}"
    shell:
        "perl {input.edit_script} {input.ini} {params.path} {input.sam} {input.tbit} > {output.ini} 2> {log}"


rule paralyzer:
    input:
        ini="results/paralyzer_params_{sample}.ini"
    output: 
        multiext("results/call_peaks/{sample}", config["call_peaks"]["OUTPUT_DISTRIBUTIONS_FILE"], 
                                                config["call_peaks"]["OUTPUT_GROUPS_FILE"], 
                                                config["call_peaks"]["OUTPUT_CLUSTERS_FILE"],
                                                config["call_peaks"]["OUTPUT_READS_FILE"])
    container:
        "docker://quay.io/biocontainers/paralyzer:1.5--hdfd78af_3"
    log: 
        "results/logs/paralyzer_{sample}.log"
    resources:
        mem_mb_per_cpu='12G'
    shell:
        "PARalyzer {resources.mem_mb_per_cpu} {input.ini} 2> {log}"

