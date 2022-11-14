rule paralyzer_ini:
    input:
        #bam="test_data/{sample}.sorted.bam",
        sam="test_data/{sample}_aligned_filtered.sam",
        edit_file="workflow/scripts/editPARalyzerINIfile.pl",
        ini_file="config/Default_PARalyzer_Parameters.ini"
    output:
        ini_file="results/paralyzer_params_{sample}.ini"
    singularity:
        "docker://perl:5.36.0-threaded-bullseye"
    log:
        "results/logs/paralyzer_ini_{sample}.log"
    params:
        path="test_data/peaks/{sample}"
    shell:
        "perl {input.edit_file} {input.ini_file} {params.path} {input.sam} > {output.ini_file} 2> {log}"


rule paralyzer:
    input:
        #bam="test_data/{sample}.sorted.bam",
        sam="test_data/{sample}_aligned_filtered.sam",
        ini_file="results/paralyzer_params_{sample}.ini"
    output: 
        multiext("test_data/peaks/{sample}", config["call_peaks"]["OUTPUT_DISTRIBUTIONS_FILE"], 
                                            config["call_peaks"]["OUTPUT_GROUPS_FILE"], 
                                            config["call_peaks"]["OUTPUT_CLUSTERS_FILE"],
                                            config["call_peaks"]["OUTPUT_READS_FILE"])
    singularity:
        "docker://quay.io/biocontainers/paralyzer:1.5--hdfd78af_3"
    log: 
        "results/logs/paralyzer_{sample}.log"
    params:
        path="test_data/peaks/{sample}",
    resources:
        mem_gb=12
    shell:
        "PARalyzer {resources.mem_gb}G {input.ini_file} 2> {log}"

