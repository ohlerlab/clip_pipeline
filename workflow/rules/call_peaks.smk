"""
rule paralyzer_parameters:
    input: 
        params_file="Default_PARalyzer_Parameters.ini"
        edit_file="editPARalyzerINIfile.pl"
    output: 
        ini_path="params.ini"
    shell:
        "perl {input.edit_file} {input.params_file} {output.peaks_bed} {input.bam} > {output.ini_path}"
"""

rule paralyzer:
    input:
        bam="test_data/{sample}.sorted.bam",
    output: 
        multiext("results/peaks/{sample}", config["call_peaks"]["OUTPUT_DISTRIBUTIONS_FILE"], 
                                            config["call_peaks"]["OUTPUT_GROUPS_FILE"], 
                                            config["call_peaks"]["OUTPUT_CLUSTERS_FILE"],
                                            config["call_peaks"]["OUTPUT_READS_FILE"])
    params:
        reads="BAM_FILE={input.bam}=COLLAPSED",
        others=config["call_peaks"]
    log: 
        "results/logs/paralyzer_{sample}.log"
    resources:
        mem_mb=4000
    shell:
        "PARalyzer {resources.mem_gb}G {params.reads} {params.others} > {log}"