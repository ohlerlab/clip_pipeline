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


rule omniclip_db:
    input:
        gene_anot_gff=config["GFF"]
    output:
        db_sql="results/omniclip/{sample}_db.db"
    container:
        "docker://zavolab/omniclip:master_83e44f8a044227afb3d457ee33a8da7e07964d91_9"
    log:
        "results/logs/omniclip/omniclip_{sample}_db.log"
    resources:
    shell:
        """
        omniCLIP generateDB --gff-file {input.gene_anot_gff} --db-file {output.db_sql} 2> {log}
        """


rule omniclip_parse_bg:
    input:
        db_sql="results/omniclip/{sample}_db.db",
        background_bam=config["BACKGROUND"],
        ref=config["GENOME_DIR"]
    output:
        bg_dat="results/omniclip/{sample}_bg.dat"
    container:
        "docker://zavolab/omniclip:master_83e44f8a044227afb3d457ee33a8da7e07964d91_9"
    log:
        "results/logs/omniclip/omniclip_{sample}_parse_bg.log"
    resources:
    shell:
        """
        omniCLIP parsingBG --db-file {input.db_sql} --genome-dir {input.ref} \
        --bg-files {input.background_bam} --out-file {output.bg_dat} 2> {log}
        """


rule omniclip_parse_clip:
    input:
        clip_bam="results/prepare_aligned/{sample}_sorted_deduplicated.bam",
        db_sql="results/omniclip/{sample}_db.db",
        ref=config["GENOME_DIR"]
    output:
        clip_dat="results/omniclip/{sample}_clip.dat"
    container:
        "docker://zavolab/omniclip:master_83e44f8a044227afb3d457ee33a8da7e07964d91_9"
    log:
        "results/logs/omniclip/omniclip_{sample}_parse_clip.log"
    resources:
        mem_mb_per_cpu='12G'
    shell:
        """
        omniCLIP parsingCLIP --db-file {input.db_sql} --genome-dir {input.ref} \
        --clip-files {input.clip_bam} --out-file {output.clip_dat} 2> {log}
        """


rule omniclip_run:
    input:
        db_sql="results/omniclip/{sample}_db.db",
        bg_dat="results/omniclip/{sample}_bg.dat",
        clip_dat="results/omniclip/{sample}_clip.dat"
    output:
        outdir=directory("results/omniclip/output_{sample}/")
    container:
        "docker://zavolab/omniclip:master_83e44f8a044227afb3d457ee33a8da7e07964d91_9"
    log:
        "results/logs/omniclip/omniclip_{sample}_run.log"
    resources:
        mem_mb_per_cpu='12G'
    shell:
        """
        mkdir -p {output.outdir} &&
        omniCLIP run_omniCLIP --db-file {input.db_sql} --bg-dat {input.bg_dat} \
        --clip-dat {input.clip_dat} --out-dir {output.outdir} 2> {log}
        """
