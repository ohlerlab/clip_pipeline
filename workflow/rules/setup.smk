
###############
# Setup rules #
###############

# Below are some rules that demonstrate the functionality offered by snakemake, including

# - working with a sample-sheet
# - wildcards
# - logfiles
# - the "run:" and  "shell:" directives
# - conda and singularity integration


rule link_input_files:
    input:
        # the input files
        samples.path
    output:
        # the output files
        expand("resources/input/{sample}.bed.gz", sample=samples.index)
    run:
        # the "run" directive can execute python code
        for sample in samples.index:
            infile = samples.loc[sample,'path']
            outfile = 'resources/input/{}.bed.gz'.format(sample)
            # the "shell()" function can call shell commands within python
            shell('ln -s -r {infile} {outfile}'.format(infile=infile, outfile=outfile))

rule index_reference:
    input:
        config['reference']
    output:
        ref_fa="resources/reference/reference.fa",
        ref_fai="resources/reference/reference.fa.fai"
    log:
        "logs/index_reference.log"
    conda:
        # a path to a conda environment to use
        # the path is relative to the file the rule is defined in.
        "../envs/samtools.yaml"
    singularity:
        # alternatively, a docker container:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        # the "shell" directive allows to run shell commands
        "("
        "ln -s -r {input} {output[ref_fa]} && "
        "samtools faidx {output[ref_fa]} "
        ") &> {log}"

rule sort_bed:
    input:
        "resources/input/{sample}.bed.gz"
    output:
        # we can make certain output files temporary with "temp(...)"
        temp("resources/input/{sample}.srt.bed")
    log:
        'logs/sort_bed/{sample}.log'
    shell:
        "("
        "zcat {input} | sort -k1,1 -k2,2n -k3,3n > {output} "
        ") &> {log}"
    
rule intersect_bed:
    input:
        expand(rules.sort_bed.output, sample=samples.index.tolist())
    output:
        "results/intersect_bed/intersected.bed.gz"
    params:
        # we can set additional parameters with "params:"
        samples_list = lambda wc, input: ' '.join(input)
    log:
        'logs/intersect_bed.log'
    conda:
        '../envs/bedtools.yaml'
    singularity:
        "docker://biocontainers/bedtools:v2.28.0_cv2"
    shell:
        "("
        "multiIntersectBed -i {params[samples_list]} | gzip > {output} "
        ") &> {log}"

