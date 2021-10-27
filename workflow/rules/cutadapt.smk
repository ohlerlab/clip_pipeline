rule cutadapt:
    input:
      fastq=lambda wc : samples.loc[wc.sample, "path"],
    output:
      fastq="test_data/{sample}_trim.fastq.gz",
    conda:
      "../env/cutadapt.yaml" 
    shell:
      """
     cutadapt -a {THREE_PRIME_ADAPTER_SEQUENCE} -g {FIVE_PRIME_ADAPTER_SEQUENCE}  \
                -b {NINETEEN_NUCLEOTIDE_MARKER_SEQUENCE} -b {TWENTY_FOUR_NUCLEOTIDE_MARKER_SEQUENCE}  \
                -b {KIT_5ADAPTER} -b {ICLIP_ADAPTER} -b {HITSCLIP_3ADAPTER} -b {HITSCLIP_5ADAPTER}  \
                -b {SPECIAL_AD_1} -b {SPECIAL_AD_2} -b {SPECIAL_AD_3}  \
                -j {THREADS} -n 2 -m 18 -o {output.fastq} {input.fastq} 
      """

