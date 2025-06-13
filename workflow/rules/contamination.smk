import os
import glob

contam_files = glob.glob("resources/contamination/FASTA/*.[fF][aAnN]*")

rule combine_contaminants:
    input:
        contam_files
    output:
        "resources/contamination/contaminants.fa"
    shell:
        """
        cat {input} > {output}
        """

rule bowtie2_build_contam:
    input:
        fasta = "resources/contamination/contaminants.fa"
    output:
        index = "resources/contamination/contaminants.fa.1.bt2"
    threads: 1
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build {input.fasta} resources/contamination/contaminants.fa
        """

rule remove_contamination:
    input:
        r1 = "results/trimmed/{sample}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.fastq.gz",
        index = "resources/contamination/contaminants.fa.1.bt2"
    output:
        r1 = temp("results/filtered/{sample}_R1.fastq.gz"),
        r2 = temp("results/filtered/{sample}_R2.fastq.gz")
    threads: 4
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -x resources/contamination/contaminants.fa -1 {input.r1} -2 {input.r2} --un-conc-gz results/filtered/{wildcards.sample}_R%.fastq.gz -p {threads} -S /dev/null
        """

