import os
import glob

contam_fastas = config.get("contamination_fastas", [])

rule bowtie2_build_contam:
    input:
        fasta=contam_fastas
    output:
        "resources/contamination/contaminants.bt2.index.done"
    threads: 4
    params:
        fasta_list = lambda wildcards, input: ",".join(input.fasta)
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2-build {params.fasta_list} resources/contamination/contaminants
        touch {output}
        """

rule remove_contamination:
    input:
        r1 = "results/trimmed/{sample}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.fastq.gz",
        index = rules.bowtie2_build_contam.output
    output:
        r1 = protected("results/filtered/{sample}_R1.fastq.gz"),
        r2 = protected("results/filtered/{sample}_R2.fastq.gz"),
        sam = temp("results/filtered/{sample}.contam.sam")
    log:
        "logs/contamination/{sample}.log"
    threads: 8
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -x resources/contamination/contaminants \
                -1 {input.r1} -2 {input.r2} \
                -S {output.sam} \
                -p {threads} 2> {log}

        samtools fastq -f 12 -F 256 \
                -1 {output.r1} -2 {output.r2} \
                {output.sam}
        """