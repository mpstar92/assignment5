rule bwa_index:
    input:
        ref = lambda wc: get_final_assembly(wc.sample)
    output:
        ref = "results/assembly/{sample}/spades/contigs.fasta.bwt"
    log:
        "logs/bwa/index_{sample}.log"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa index {input.ref} &> {log}
        """

rule bwa_map:
    input:
        ref = lambda wc: get_final_assembly(wc.sample),
        index = lambda wc: get_final_assembly(wc.sample) + ".bwt",
        fq1 = lambda wc: get_reads(wc)[0],
        fq2 = lambda wc: get_reads(wc)[1]
    output:
        temp("results/mapping/{sample}.sam")
    log:
        "logs/bwa/mem_{sample}.log"
    threads: 4
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.fq1} {input.fq2} > {output} 2> {log}
        """

rule sam_to_bam:
    input:
        "results/mapping/{sample}.sam"
    output:
        "results/mapping/{sample}.bam"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -bS {input} > {output}
        """

rule sort_and_index:
    input:
        bam = "results/mapping/{sample}.bam"
    output:
        sorted = "results/mapping_sorted/{sample}.bam",
        bai = "results/mapping_sorted/{sample}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort {input.bam} -o {output.sorted}
        samtools index {output.sorted}
        """
