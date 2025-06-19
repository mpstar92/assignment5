
rule spades:
    input:
        fq1 = lambda wc: get_reads(wc)[0],
        fq2 = lambda wc: get_reads(wc)[1]
    output:
        fasta = "results/assembly/{sample}/spades/contigs.fasta"
    log:
        "logs/assembly/{sample}.log"
    threads: 8
    conda:
        "../envs/assembly.yaml"
    #shell:
     #   """
      #  spades.py --phred-offset 33 -1 {input.fq1} -2 {input.fq2} -o results/assembly/{wildcards.sample}/spades -t {threads} > {log} 2>&1
       # """
    shell:
        """
        spades.py --phred-offset 33 -1 {input.fq1} -2 {input.fq2} -o results/assembly/{wildcards.sample}/spades -t {threads} {config[spades][extra]} > {log} 2>&1
        """

rule quast:
    input:
        fasta = lambda wildcards: get_final_assembly(wildcards.sample)
    output:
        "results/quast/{sample}/report.txt"
    log:
        "logs/quast/{sample}.log"
    threads: 4
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        quast --threads {threads} -o results/quast/{wildcards.sample} {input.fasta} > {log} 2>&1
        """