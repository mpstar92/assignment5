rule pilon_polish:
    input:
        bam = "results/mapping_sorted/{sample}.bam",
        bai = "results/mapping_sorted/{sample}.bam.bai",
        ref = "results/assembly/{sample}/spades/contigs.fasta"  
    output:
        polished = "results/assembly/{sample}/polished.fasta"
    log:
        "logs/pilon/{sample}.log"
    threads: 4
    conda:
        "../envs/pilon.yaml"
    shell:
        """
        pilon --genome {input.ref} --frags {input.bam} \
              --output {wildcards.sample}_polished \
              --outdir results/assembly/{wildcards.sample} \
              --threads {threads} > {log} 2>&1
        mv results/assembly/{wildcards.sample}/{wildcards.sample}_polished.fasta {output.polished}
        """