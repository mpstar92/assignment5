rule generate_consensus:
    input:
        bam="results/mapping_sorted/{sample}.bam",
        bai="results/mapping_sorted/{sample}.bam.bai",
        ref="results/patching/{sample}/ragtag.patch.fasta"
    output:
        consensus="results/consensus/{sample}.consensus.fasta"
    log:
        "logs/consensus/{sample}.log"
    threads: 8
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        mkdir -p results/consensus
        bcftools mpileup -f {input.ref} {input.bam} | \
          bcftools call -mv -Oz -o results/consensus/{wildcards.sample}.vcf.gz
        tabix -p vcf results/consensus/{wildcards.sample}.vcf.gz
        bcftools consensus -f {input.ref} results/consensus/{wildcards.sample}.vcf.gz > {output.consensus} 2> {log}
        """