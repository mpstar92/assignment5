rule nextpolish:
    input:
        fasta = "results/assembly/{sample}/spades/contigs.fasta",
        fq1 = "results/filtered/{sample}_R1.fastq.gz",
        fq2 = "results/filtered/{sample}_R2.fastq.gz"
    output:
        polished = "results/assembly/{sample}/spades/contigs.nextpolish.fa"
    log:
        "logs/nextpolish/{sample}.log"
    threads: 8
    conda:
        "../envs/nextpolish.yaml"
    shell:
        """
        rm -rf results/nextpolish/{wildcards.sample}
        mkdir -p results/nextpolish/{wildcards.sample}
        cp {input.fasta} results/nextpolish/{wildcards.sample}/genome
        cd results/nextpolish/{wildcards.sample}
        bwa index genome
        bwa mem -t {threads} genome ${{PWD}}/../../../{input.fq1} ${{PWD}}/../../../{input.fq2} | samtools sort -@ {threads} -o {wildcards.sample}.bam
        samtools index {wildcards.sample}.bam
        echo -e "${{PWD}}/../../../{input.fq1}\n${{PWD}}/../../../{input.fq2}" > {wildcards.sample}.sgs.fofn
        echo -e "[General]\njob_type = local\njob_prefix = nextpolish\ntask = best\nrewrite = yes\nrerun = 3\nparallel_jobs = {threads}\nmultithread_jobs = {threads}\nmin_read_len = 40\npolish_options = -p {threads}\n\n[sgs]\ngenome = genome\nsgs_fofn = {wildcards.sample}.sgs.fofn\nbam = {wildcards.sample}.bam\n" > {wildcards.sample}.nextpolish.cfg
        mkdir -p ../../logs/nextpolish
        nextPolish {wildcards.sample}.nextpolish.cfg > ../../logs/nextpolish/{wildcards.sample}.log 2>&1
        mv genome.nextpolish.fasta ../../assembly/{wildcards.sample}/spades/contigs.nextpolish.fa
        cd ../../../..
        rm -rf results/nextpolish/{wildcards.sample}
        """