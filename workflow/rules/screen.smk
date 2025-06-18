KRAKEN_DB = config["kraken2_db"]

rule kraken2:
    input:
        r1 = "results/trimmed/{sample}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.fastq.gz"
    output:
        report = "results/kraken2/{sample}.report",
        output = "results/kraken2/{sample}.kraken"
    log:
        "logs/kraken2/{sample}.log"
    conda:
        "../envs/kraken2.yaml"
    shell:
        """
        kraken2 --paired --report {output.report} --output {output.output} \
            --db {KRAKEN_DB} {input.r1} {input.r2} &> {log}
        """

rule multiqc_kraken2:
    input:
        expand("results/kraken2/{sample}.report", sample=SAMPLES)
    output:
        "results/kraken2/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        multiqc results/kraken2 -o results/kraken2
        """