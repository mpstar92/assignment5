import os

KRAKEN2_DB_URL = config["kraken2_db"]
KRAKEN2_DB_DIR = "resources/kraken2"

rule download_kraken2_db:
    output:
        tarball = temp("resources/kraken2_db.tar.gz")
    shell:
        """
        wget -O {output.tarball} {KRAKEN2_DB_URL}
        """

rule extract_kraken2_db:
    input:
        tarball = rules.download_kraken2_db.output.tarball
    output:
        directory = directory(KRAKEN2_DB_DIR)
    shell:
        """
        mkdir -p {output.directory}
        tar -xzf {input.tarball} -C {output.directory} --strip-components=1
        """

rule kraken2:
    input:
        r1 = "results/trimmed/{sample}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.fastq.gz",
        db = rules.extract_kraken2_db.output.directory
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
            --db {input.db} {input.r1} {input.r2} &> {log}
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
        multiqc results/kraken2 -o results/kraken2 --force
        """