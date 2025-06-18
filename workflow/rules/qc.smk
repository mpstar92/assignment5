# rules/qc.smk

rule fastqc_raw:
    input:
        fq = lambda wc: samples.at[wc.sample, f"fq{wc.idx}"]
    output:
        html = temp("results/fastqc_raw/{sample}_{idx}_fastqc.html"),
        zip = temp("results/fastqc_raw/{sample}_{idx}_fastqc.zip")
    log:
        "logs/fastqc_raw/{sample}_{idx}.log"
    threads: 2
    conda:
        "../envs/qc.yaml"
    wrapper:
        "v6.2.0/bio/fastqc"

rule fastp:
    input:
        fq1 = lambda wc: samples.at[wc.sample, "fq1"],
        fq2 = lambda wc: samples.at[wc.sample, "fq2"]
    output:
        fq1 = protected("results/trimmed/{sample}_R1.fastq.gz"),
        fq2 = protected("results/trimmed/{sample}_R2.fastq.gz"),
        html = "results/qc/fastp/{sample}.html",
        json = "results/qc/fastp/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    threads: 4
    conda:
        "../envs/qc.yaml"
    shell:
        """
        fastp -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} \
              --html {output.html} --json {output.json} -w {threads} {config[fastp][extra]} > {log} 2>&1
        """

rule fastqc_trimmed:
    input:
        fq = "results/trimmed/{sample}_R{idx}.fastq.gz"
    output:
        html = temp("results/fastqc_trimmed/{sample}_{idx}_fastqc.html"),
        zip = temp("results/fastqc_trimmed/{sample}_{idx}_fastqc.zip")
    log:
        "logs/fastqc_trimmed/{sample}_{idx}.log"
    threads: 2
    conda:
        "../envs/qc.yaml"
    wrapper:
        "v6.2.0/bio/fastqc"


rule multiqc:
    input:
        expand("results/fastqc_raw/{sample}_{idx}_fastqc.zip", sample=samples.index, idx=["1", "2"]) +
        expand("results/fastqc_trimmed/{sample}_{idx}_fastqc.zip", sample=samples.index, idx=["1", "2"]) +
        expand("results/qc/fastp/{sample}.json", sample=samples.index) +
        expand("results/quast/{sample}/report.txt", sample=samples.index)
    output:
        html = "results/multiqc/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log"
    threads: 1
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc results/fastqc_raw results/fastqc_trimmed results/qc/fastp results/quast -o results/multiqc > {log} 2>&1
        """