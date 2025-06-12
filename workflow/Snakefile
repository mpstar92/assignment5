import pandas as pd
import os
import glob
from snakemake.utils import min_version
min_version("7.0")

wrapper_prefix = "https://github.com/snakemake/snakemake-wrappers/raw/master/bio"

configfile: "config/config.yaml"


# Use filtered reads if contamination_fasta is set, else use trimmed reads
contam_files = glob.glob("resources/contamination/FASTA/*.[fF][aAnN]*")

contam_fasta = "resources/contamination/contaminants.fa" if contam_files else None

if contam_files:
    include: "rules/contamination.smk"
def get_reads(wildcards):
    if contam_fasta:
        return [
            f"results/filtered/{wildcards.sample}_R1.fastq.gz",
            f"results/filtered/{wildcards.sample}_R2.fastq.gz"
        ]
    else:
        return [
            f"results/trimmed/{wildcards.sample}_R1.fastq.gz",
            f"results/trimmed/{wildcards.sample}_R2.fastq.gz"
        ]
        
samples = pd.read_csv(config["samples"], sep="\t", index_col="sample")
SAMPLES = samples.index.tolist()

include: "rules/qc.smk"
include: "rules/assembly.smk"
include: "rules/scaffold.smk"
include: "rules/phylogeny.smk"
include: "rules/variability.smk"
include: "rules/screen.smk"


rule all:
    input:
        expand("results/quast/{sample}/report.txt", sample=SAMPLES),
        "results/alignment/mafft_alignment.fasta",
        "results/tree/iqtree.treefile",
        "results/variability/entropy.tsv",
        "results/variability/entropy_plot.png",
        "results/multiqc/multiqc_report.html" 

rule screen:
    input:
        expand("results/kraken2/{sample}.report", sample=SAMPLES),
        "results/kraken2/multiqc_report.html"