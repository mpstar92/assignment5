rule download_reference_patch:
    output:
        "resources/reference_patch/AC_000007.1.fasta"
    conda:
        "../envs/patching.yaml"
    shell:
        """
        mkdir -p resources/reference_patch
        wget -O resources/reference_patch/AC_000007.1.fasta "https://www.ncbi.nlm.nih.gov/search/api/sequence/AC_000007.1/?report=fasta"
        """

rule patch_with_reference:
    input:
        reference="resources/reference_patch/AC_000007.1.fasta",
        assembly=lambda wildcards: get_patching_input(wildcards.sample)
    output:
        patched="results/patching/{sample}/ragtag.patch.ctg.fasta"
    log:
        "logs/patching/{sample}.log"
    threads: 4
    conda:
        "../envs/patching.yaml"
    shell:
        """
        mkdir -p results/patching/{wildcards.sample}
        ragtag.py patch -t {threads} -o results/patching/{wildcards.sample} {input.reference} {input.assembly} &> {log}
        """