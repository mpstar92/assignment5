rule ragtag_scaffold:
    input:
        reference="results/patching/{sample}/selected_reference.fasta",
        assembly="results/patching/{sample}/ragtag.patch.fasta"
    output:
        scaffold="results/scaffold/{sample}/ragtag.scaffold.fasta"
    log:
        "logs/ragtag/{sample}.log"
    threads: 8
    conda:
        "../envs/patching.yaml"
    shell:
        """
        mkdir -p results/scaffold/{wildcards.sample}
        ragtag.py scaffold -t {threads} -o results/scaffold/{wildcards.sample} \
            {input.reference} {input.assembly} &> {log}
        """