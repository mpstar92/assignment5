rule ragtag_scaffold:
    input:
        contigs = lambda wc: get_final_assembly(wc.sample),
        ref = config["ragtag_ref"]
    output:
        scaffold = "results/scaffold/{sample}/ragtag.scaffold.fasta"
    log:
        "logs/ragtag/{sample}.log"
    threads: 2
    conda:
        "../envs/scaffold.yaml"
    shell:
        """
        ragtag.py scaffold -t {threads} -o results/scaffold/{wildcards.sample} {input.ref} {input.contigs} &> {log}
        # Overwrite output with only scaffolded contigs (exclude headers with 'unplaced')
        awk '/^>/ {{p = ($0 !~ /unplaced/)}} p' results/scaffold/{wildcards.sample}/ragtag.scaffold.fasta > tmp
        if [ -s tmp ]; then mv tmp {output.scaffold}; else : > {output.scaffold}; rm -f tmp; fi
        """