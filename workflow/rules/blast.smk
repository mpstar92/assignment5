rule makeblastdb_candidates:
    input:
        fasta=config["candidate_refs"]
    output:
        db_done="reference/candidate_refs.nsq"
    params:
        db_name="reference/candidate_refs"
    conda:
        "../envs/blast.yaml"
    shell:
        """
        mkdir -p reference
        makeblastdb -in {input.fasta} -dbtype nucl -out {params.db_name}
        """

rule blastn_viral:
    input:
        fasta=lambda wc: get_final_assembly(wc.sample),
        db_done="reference/candidate_refs.nsq"
    output:
        txt="results/blast/{sample}_viral.txt"
    params:
        db="reference/candidate_refs"
    threads: 8
    conda:
        "../envs/blast.yaml"
    shell:
        """
        blastn -query {input.fasta} \
               -db {params.db} \
               -out {output.txt} \
               -outfmt 6 \
               -max_target_seqs 5 \
               -evalue 1e-10 \
               -num_threads {threads}
        """

rule blast_candidates:
    input:
        assembly=lambda wc: get_patching_input(wc),
        db_done="reference/candidate_refs.nsq"
    output:
        tsv="results/blast/{sample}_candidates.tsv"
    params:
        db="reference/candidate_refs"
    threads: 8
    conda:
        "../envs/blast.yaml"
    shell:
        """
        blastn -query {input.assembly} -db {params.db} \
            -out {output.tsv} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
            -max_target_seqs 10 -evalue 1e-5 -num_threads {threads}
        """

rule select_reference:
    input:
        blast="results/blast/{sample}_candidates.tsv",
        fasta=config["candidate_refs"]
    output:
        selected="results/patching/{sample}/selected_reference.fasta"
    conda:
        "../envs/helpers.yaml"
    run:
        import sys, pandas as pd
        df = pd.read_csv(input.blast, sep="\t", header=None, \
                         names=["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"])
        if df.empty:
            sys.exit(f"No BLAST hits for sample {wildcards.sample}")
        best = df.loc[df.bitscore.idxmax(), "sseqid"]
        shell("mkdir -p results/patching/{wildcards.sample}")
        shell("seqtk subseq {input.fasta} <(echo {best}) > {output.selected}")