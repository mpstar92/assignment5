rule download_ref_viruses_rep_genomes:
    output:
        "reference/blast_db/ref_viruses_rep_genomes.nsq"
    shell:
        """
        mkdir -p reference/blast_db
        cd reference/blast_db
        wget -N https://ftp.ncbi.nlm.nih.gov/blast/db/ref_viruses_rep_genomes.tar.gz
        tar -xzvf ref_viruses_rep_genomes.tar.gz
        """

rule blastn_viral:
    input:
        fasta="results/assembly/{sample}/polished.fasta",
        db_ready="reference/blast_db/ref_viruses_rep_genomes.nsq"
    output:
        txt="results/blast/{sample}_viral.txt"
    params:
        db_path="reference/blast_db/ref_viruses_rep_genomes"
    threads: 4
    conda:
        "../envs/blast.yaml"
    shell:
        """
        blastn -query {input.fasta} \
               -db {params.db_path} \
               -out {output.txt} \
               -outfmt 6 \
               -max_target_seqs 5 \
               -evalue 1e-10 \
               -num_threads {threads}
        """