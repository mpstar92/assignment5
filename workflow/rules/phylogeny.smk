rule rename_scaffolds:
    input:
        "results/scaffold/{sample}/ragtag.scaffold.fasta"
    output:
        "results/scaffold/{sample}/renamed.fasta"
    conda:
        "../envs/scaffold.yaml"
    script:
        "../scripts/rename_fasta.py"

rule mafft:
    input:
        lambda wildcards: expand("results/scaffold/{sample}/renamed.fasta", sample=SAMPLES)
    output:
        "results/alignment/mafft_alignment.fasta"
    log:
        "logs/mafft/alignment.log"
    threads: 12
    conda:
        "../envs/phylo.yaml"
    shell:
        """
        cat {input} > results/alignment/tmp_combined.fasta
        mafft --thread {threads} --auto results/alignment/tmp_combined.fasta > {output} 2> {log}
        """

rule rtree:
    input:
        alignment = "results/alignment/mafft_alignment.fasta"
    output:
        tree = "results/tree/phylo_tree.nwk"
    threads: 8
    conda:
        "../envs/phylo_r.yaml"
    script:
        "../scripts/tree_builder.R"


#rule iqtree:
 #   input:
  #      "results/alignment/mafft_alignment.fasta"
   # output:
    #    tree = "results/tree/iqtree.treefile",
    #log:
    #    "logs/iqtree/iqtree.log"
    #threads: 8
    #conda:
    #    "../envs/phylo.yaml"
    #shell:
    #    """
    #    iqtree2 -s {input} -nt {threads} {config[iqtree][extra]} -pre results/tree/iqtree &> {log}
    #    """