# rule download_reference_patch:
#    output:
#        "resources/reference_patch/AC_000007.1.fasta"
#    conda:
#        "../envs/patching.yaml"
 #   shell:
  #      """
   #     mkdir -p resources/reference_patch
    #    wget -O resources/reference_patch/AC_000007.1.fasta "https://www.ncbi.nlm.nih.gov/search/api/sequence/AC_000007.1/?report=fasta"
     #   """

rule patch_with_reference:
    input:
        reference="results/patching/{sample}/selected_reference.fasta",
        assembly=lambda wc: get_patching_input(wc)
    output:
        patched="results/patching/{sample}/ragtag.patch.fasta"
    log:
        "logs/patching/{sample}.log"
    threads: 8
    conda:
        "../envs/patching.yaml"
    shell:
        """
        mkdir -p results/patching/{wildcards.sample}
        ragtag.py patch -t {threads} -o results/patching/{wildcards.sample} \
            {input.reference} {input.assembly} &> {log}
        test -f {output.patched}
        """