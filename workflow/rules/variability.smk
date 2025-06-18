rule compute_entropy:
    input:
        "results/alignment/mafft_alignment.fasta"
    output:
        "results/variability/entropy.tsv"
    log:
        "logs/variability/entropy.log"
    conda:
        "../envs/variability.yaml"
    script:
        "../scripts/variability.py"

rule plot_entropy:
    input:
        "results/variability/entropy.tsv"
    output:
        "results/variability/entropy_plot.png"
    log:
        "logs/variability/plot.log"
    conda:
        "../envs/variability.yaml"
    shell:
        """
        python -c "import pandas as pd; import matplotlib.pyplot as plt; \
        df = pd.read_csv('{input}', sep='\t'); \
        plt.plot(df['start'], df['entropy']); \
        plt.xlabel('Window Start'); \
        plt.ylabel('Avg. Entropy'); \
        plt.title('Sequence Variability'); \
        plt.savefig('{output}')" &> {log}
        """