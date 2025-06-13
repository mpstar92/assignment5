# workflow/scripts/rename_fasta.py

from Bio import SeqIO

# Statt sys.argv verwendet man das snakemake-Objekt:
input_file  = snakemake.input[0]
output_file = snakemake.output[0]
sample_name = snakemake.wildcards.sample

with open(output_file, "w") as out_f:
    for record in SeqIO.parse(input_file, "fasta"):
        record.id = f"{sample_name}|{record.id}"
        record.description = ""
        SeqIO.write(record, out_f, "fasta")
