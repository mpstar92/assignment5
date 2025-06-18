library(ape)

alignment <- read.dna(snakemake@input[[1]], format = "fasta")
dist_matrix <- dist.dna(alignment, model = "raw", pairwise.deletion = TRUE)
tree <- njs(dist_matrix)
write.tree(tree, file = snakemake@output[[1]])