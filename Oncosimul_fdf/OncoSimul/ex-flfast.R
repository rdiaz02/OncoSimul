

df1 <- data.frame(Genotype = c("WT", "A", "B, C"), Fitness = c(5, 1.3, 2),
                  stringsAsFactors = FALSE)
OncoSimulR:::allGenotypes_to_matrix(df1)

set.seed(4)
(r4 <- rfitness(4))
