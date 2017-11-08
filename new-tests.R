df1 <- data.frame(Genotype = c("A", "B, C"), Fitness = c(1.3, 2),
                  stringsAsFactors = FALSE)
expect_warning(OncoSimulR:::allGenotypes_to_matrix(df1),
               "No WT genotype. Setting its fitness to 1.", fixed = TRUE)

df2 <- data.frame(Genotype = c("WT", "A", "B, C"), Fitness = c(5, 1.3, 2))
expect_warning(OncoSimulR:::allGenotypes_to_matrix(df2),
               "First column of genotype fitness is a factor.",
               fixed = TRUE)




df1 <- data.frame(Genotype = c("A", "B, C"), Fitness = c(1.3, 2),
                  stringsAsFactors = FALSE)
expect_warning(OncoSimulR:::to_genotFitness_std(df1),
               "No WT genotype. Setting its fitness to 1.", fixed = TRUE)


df2 <- data.frame(Genotype = c("WT", "A", "B, C"), Fitness = c(5, 1.3, 2))
expect_warning(OncoSimulR:::to_genotFitness_std(df2),
               "First column of genotype fitness is a factor.",
               fixed = TRUE)

df3 <- data.frame(Genotype = c("WT", "A", "B, C"), Fitness = c(1.3, 2, 0),
                  stringsAsFactors = FALSE)

m3 <- rbind(c(0, 0, 0, 1.3), c(1, 0, 0, 2.0))
colnames(m3) <- c("A", "B", "C", "Fitness")
m4 <- rbind(c(0, 0, 0, 1.3), c(1, 0, 0, 2.0), c(0, 1, 1, 0))
colnames(m4) <- c("A", "B", "C", "Fitness")
expect_equal(OncoSimulR:::to_genotFitness_std(df3), m3)
expect_equal(OncoSimulR:::to_genotFitness_std(df3, simplify = FALSE), m4)

for(i in 1:10) {
    rxx <- rfitness(7)
    expect_equal(
        allFitnessEffects(genotFitness = rxx)$fitnessLandscape,
        OncoSimulR:::to_genotFitness_std(rxx))
}


for(i in 1:10) {
    ng <- 7
    rxx <- rfitness(ng)
    cnn <-
        replicate(ng,
                  paste(sample(c(LETTERS,
                                 letters,
                                 0:9), 5, replace = TRUE),
                        collapse = ""))
    if(any(duplicated(cnn))) cnn <- LETTERS[1:ng]
    colnames(rxx)[1:ng] <- cnn
    fex <- allFitnessEffects(genotFitness = rxx)
    gn <- OncoSimulR:::allNamedGenes(fex)
    m1x <- data.frame(Gene = cnn, GeneNumID = 1:ng,
                      stringsAsFactors = FALSE)
    expect_identical(m1x, gn)
}

## FIXME: to do
## taken an rT, convert to fitness landscape, and verify we get
## same fitnesses

## this should work!
allFitnessEffects(genotFitness = rxx, drvNames = LETTERS[1:4])
