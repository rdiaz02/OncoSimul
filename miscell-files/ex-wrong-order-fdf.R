df1 <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                  Fitness = c("1",
                              "n_1_3", ## BA
                              "n_2_3", ## CA
                              "n_1_2",
                              "n_1",
                              "n_2"
                              ))
adf1 <- allFitnessEffects(genotFitness = df1,
                          frequencyDependentFitness = TRUE)

(adf1)
evalAllGenotypes(adf1, spPopSizes = 1:6) ## Breaks in old: n_2_3


df2 <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                  Fitness = c("1",
                              "n_1_3", ## you want BA
                              "n_3", 
                              "n_1_2",
                              "n_1",
                              "n_2"
                              ))
adf2 <- allFitnessEffects(genotFitness = df2,
                          frequencyDependentFitness = TRUE)

(adf2)
evalAllGenotypes(adf2, spPopSizes = 1:6) 
## Wrong: gives for B fitness of using CA, not BA
