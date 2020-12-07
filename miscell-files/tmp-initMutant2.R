

## Fitness = c("0",
        ##             "1.0 + 1.5 * f_1_2",
        ##             "1.3 + 1 * f_1_2",
        ##             "1.2",
        ##             "1.3 + 1.2 * f_3_5",
        ##             "1.4 + 1.3 * f_3"),

gffd <- data.frame(
    Genotype = c("WT", "A", "C", "A, B"), 
    Fitness = c("0",
                "1.3 + 1.5 * f_1_2",
                "1.1 + 0.7*((f_1 + f_2) > 0.3) + f_1_2",
                "1.2 + sqrt(f_1 + f_3 + f_2) - 0.3 * (f_1_2 > 0.5)"),
    stringsAsFactors = FALSE)
    
afd <- allFitnessEffects(genotFitness = gffd, 
                         frequencyDependentFitness = TRUE)
evalAllGenotypes(afd, spPopSizes = rep(10, 4))


## this breaks
gffd <- data.frame(
    Genotype = c("WT", "A", "C", "A, B"), 
    Fitness = c("1",
                "1.3 + 1.2 * n_1_2",
                "1.1 + 0.7*((n_1 + n_2) > 0.3)",
                "1.2 + sqrt(n_1 + n_3 + n_2)"),
    stringsAsFactors = FALSE)
afd <- allFitnessEffects(genotFitness = gffd, 
                         frequencyDependentFitness = TRUE)
evalAllGenotypes(afd, spPopSizes = rep(1:5))


## this works
gffd2 <- data.frame(
    Genotype = c("WT", "A", "B", "C", "A, B"), 
    Fitness = c("1",
                "1.3 + 1.2 * n_1_2",
                "0",
                "1.1 + 0.7*((n_1 + n_2) > 0.3)",
                "1.2 + sqrt(n_1 + n_3 + n_2)"),
    stringsAsFactors = FALSE)
afd2 <- allFitnessEffects(genotFitness = gffd2, 
                         frequencyDependentFitness = TRUE)
evalAllGenotypes(afd2, spPopSizes = rep(1:5))



gffdm2 <- afd2$fitnessLandscape
afd22 <- allFitnessEffects(genotFitness = gffdm2,
                           frequencyDependentFitness = TRUE)
evalAllGenotypes(afd22, spPopSizes = rep(1:5))

## This crashes too
gffdm <- afd2$fitnessLandscape[-3, ]
gffdm
rownames(gffdm) <- 1:4
afd12 <- allFitnessEffects(genotFitness = gffdm,
                           frequencyDependentFitness = TRUE)

evalAllGenotypes(afd12, spPopSizes = rep(1:4))


## But this works. So the problem is in evalRGenotype?
set.seed(1)
oncoSimulIndiv(afd2, onlyCancer = FALSE, finalTime = 500, seed = NULL)
set.seed(1)
oncoSimulIndiv(afd, onlyCancer = FALSE, finalTime = 500, seed = NULL)

## But it blows up here
set.seed(1)
oncoSimulIndiv(afd2, initMutant = "A", onlyCancer = FALSE, finalTime = 500, seed = NULL)
set.seed(1)
oncoSimulIndiv(afd, initMutant = "A", onlyCancer = FALSE, finalTime = 500, seed = NULL)












######################################################################
######################################################################


## This breaks. And it is wrong: It thinks that A, B is B
## In the C++ code, n_2 is taken as 6 and there is no n_1_2










## The problem is that in the output of allFitnessEffects
## $fitnessLandscape_df
## $fitnessLandscape
## have smaller dimension than
## $fitnessLandscapeVariables








## Conjecture:
## It takes the "Genotype" as mapping to n, n_1, n_2, ...

## Genotype is NOT genotype here, it is gene??!!!
g3 <- data.frame(
    Genotype = c("WT", "A", "caracol", "hola"), 
    Fitness = c("1",
                "1.3 + 1.2 * n_1_2",
                "1.1",
                "2 * n_3"),
    stringsAsFactors = FALSE)

fg3 <- allFitnessEffects(genotFitness = g3, 
                         frequencyDependentFitness = TRUE)

evalAllGenotypes(fg3, spPopSizes = c(9, 2, 6, 8))


## spPopSizes are in the same order as given in "Genotype": NO
