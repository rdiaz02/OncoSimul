library(OncoSimulR)

r1 <- rfitness(1)

N <- 100
initS <- 100
finalTime <- 1000
s <- 0.1
ngenes <- 5
ni <- rep(s, ngenes)
names(ni) <- letters[1:ngenes]
fe <- allFitnessEffects(noIntGenes = ni, drvNames = names(ni))
op <- oncoSimulPop(N,
                   fe,
                   model = "McFL",
                   onlyCancer = FALSE,
                   initSize = initS,
                   finalTime = finalTime)







## r1 <- data.frame(Genotype = c("WT", "A"), Fitness = c(1, 1 + s),
##                  stringsAsFactors = FALSE)

## fe <- allFitnessEffects(genotFitness = data.frame(r1))
