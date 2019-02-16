## Combine the genotype-fitness mapping with other specifications.


## Allow this/document this in OncoSimulR?
ee <- OncoSimulR:::from_genotype_fitness(fi1)
numnoint <- 7
noint <- runif(numnoint, -.01, .01)
names(noint) <- paste0("n", seq.int(numnoint))

feee <- allFitnessEffects(epistasis = ee, noIntGenes = noint)
