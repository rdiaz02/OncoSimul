library(OncoSimulR)

set.seed(456)
nd <- 70  
np <- 2000 
s <- 0.1  
sp <- 1e-3 
spp <- -sp/(1 + sp)
mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                          drv = seq.int(nd))

mcf1s <-  oncoSimulIndiv(mcf1,
                         model = "McFL", 
                         mu = 1e-7,
                         detectionSize = 1e8, 
                         detectionDrivers = 100,
                         sampleEvery = 0.02,
                         keepEvery = 2,
                         initSize = 2000,
                         finalTime = 1000,
                         keepPhylog = TRUE,
                         onlyCancer = FALSE)

plot(mcf1s, type = "stacked", show = "genotypes")


## Fishplot: can it use drivers plot? Si. Pero más lioso. Muuuuuuucho más
## lioso. Olvidar.

mcf1s$other$PhylogDF
