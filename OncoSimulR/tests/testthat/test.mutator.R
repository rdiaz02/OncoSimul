fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                      "b : c" = 0.5),
                        noIntGenes = c("e" = 0.1))

evalAllGenotypes(fe, order = FALSE)

fm <- OncoSimulR:::allMutatorEffects(fe, noIntGenes = c("a" = 10,
                                                        "c" = 5))

OncoSimulR:::evalAllGenotypesMut(fm) ## OK

## should fail
fm <- OncoSimulR:::allMutatorEffects(fe, noIntGenes = c("a" = 10,
                                                        "d" = 5))



OncoSimulR:::evalGenotypeMut("a", fm)





evalAllGenotypes(
    allFitnessEffects(epistasis = c("a : b" = 0.3,
                                    "b : c" = 0.5)))


## Fitness of 0, but mutator effects.

## Modules same and different from fitness effects.



## oe <- c("C > F" = -0.1, "H > I" = 0.12)
sm <- c("I:J"  = -1)
sv <- c("-K:M" = -.5, "K:-M" = -.5)
epist <- c(sm, sv)

modules <- c("Root" = "Root",
             "I" = "i1", "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")

set.seed(1) ## for repeatability
noint <- rexp(5, 10)
names(noint) <- paste0("n", 1:5)

fea <- allFitnessEffects(epistasis = epist, 
                         noIntGenes = noint,
                         geneToModule = modules)
