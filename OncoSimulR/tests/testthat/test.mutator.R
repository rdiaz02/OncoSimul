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

OncoSimulR:::evalGenotypeMut("c", fm)





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



fea2 <- allFitnessEffects(noIntGenes = c("a" = .1, "b" = .2))

## Modules, but not interactions. Suppose I want to have a model where A
## and B are modules, and they have the usual multiplicative fitness.

## It is complicated, because noInts do not allow for modules. Well, not
## really. See below.


## This will not work
fe2 <- allFitnessEffects(epistasis = c("A" = 0.3,
                                      "B" = 0.5),
                        geneToModule = c("Root" = "Root",
                                         "A" = "a1, a2",
                                         "B" = "b1"))

## This does work, but A and B have wt
sa <- 0.1
sb <- 0.2


fe3 <- allFitnessEffects(epistasis = c("A : B" = sa * sb),
                        geneToModule = c("Root" = "Root",
                                         "A" = "a1, a2",
                                         "B" = "b1"))
evalAllGenotypes(fe3, order = FALSE)

## works, but not what we want
fe4 <- allFitnessEffects(epistasis = c(
                             "A : -B" = sa,
                             "-A : B" = sb,
                             "A : B" = sa * sb),
                        geneToModule = c("Root" = "Root",
                                         "A" = "a1, a2",
                                         "B" = "b1"))
evalAllGenotypes(fe4, order = FALSE)
                        

## works, but not what we want
fe5 <- allFitnessEffects(epistasis = c(
                             "A" = sa,
                             "B" = sb,
                             "A : B" = sa * sb),
                        geneToModule = c("Root" = "Root",
                                         "A" = "a1, a2",
                                         "B" = "b1"))
evalAllGenotypes(fe5, order = FALSE)


## works, and does what we want
fe6 <- allFitnessEffects(epistasis = c(
                             "A" = sa,
                             "B" = sb,
                             "A : B" = 0),
                        geneToModule = c("Root" = "Root",
                                         "A" = "a1, a2",
                                         "B" = "b1"))
evalAllGenotypes(fe6, order = FALSE)


## works, but not what we want
fe7 <- allFitnessEffects(epistasis = c(
                             "a1" = sa,
                             "a2" = sa,
                             "b1" = sb,
                             "a1 : a2" = -sa,
                             "a1 : b1" = sb,
                             "a2: b1" = sb))
evalAllGenotypes(fe7, order = FALSE)







## works, and does what we want
sc <- 0.05
fe8 <- allFitnessEffects(epistasis = c(
                             "A" = sa,
                             "B" = sb,
                             "C" = sc,
                             "A : B" = 0),
                        geneToModule = c("Root" = "Root",
                                         "A" = "a1, a2",
                                         "B" = "b1",
                                         "C" = "c1"))
evalAllGenotypes(fe8, order = FALSE)


## So explain why noIntGenes have no modules: because a module implicitly
## means an interaction: two or more things have the same effect.
