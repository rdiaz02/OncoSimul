fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                      "b : c" = 0.5),
                        noIntGenes = c("e" = 0.1))
fm <- OncoSimulR:::allMutatorEffects(noIntGenes = c("a" = 10,
                                                    "c" = 5))

OncoSimulR:::evalGenotypeFitAndMut("a", fe, fm)
OncoSimulR:::evalGenotypeFitAndMut("b", fe, fm) 

OncoSimulR:::evalGenotypeFitAndMut("e", fe, fm)
OncoSimulR:::evalGenotypeFitAndMut("b, e", fe, fm)
OncoSimulR:::evalGenotypeFitAndMut("a, b, e", fe, fm)
OncoSimulR:::evalGenotypeFitAndMut("a, b, c, e", fe, fm)

fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                      "b : c" = 0.5),
                        noIntGenes = c("e" = 0.1))

fe <- allFitnessEffects(noIntGenes = c("a" = 0.2, "c" = 0.4, "d" = 0.6, "e" = 0.1))
fm <- allMutatorEffects(noIntGenes = c("a" = 10,
                                       "c" = 5))

oncoSimulIndiv(fe, muEF = fm)

oncoSimulIndiv(fe)


## test with var mut rate,
## run all tests
## create new tests

## docs:
##    - help
##  -fignete

## Fitness of 0, but mutator effects.

## Modules same and different from fitness effects.




fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                      "b : c" = 0.5),
                        noIntGenes = c("e" = 0.1))

evalAllGenotypes(fe, order = FALSE)

fm <- OncoSimulR:::allMutatorEffects(noIntGenes = c("a" = 10,
                                                    "c" = 5))

OncoSimulR:::evalAllGenotypesMut(fm) ## OK

## ## should fail
## fm <- OncoSimulR:::allMutatorEffects(noIntGenes = c("a" = 10,
##                                                     "d" = 5))

OncoSimulR:::evalGenotypeMut("a", fm)
OncoSimulR:::evalGenotypeMut("c", fm)
## should fail
OncoSimulR:::evalGenotypeMut("b", fm)














## Is fitness of wildtype always 0? Really? Evaluate it.
## It is: see evalGenotypeFitness
OncoSimulR:::evalRGenotype(vector(mode = "integer", length = 0), fe, TRUE, FALSE, "evalGenotype")



## here, we DO evaluate the length 0 genotype. Turn into test.
fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                      "b : c" = 0.5),
                        noIntGenes = c("e" = 0.1))
fm <- OncoSimulR:::allMutatorEffects(noIntGenes = c("a" = 10,
                                                    "c" = 5))
OncoSimulR:::evalRGenotypeAndMut(vector(mode = "integer", length = 0),
                                 fe,
                                 fm,
                                 OncoSimulR:::matchGeneIDs(fm, fe)$Reduced,
                                 TRUE, FALSE)

