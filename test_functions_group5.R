library(OncoSimulR)
library(ggplot2)
library(testthat)

## Summary: In this script we are going to test that the functions "compositionPop2" receives an
## oncosimulpop object and works properly giving a data frame to "simul_boxplot2"
## Code for these functions is available at "Funciones_plot_markdown.R" 

############################ Define variables #######################

## Fitness equations
avc <- function (a, v, c) {
  data.frame(Genotype = c("WT", "G", "V", "A"),
             Fitness = c("1",
                         paste0("1 + ", a, " * (f_1 + 1)"),
                         paste0("1 + ", a, " * f_1 + ", v, " * (f_2 + 1) - ", c),
                         paste0("1 + ", a, " * f_1 + ", v, " * f_2")))
}

afavc <- allFitnessEffects(genotFitness = avc(37.5, 2, 1),
                           frequencyDependentFitness = TRUE,
                           frequencyType = "rel")

## Next step is to run several simulation, plot the trajectorie and box plot.
## We are going to check that the functions "compositionPop2" receives an
## oncosimulpop object and works properly giving a data frame to "simul_boxplot2"


################################# Tests #################################

## Function to check if we obtain a oncosimulpop object
is.oncosimulpop <- function(x) inherits(x, "oncosimulpop")

test_that("simuloncoPop gives an oncosimulpop object", {
  simulation <- oncoSimulPop(5,
                             mc.cores = 6,
                             afavc,
                             model = "McFL",
                             onlyCancer = FALSE,
                             finalTime = 25,
                             mu = 1e-2,
                             initSize = 4000,
                             keepPhylog = TRUE,
                             seed = NULL,
                             errorHitMaxTries = FALSE,
                             errorHitWallTime = FALSE)
  
  expect_true(is.oncosimulpop(simulation))
})

##################### Testing that "compositionPop2" works properly

test_that("sapply should give me a matrix (more friendly) to plot it", {
  simulation <- oncoSimulPop(3,
                             mc.cores = 6,
                             afavc,
                             model = "McFL",
                             onlyCancer = FALSE,
                             finalTime = 25,
                             mu = 1e-2, ## Notice the high value for mu
                             initSize = 4000,
                             keepPhylog = TRUE,
                             seed = NULL,
                             errorHitMaxTries = FALSE,
                             errorHitWallTime = FALSE)
  
  ## sapply should give me a matrix (more friendly)
  listPop <- sapply(simulation, function(x) tail(x[[1]], 1)[1, -1])
  expect_true(is.matrix(listPop))
})

## Check we get one value for each genotype in each simulation
test_that("Matrix should have one value for each genotype (rectangular)", {
  simulation <- oncoSimulPop(3,
                             mc.cores = 6,
                             afavc,
                             model = "McFL",
                             onlyCancer = FALSE,
                             finalTime = 25,
                             mu = 1e-2, ## Notice the high value for mu
                             initSize = 4000,
                             keepPhylog = TRUE,
                             seed = NULL,
                             errorHitMaxTries = FALSE,
                             errorHitWallTime = FALSE)
  
  ## sapply should give me a matrix (more friendly)
  listPop <- sapply(simulation, function(x) tail(x[[1]], 1)[1, -1])
  
  num_genotypes <- c("WT", simulation[[1]]$geneNames)
  
  # Each column should contains results for all possible genotype (4 rows)
  expect_true(dim(listPop)[1] == length(num_genotypes))
})

## Unfortunately this fails...
## Check that everything is correct when we change mu value
test_that("sapply should give me a matrix (more friendly) to plot it", {
  simulation <- oncoSimulPop(10,
                             mc.cores = 6,
                             afavc,
                             model = "McFL",
                             onlyCancer = FALSE,
                             finalTime = 25,
                             mu = 1e-6,  ## Notice the low value for mu
                             initSize = 4000,
                             keepPhylog = TRUE,
                             seed = NULL,
                             errorHitMaxTries = FALSE,
                             errorHitWallTime = FALSE)
  
  ## sapply should give me a matrix (more friendly)
  listPop <- sapply(simulation, function(x) tail(x[[1]], 1)[1, -1])
  expect_true(is.matrix(listPop))
})


## Something is wrong. We are recovering the final result from each simulation
## sapply should receive 4 values (from 4 genotypes) and build a matrix
## What are we receiving?
## A list...?
test_that("sapply gives a list with low mu values", {
  simulation <- oncoSimulPop(10,
                             mc.cores = 6,
                             afavc,
                             model = "McFL",
                             onlyCancer = FALSE,
                             finalTime = 25,
                             mu = 1e-6, ## Notice the low value for mu
                             initSize = 4000,
                             keepPhylog = TRUE,
                             seed = NULL,
                             errorHitMaxTries = FALSE,
                             errorHitWallTime = FALSE)
  
  ## sapply should give me a matrix (more friendly)
  listPop <- sapply(simulation, function(x) tail(x[[1]], 1)[1, -1])
  expect_true(is.list(listPop))
})

## If sapply does not give us a matrix (rectangular) is because there 
## are not the same number of genotypes in each simulation 
test_that("There is the same number of genotypes results for each
          simulation", {
  simulation <- oncoSimulPop(10,
                             mc.cores = 6,
                             afavc,
                             model = "McFL",
                             onlyCancer = FALSE,
                             finalTime = 25,
                             mu = 1e-6, ## Notice the low value for mu
                             initSize = 4000,
                             keepPhylog = TRUE,
                             seed = NULL,
                             errorHitMaxTries = FALSE,
                             errorHitWallTime = FALSE)
  
  listPop <- sapply(simulation, function(x) tail(x[[1]], 1)[1, -1])
  
  lengths_simulation <- sapply(listPop, function(x){length(x)})
  num_genotypes <- c("WT", simulation[[1]]$geneNames)
  
  expect_true(all(lengths_simulation == length(num_genotypes)))
})


## ATTENTION!!!!!
## Be careful. When the mutation rate is too low, the WT disappears 
## too quickly and the simulation does not last long enough, there 
## may not be time for all the genotypes to appear.  Raising the 
## mutation rate may be a solution, or increasing the birth rate of the WT.


###############################################################################











