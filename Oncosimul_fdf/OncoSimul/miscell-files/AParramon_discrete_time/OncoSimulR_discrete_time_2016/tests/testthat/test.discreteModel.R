cat(paste("\n Starting DiscreteModel tests", date(), "\n"))

##alternative definitions of  fitness effects:
orderFE <- allFitnessEffects(orderEffects = c("F > D" = -0.2, "D > F" = 0.2))
epistFE <- allFitnessEffects(epistasis = c("A:-B" = 0.1, "B:-A" = 0.6, "A : B" = 0.2))
noIntFE <- allFitnessEffects(noIntGenes = c("A" = -0.03, "B" = 0.5, "C" = 0.8, 
                                         "D" = -0.07, "E" = 0.21, "F" = 0.35,
                                         "G" = -0.06, "H" = -0.12, "I" = 0.19))
mixFE <- allFitnessEffects(orderEffects = c("F > D" = -0.2, "D > F" = 0.2),
                            epistasis = c("A:-B" = 0.1, "B:-A" = 0.4, "A : B" = 0.2),
                            noIntGenes = c("X" = -0.03, "Y" = 0.05, "Z" = 0.0))
##rest of variables
mu <- 1e-4
popIni <- 100000
tMax <- 1000
verbosity <- 0
keepEvery <- 0
sampleEvery <- 0


#test_that("error if any mu < 0", {
#    muAux <- runif(2,-0.001,0)
#    expect_error(simulDiscrete(rFE = orderFE, 
#                                    mu = muAux, 
#                                    popIni = popIni, 
#                                    tMax = tMax, 
#                                    seed = seed, 
#                                    itInfo = itInfo, 
#                                    verbosity = verbosity),
#                  "(at least one) mutation rate (mu) is negative", fixed = TRUE)
#})
#
#test_that("error if mu and genotype don't have the same length", {
#  muAux <- runif(3,0,0.001)
#  expect_error(simulDiscrete(rFE = orderFE, 
#                             mu = muAux, 
#                             popIni = popIni, 
#                             tMax = tMax, 
#                             seed = seed, 
#                             itInfo = itInfo, 
#                             verbosity = verbosity),
#               "When using per-gene mutation rates, ",
#               "there must be the same number of genes in the ",
#               "mu vector and the fitness effects object.", fixed = TRUE)
#})
#
#test_that("error if initial population < 1", {
#  expect_error(simulDiscrete(rFE = orderFE, 
#                             mu = mu, 
#                             popIni = 0, 
#                             tMax = tMax, 
#                             seed = seed, 
#                             itInfo = itInfo, 
#                             verbosity = verbosity),
#               "Max population < 1", fixed = TRUE)
#})
#
#
#test_that("error if max time < 1", {
#  expect_error(simulDiscrete(rFE = orderFE, 
#                             mu = mu, 
#                             popIni = popIni, 
#                             tMax = 0, 
#                             seed = seed, 
#                             itInfo = itInfo, 
#                             verbosity = verbosity),
#               "Max time < 1", fixed = TRUE)
#})

test_that("The number of cells is maintained around the maximum population",{
  popAcc <- 0
  popIniAux <- 100
  flag <- TRUE
  it <- 0

  while(flag && it<7){
    for(i in 1:400){
      s <- oncoSimulIndiv(orderFE,
                    model = 'DM',
                    mu = mu,
                    sampleEvery = sampleEvery,
                    keepEvery = keepEvery,
                    initSize = popIniAux, 
                    finalTime = 1, 
                    verbosity = 0)
      popAcc <- popAcc + s$TotalPopSize
    }
    popAcc = popAcc / 400
    if(popAcc >= 99.02 && popAcc <= 100.98){
      flag <- FALSE
    }
    it <- it+1
  }
  expect_true(!flag)
})

test_that("parameter verbosity is ok", {
  muAux <- 1e-2
  expect_output(oncoSimulIndiv(mixFE,
                            model = 'DM', 
                            mu = muAux,
                            sampleEvery = 20,
                            keepEvery = keepEvery,
                            initSize = 100, 
                            finalTime = 50, 
                            verbosity = 1),
               "Start", "End", fixed = TRUE)
  
  expect_output(oncoSimulIndiv(mixFE,
                              model = 'DM',
                              mu = muAux,
                              sampleEvery = 20,
                              keepEvery = keepEvery,
                              initSize = 100, 
                              finalTime = 50,  
                              verbosity = 2),
                "After mutation", fixed = TRUE)
})

test_that("parameter sampleEvery is ok", {
  expect_output(oncoSimulIndiv(orderFE,
                              model = 'DM', 
                              mu = mu, 
                              sampleEvery = 11,
                              keepEvery = keepEvery,
                              initSize = 100,
                              finalTime = 12,
                              verbosity = 2),
                "After mutation 11", fixed = TRUE)
})

#test_that("example with order effects works well", {
#  expect_true(simulDiscrete(rFE = orderFE, 
#                              mu = mu, 
#                              popIni = popIni, 
#                              tMax = tMax, 
#                              seed = seed, 
#                              itInfo = itInfo, 
#                              verbosity = verbosity)$LargestCloneGenotype == "D > F _ ")
#})
#
#test_that("example with epistastis works well", {
#  expect_true(simulDiscrete(rFE = epistFE, 
#                            mu = mu, 
#                            popIni = popIni, 
#                            tMax = tMax, 
#                            seed = seed, 
#                            itInfo = itInfo, 
#                            verbosity = verbosity)$LargestCloneGenotype == "B")
#})
#
#test_that("example with noInt works well", {
#  muAux <- runif(9,0,0.001)
#  expect_true(simulDiscrete(rFE = noIntFE, 
#                            mu = muAux, 
#                            popIni = popIni, 
#                            tMax = tMax, 
#                            seed = seed, 
#                            itInfo = itInfo, 
#                            verbosity = verbosity)$LargestCloneGenotype == "B, C, E, F, I")
#})
