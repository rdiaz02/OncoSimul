inittime <- Sys.time()
cat(paste("\n Starting plotMuller at", date()))

test_that("oncosimul v.2 objects and genotype plotting", {
   data(examplesFitnessEffects)
      p1 <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL",
                       mu = 5e-5,
                       detectionSize = 1e8,
                       detectionDrivers = 3,
                       sampleEvery = 0.025,
                       max.num.tries = 10,
                       keepEvery = 5,
                       initSize = 2000,
                       finalTime = 3000,
                       onlyCancer = FALSE,
                       keepPhylog = TRUE)
 class(p1)
 plot(p1, type = "muller")
})

test_that("only recognized arguments", {
  data(examplesFitnessEffects)
  simulWithoutPhyLog<-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                         model = "McFL",
                         mu = 5e-5,
                         detectionSize = 1e8,
                         detectionDrivers = 3,
                         sampleEvery = 0.025,
                         max.num.tries = 10,
                         keepEvery = 5,
                         initSize = 2000,
                         finalTime = 3000,
                         onlyCancer = FALSE,
                         keepPhylog = FALSE)

  expect_error(plot(simulWithoutPhyLog, type = "muller"),
               "Object simulation must has property: other$PhylogDF", fixed = TRUE)
})

test_that("OncoSimul class", {
  data(examplePosets)
  p705 <- examplePosets[["p705"]]
  simulClassOncosimul1 <- oncoSimulIndiv(p705, model = "McFL",
                         mu = 5e-6,
                         sampleEvery = 0.02,
                         keepEvery = 10,
                         initSize = 2000,
                         finalTime = 3000,
                         max.num.tries = 100,
                         onlyCancer = FALSE)

  expect_error(plot(simulClassOncosimul1, type = "muller"),
               "Type of object class must be: oncosimul2", fixed = TRUE)
})

test_that("only recognized arguments muller type", {
  data(examplesFitnessEffects)
  simulWithoutPhyLog<-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                                       model = "McFL",
                                       mu = 5e-5,
                                       detectionSize = 1e8,
                                       detectionDrivers = 3,
                                       sampleEvery = 0.025,
                                       max.num.tries = 10,
                                       keepEvery = 5,
                                       initSize = 2000,
                                       finalTime = 3000,
                                       onlyCancer = FALSE,
                                       keepPhylog = TRUE)
  
  expect_error(plot(simulWithoutPhyLog, type = "muller", muller.type="invent"),
               "Type of muller.plot unknown: it must be one offrequency or population", fixed = TRUE)
})


