## when dummy can be very low, say 1e-15, many things
## do not finish. Otherwise, they do, even if we get many of
## Note: mutation = 0; no positions left for mutation
## AND (and this is weird now to me) none of the
## "Note: updating in null mutation\n"

data(examplesFitnessEffects)
i <- 3
set.seed(i)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.03,
                       max.num.tries = 1,
                       keepEvery = 100,
                       initSize = 2000,
                       finalTime = 3000,
                       onlyCancer = TRUE,
                       keepPhylog = TRUE,
                       minDetectDrvCloneSz = 2,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp

## at 2, they have diverged

## with the mindetection of 259 they are already very different.



data(examplesFitnessEffects)
i <- 24232
set.seed(i)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.03,
                       max.num.tries = 1,
                       keepEvery = 100,
                       initSize = 2000,
                       finalTime = 3000,
                       onlyCancer = TRUE,
                       keepPhylog = TRUE,
                       minDetectDrvCloneSz = 159,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp






i <- 17
set.seed(i)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.03,
                       max.num.tries = 50,
                       keepEvery = 100,
                       initSize = 2000,
                       finalTime = 3000,
                       onlyCancer = TRUE,
                       keepPhylog = TRUE,
                       minDetectDrvCloneSz = 3000,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp













data(examplesFitnessEffects)
i <- 3
set.seed(i)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.03,
                       max.num.tries = 1,
                       keepEvery = 100,
                       initSize = 2000,
                       finalTime = 3000,
                       onlyCancer = TRUE,
                       keepPhylog = TRUE,
                       minDetectDrvCloneSz = 5000,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp
















data(examplesFitnessEffects)
i <- 3
set.seed(i)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.01,
                       max.num.tries = 1,
                       keepEvery = 100,
                       initSize = 2000,
                       finalTime = 3000,
                       onlyCancer = TRUE,
                       keepPhylog = TRUE,
                       minDetectDrvCloneSz = 2,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp





data(examplesFitnessEffects)
i <- 593424
set.seed(i)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = .1,
                       max.num.tries = 20,
                       keepEvery = 100,
                       initSize = 2000,
                       finalTime = 3000,
                       onlyCancer = TRUE,
                       keepPhylog = TRUE,
                       minDetectDrvCloneSz = 500,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp














## significant that it finishes later often?

##mindummy of 1e-15 is way too low probably and I think leads to strange
## numerical issues when updating pops.

## So try 1e-13 1e-17 1e-11 and compare
set.seed(135)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                               model = "McFL", 
                               mu = 5e-5,
                               detectionSize = 1e8, 
                               detectionDrivers = 3,
                               sampleEvery = 0.03,
                               max.num.tries = 100,
                               keepEvery = 100,
                               initSize = 2000,
                               finalTime = 3000,
                               onlyCancer = TRUE,
                       keepPhylog = TRUE,
                       verbosity = 1)
