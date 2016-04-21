## WIth low dummy mutations, the W and R and, then, eventually, the pE can
## be screwed up. See notes in the C++ code. So this is done and settled.


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
                               max.num.tries = 10,
                               keepEvery = 100,
                               initSize = 2000,
                               finalTime = 3000,
                               onlyCancer = TRUE,
                       keepPhylog = TRUE,
                       verbosity = 1)
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
                       minDetectDrvCloneSz = 2,
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
                       minDetectDrvCloneSz = 100,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp



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
                       minDetectDrvCloneSz = 500,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp


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
                       minDetectDrvCloneSz = 1000,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp


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
                       minDetectDrvCloneSz = 2000,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp

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
                       minDetectDrvCloneSz = 3000,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp


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
                       minDetectDrvCloneSz = 3500,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp




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
                       minDetectDrvCloneSz = 3200,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp

## nietieher -13 nor -11
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
                       minDetectDrvCloneSz = 3480,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp


## not -13
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
                       minDetectDrvCloneSz = 3470,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp


## not -13
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
                       minDetectDrvCloneSz = 3460,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)         
tmp



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
                       minDetectDrvCloneSz = 3477,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)
tmp


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
                       minDetectDrvCloneSz = 3478,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)
tmp


## with -11 we suddenly cannot get here and stop at 3218
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
                       minDetectDrvCloneSz = 3479,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)
tmp


## no number limit


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
                       finalTime = 1076.03,
                       onlyCancer = FALSE,
                       keepPhylog = TRUE,
                       minDetectDrvCloneSz = 3479,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)
tmp






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
                       finalTime = 1078,
                       onlyCancer = FALSE,
                       keepPhylog = TRUE,
                       minDetectDrvCloneSz = 3479,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)
tmp

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
                       finalTime = 1077,
                       onlyCancer = FALSE,
                       keepPhylog = TRUE,
                       minDetectDrvCloneSz = 3479,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)
tmp


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
                       finalTime = 1076,
                       onlyCancer = FALSE,
                       keepPhylog = TRUE,
                       minDetectDrvCloneSz = 3479,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)
tmp


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
                       finalTime = 1075,
                       onlyCancer = FALSE,
                       keepPhylog = TRUE,
                       minDetectDrvCloneSz = 3479,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)
tmp
















## with -11 we suddenly cannot get here and stop at 3218
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
                       minDetectDrvCloneSz = 3700,
                       errorHitMaxTries = FALSE,
                       verbosity = 1,
                       max.wall.time = 1000)
tmp





tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                               model = "McFL", 
                               mu = 5e-5,
                               detectionSize = 1e8, 
                               detectionDrivers = 3,
                               sampleEvery = 0.03,
                               max.num.tries = 10,
                               keepEvery = 100,
                               initSize = 2000,
                               finalTime = 3000,
                               onlyCancer = TRUE,
                               keepPhylog = TRUE)
