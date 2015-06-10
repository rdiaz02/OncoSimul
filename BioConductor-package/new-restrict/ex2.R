cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = 0.1,
                 sh = -0.9,
                 typeDep = "MN")
fcs <- allFitnessEffects(cs)
gfs <- evalAllGenotypes(fcs, order = FALSE)
gfs


initSize <- 500
sampleEvery <- 1
## remember I am also using R's rangen
set.seed(3) 
tmp <- nr_oncoSimul.internal(fcs,
                             birth = -99,
                             death = 1,
                             mu = 1e-5,
                             initSize = initSize,
                             sampleEvery = sampleEvery,
                             detectionSize = 1e6,
                             finalTime = 365 * 25,
                             initSize_species = 2000,
                             initSize_iter = 500,
                             seed = 9,
                             verbosity = 1,  ## 2 for verbose
                             speciesFS = 40000,
                             ratioForce = 2,
                             typeFitness = "bozic1",
                             max.memory = 2000,
                             mutatorGenotype = 0,
                             initMutant = NULL,
                             max.wall.time = 200,
                             keepEvery = 1,
                             alpha = 0.0015,
                             K = 100/(exp(1) - 1),
                             detectionDrivers = 4,
                             onlyCancer = TRUE,
                             errorHitWallTime = TRUE,
                             max.num.tries = 500,
                             errorHitMaxTries = TRUE,
                             minDDrPopSize = 0,
                             extraTime = 0)


## tmp <- nr_oncoSimul.internal(p701nr,
##                              birth = -99,
##                              death = 1,
##                              mu = 1e-6,
##                              initSize = initSize,
##                              sampleEvery = sampleEvery,
##                              detectionSize = 1e8,
##                              finalTime = 0.25 * 365 * 25,
##                              initSize_species = 2000,
##                              initSize_iter = 500,
##                              seed = NULL,
##                              verbosity = 1,  ## 2 for verbose
##                              speciesFS = 40000,
##                              ratioForce = 2,
##                              typeFitness = "bozic1",
##                              max.memory = 2000,
##                              mutatorGenotype = 0,
##                              initMutant = NULL,
##                              max.wall.time = 200,
##                              keepEvery = 1,
##                              alpha = 0.0015,
##                              K = 100/(exp(1) - 1),
##                              detectionDrivers = 4,
##                              onlyCancer = FALSE,
##                              errorHitWallTime = TRUE,
##                              max.num.tries = 500,
##                              errorHitMaxTries = TRUE,
##                              minDDrPopSize = 0,
##                              extraTime = 0)
## dim(tmp[[1]])



p701nr <- allFitnessEffects(data.frame(parent = c("Root", rep("1", 4), 2, 3, 4, 4, 5),
                                       child = c(1, 2, 3, 4, 5, 6, 6, 6, 7, 7),
                                       s = 0.1,
                                       sh = -Inf,
                                       typeDep = "MN"))
evalAllGenotypes(p701nr, order = FALSE)

tmp <-  oncoSimulIndiv(fitnessEffects = p701nr,
                       model = "Bozic", numPassengers = 0, mu = 1e-6,
                       detectionSize = 1e8, detectionDrivers = 4,
                       sampleEvery = sampleEvery,
                       initSize = initSize, s = 5, sh = -1,
                       K = initSize/(exp(1) - 1), keepEvery = sampleEvery,
                       minDDrPopSize = "auto",
                       extraTime = 0,
                       finalTime = 0.25 * 25 * 365, onlyCancer = TRUE,
                       max.memory = 2000, max.wall.time = 200,
                       max.num.tries = 500,
                       errorHitWallTime = TRUE,
                       errorHitMaxTries = TRUE,
                       verbosity = 0,
                       seed = NULL)
dim(tmp[[1]]); tmp[[1]][nrow(tmp[[1]]), ]




mc1 <- oncoSimulIndiv(fitnessEffects = p701nr,
                      model = "McFL",
                      mu = 5e-7,
                      initSize = 4000,
                      sampleEvery = 0.025,
                      finalTime = 15000,
                      keepEvery = 5)



data(examplePosets)
p701 <- examplePosets[["p701"]]


tmpo <-  oncoSimulIndiv(p701, model = "Bozic", numPassengers = 0, mu = 1e-6,
                        detectionSize = 1e8, detectionDrivers = 4,
                        sampleEvery = sampleEvery,
                        initSize = initSize, s = 0.1, sh = -1,
                        K = initSize/(exp(1) - 1), keepEvery = sampleEvery,
                        minDDrPopSize = "auto",
                        extraTime = 0,
                        finalTime = 0.25 * 25 * 365, onlyCancer = FALSE,
                        max.memory = 2000, max.wall.time = 200,
                        max.num.tries = 500,
                        errorHitWallTime = TRUE,
                        errorHitMaxTries = TRUE,
                        verbosity = 2)

