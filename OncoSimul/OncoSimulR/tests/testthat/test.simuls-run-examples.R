inittime <- Sys.time()
cat(paste("\n Starting simuls-runs-examples tests", date(), "\n"))

## None should crash or give an uncaught error
## Just a minimal set. Will later check warnings when they should, etc.

data(examplesFitnessEffects)
## RNGkind("Mersenne-Twister")
## sometimes cancer is not reached. No problem.

## Very rarely, popSize > 1e15, and we get an exception. Decrease
## sampleEvery. And e2 only has two genes.

## We take a sample. All of the 22 are run in the long tests.
nex <- 5
examplesFitnessEffects <- examplesFitnessEffects[
    sample(length(examplesFitnessEffects), nex)]

for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    if (names(examplesFitnessEffects)[i] == "e2") {
        detectionDrv <- 2
        sE <- 0.05
    } else {
        detectionDrv <- 4
        sE <- 0.02 ## smaller than usual but very rarely (< 1/1000) I can
                  ## get crashes of Algo 2 as popSize > 1e15
    }
    tmp <-  oncoSimulIndiv(examplesFitnessEffects[[i]],
                           model = "Bozic", 
                           mu = 1e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = detectionDrv,
                           sampleEvery = sE, keepEvery = 1,
                           max.num.tries = 100,
                           initSize = 2000,
                           onlyCancer = FALSE, detectionProb = NA)
    expect_true(inherits(tmp, "oncosimul2"))
}

for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    if (names(examplesFitnessEffects)[i] == "e2") {
        detectionDrv <- 2
        sE <- 0.05
    } else {
        detectionDrv <- 4
        sE <- 0.02
    }
    tmp <-  oncoSimulIndiv(examplesFitnessEffects[[i]],
                           model = "Exp", 
                           mu = 1e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = detectionDrv,
                           sampleEvery = sE, keepEvery = 1,
                           max.num.tries = 100,
                           initSize = 2000,
                           onlyCancer = FALSE, detectionProb = NA)
    expect_true(inherits(tmp, "oncosimul2"))
}


for(i in 1:length(examplesFitnessEffects)) {
    tmp <-  oncoSimulIndiv(examplesFitnessEffects[[i]],
                           model = "McFL", 
                           mu = 5e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = 2,
                           sampleEvery = 0.025,
                           keepEvery = 1,
                           max.num.tries = 10,
                           initSize = 2000,
                           finalTime = 15000,
                           onlyCancer = FALSE, detectionProb = NA)
    expect_true(inherits(tmp, "oncosimul2"))
}



for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
        cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    if (names(examplesFitnessEffects)[i] == "e2") {
        sE <- 0.05
    } else {
        sE <- 0.1
    }
    tmp <-  oncoSimulSample(4, examplesFitnessEffects[[i]],
                            onlyCancer = FALSE, detectionProb = NA,
                            sampleEvery = sE)
    expect_true(inherits(tmp, "list"))
}



for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    if (names(examplesFitnessEffects)[i] == "e2") {
        detectionDrv <- 2
        sE <- 0.05
    } else {
        detectionDrv <- 4
        sE <- 0.02
    }
    tmp <-  oncoSimulPop(4, examplesFitnessEffects[[i]],
                         onlyCancer = FALSE, detectionProb = NA,
                         detectionDrivers = detectionDrv,
                         sampleEvery = sE, keepEvery = 1,
                         mc.cores = 2)
    expect_true(inherits(tmp, "oncosimulpop"))
    tmp2 <- samplePop(tmp)
    expect_true(inherits(tmp2, "matrix"))
}

cat(paste("\n Ending simuls-runs-examples tests", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
