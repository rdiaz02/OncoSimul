## None should crash or give an uncaught error
## Just a minimal set. Will later check warnings when they should, etc.

## Like the one in tests, but using all of them. Note also the old
## indexing of [16], because we are using all. Could have [i], but this is
## telling me we are using the original set.

data(examplesFitnessEffects)
## RNGkind("Mersenne-Twister")

## sometimes cancer is not reached. No problem.

## Very rarely, popSize > 1e15, and we get an exception. Decrease
## sampleEvery. And e2 only has two genes.

for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    if (names(examplesFitnessEffects)[16] == "e2") {
        detectionDrv <- 2
        sE <- 0.01
    } else {
        detectionDrv <- 4
        sE <- .05
    }
    tmp <-  oncoSimulIndiv(examplesFitnessEffects[[i]],
                           model = "Bozic", 
                           mu = 1e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = detectionDrv,
                           sampleEvery = sE,
                           max.num.tries = 100,
                           initSize = 2000,
                           onlyCancer = FALSE, detectionProb = NA)
    expect_true(inherits(tmp, "oncosimul2"))
}

for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    if (names(examplesFitnessEffects)[16] == "e2") {
        detectionDrv <- 2
        sE <- 0.01
    } else {
        detectionDrv <- 4
        sE <- .05
    }
    tmp <-  oncoSimulIndiv(examplesFitnessEffects[[i]],
                           model = "Exp", 
                           mu = 1e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = detectionDrv,
                           sampleEvery = sE,
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
    if (names(examplesFitnessEffects)[16] == "e2") {
        sE <- 0.01
    } else {
        sE <- .05
    }
    tmp <-  oncoSimulSample(4, examplesFitnessEffects[[i]],
                            onlyCancer = FALSE, detectionProb = NA,
                            sampleEvery = sE)
    expect_true(inherits(tmp, "list"))
}



for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    if (names(examplesFitnessEffects)[16] == "e2") {
        detectionDrv <- 2
        sE <- 0.01
    } else {
        detectionDrv <- 4
        sE <- .05
    }
    tmp <-  oncoSimulPop(4, examplesFitnessEffects[[i]],
                         onlyCancer = FALSE, detectionProb = NA,
                         detectionDrivers = detectionDrv,
                         sampleEvery = sE,
                         mc.cores = 2)
    expect_true(inherits(tmp, "oncosimulpop"))
    tmp2 <- samplePop(tmp)
    expect_true(inherits(tmp2, "matrix"))
}

