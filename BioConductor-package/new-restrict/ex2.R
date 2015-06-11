load("examplesFitnessEffects.RData")

for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    ## if(i %in% c(20)) next
    tmp <-  oncoSimulIndiv(fE = examplesFitnessEffects[[i]],
                           model = "Bozic", 
                           mu = 1e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = 4,
                           sampleEvery = 2,
                           max.num.tries = 100,
                           initSize = 2000,
                           onlyCancer = TRUE)
}

for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    ## if(i %in% c(11, 20)) next
    tmp <-  oncoSimulIndiv(fE = examplesFitnessEffects[[i]],
                           model = "Exp", 
                           mu = 1e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = 4,
                           sampleEvery = 2,
                           max.num.tries = 100,
                           initSize = 2000,
                           onlyCancer = TRUE)
}


for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    ## if(i %in% c(11, 20)) next
    tmp <-  oncoSimulIndiv(fE = examplesFitnessEffects[[i]],
                           model = "McFL", 
                           mu = 5e-7,
                           detectionSize = 1e8, 
                           detectionDrivers = 2,
                           sampleEvery = 0.025,
                           max.num.tries = 10,
                           initSize = 2000,
                           onlyCancer = TRUE)
}




