load("examplesFitnessEffects.RData")
## seeds 1 and 9 do it
for(i in 1:20) {
    cat(paste(" i = ", i, "\n"))
    set.seed(i)
    tmp <-  oncoSimulIndiv(fE = examplesFitnessEffects$cbn1,
                       model = "Bozic", 
                       mu = 1e-6,
                       detectionSize = 1e8, 
                       detectionDrivers = 4,
                       sampleEvery = 2,
                       initSize = 2000,
                       onlyCancer = TRUE)
}
