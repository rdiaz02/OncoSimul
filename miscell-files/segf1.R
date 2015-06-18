source("new-restrict.R")
load("examplesFitnessEffects.RData")
## seeds 1 and 9 do it
for(i in 1:100) {
    cat("\n #######################################################\n ")
    cat(paste(" This is the i = ", i, "\n"))
    
    set.seed(i)
    tmp <-  oncoSimulIndiv(fE = examplesFitnessEffects$cbn1,
                           model = "Bozic", 
                           mu = 1e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = 4,
                           sampleEvery = 2,
                           initSize = 2000,
                           seed =  2 * i,
                           onlyCancer = TRUE)
}
