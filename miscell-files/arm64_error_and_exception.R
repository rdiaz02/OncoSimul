library(OncoSimulR)
packageVersion("OncoSimulR")
## Rethink how we deal with some exceptions. No object
## should be returned, probably.

## From fdfmutex2 but force it to fail on linux too
r1fd <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                   Fitness = c("1",
                               "1.4 + 1*(f_2)",
                               "1.4 + 1*(f_1)",
                               "1.6 + f_1 + f_2"))
afe4 <- allFitnessEffects(genotFitness = r1fd, 
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

mtfd <- allMutatorEffects(epistasis = c("A" = 0.1,
                                        "B" = 10))

## for(i in 1:1000) {
##     cat("\n ####################################### \n")
##     set.seed(i)
##     ma <- runif(1, 1e-5, 1e5)
##     mb <- runif(1, 1e-5, 1e5)
##     cat("\n i = ", i)
##     cat("\n ma = ", ma)
##     cat("\n mb = ", mb)
##     mtfd <- allMutatorEffects(epistasis = c("A" = ma,
##                                             "B" = mb))
##     s2fd <- oncoSimulIndiv(afe4,
##                            muEF = mtfd,
##                            model = "McFL", 
##                            onlyCancer = FALSE, 
##                            finalTime = 40,
##                            mu = 1e-4,
##                            initSize = 5000, 
##                            keepPhylog = TRUE,
##                            seed = NULL, 
##                            errorHitMaxTries = FALSE, 
##                            errorHitWallTime = FALSE)
##     s2fd
## }

i <- round(runif(1, 1, 1e5))
j <- 0

for(i in 1:20000) {
    cat("\n ####################################### \n")
    j <- j + 1
    i <- i + j
    
    ## ma <- runif(1, 1e-3, 1)
    ## mb <- runif(1, ma, 1000 * ma)
    ma <- 0.1
    mb <- 10
    cat("\n j, i.e., number of iterations = ", j)
    cat("\n i = ", i)
    cat("\n ma = ", ma)
    cat("\n mb = ", mb)
    cat("\n")
    mtfd <- allMutatorEffects(epistasis = c("A" = ma,
                                            "B" = mb))
    s2fd <- oncoSimulIndiv(afe4,
                           muEF = mtfd,
                           model = "McFL", 
                           onlyCancer = FALSE, 
                           finalTime = 40,
                           mu = 1e-4,
                           initSize = 5000, 
                           keepPhylog = TRUE,
                           seed = NULL, 
                           errorHitMaxTries = FALSE, 
                           errorHitWallTime = FALSE)
    s2fd
}





## plot(s2fd, show = "genotypes")
