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

seed0 <- round(runif(1, 1, 1e5))


for(iter in 1:20000) {
    cat("\n ####################################### \n")
    seed <- seed0 + iter
    set.seed(seed)
    ## ma <- runif(1, 1e-3, 1)
    ## mb <- runif(1, ma, 1000 * ma)
    ma <- 0.1
    mb <- 10
    cat("\n iter, i.e., number of iterations = ", iter)
    cat("\n seed = ", seed)
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
    print(s2fd)
    ## plot(s2fd)
}





## plot(s2fd, show = "genotypes")


## this gives an exception; from test.Z-oncoSimulIndivFDF.R

r <- data.frame(rfitness(2))

colnames(r)[which(colnames(r) == "Birth")] <- "Fitness"

r[, "Fitness"] <- c("10*f_", 
                    "10*f_1", 
                    "50*f_2", 
                    "200*(f_1 + f_2) + 50*f_1_2")


suppressWarnings(afe <- allFitnessEffects(genotFitness = r, 
                                          frequencyDependentFitness = TRUE, 
                                          frequencyType = "rel"))

set.seed(1)
osi <- oncoSimulIndiv(afe, 
                      model = "Bozic", 
                      onlyCancer = FALSE, 
                      finalTime = 5000, 
                      verbosity = 0, 
                      mu = 1e-6,
                      initSize = 500, 
                      keepPhylog = FALSE,
                      seed = NULL, 
                      errorHitMaxTries = TRUE, 
                      errorHitWallTime = TRUE)


## So we can add code to plot and summary. If uncoreverable exception,
## which is stored in $UnrecoverExcept, then
## - do not plot
## - give error message from summary

## I can also get crashes in line 12644 of vignette, plotintex2
