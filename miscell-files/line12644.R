library(OncoSimulR)
packageVersion("OncoSimulR")
seed0 <- round(runif(1, 1, 1e5))

for(iter in 1:20000) {
    cat("\n ####################################### \n")
    seed <- seed0 + iter
    set.seed(seed)
    cat("\n iter, i.e., number of iterations = ", iter)
    cat("\n seed = ", seed)

    ## eh??!!! with set.seed(1) this takes forever!
    
    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                        Fitness = c("1",
                                    "1 + 0.25 * (n_B > 0)",
                                    ".9 + 0.4 * (n_A > 0)"
                                    ))
    afd3 <- allFitnessEffects(genotFitness = gffd3,
                              frequencyDependentFitness = TRUE,
                              frequencyType = "abs")


    timerun <- system.time(osi <- oncoSimulIndiv( afd3,
                                                 model = "McFLD",
                                                 onlyCancer = FALSE,
                                                 finalTime = 200,
                                                 mu = 1e-4,
                                                 initSize = 5000,
                                                 sampleEvery = 0.001))[1]

    cat("\n Time it took was ", timerun, "\n")
    plot(osi, show = "genotypes", type = "line")
}
