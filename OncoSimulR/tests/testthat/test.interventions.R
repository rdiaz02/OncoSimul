inittime <- Sys.time()
cat(paste("\n Starting interventions tests", date(), "\n"))

test_that("1.1 Drastically reducing A population (McFL)", {

    fa1 <- data.frame(Genotype = c("WT", "A", "B"),
                    Fitness = c("1+n_*0",
                                "1",
                                "1"))

    afd3 <- allFitnessEffects(genotFitness = fa1,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

    ep1 <- oncoSimulIndiv(
                    afd3, 
                    model = "McFL",
                    mu = 1e-4,
                    #initSize = c(20000, 20000),
                    #initMutant = c("A", "B"),
                    initSize = 2000,
                    finalTime = 5.2,
                    sampleEvery = 0.01,
                    onlyCancer = FALSE
                    )

    # we plot the simulation when no interventions are specified.
    plot(ep1, show = "genotypes", type = "line")

    # now we especify intervention to drastically reduce A population
    interventions <- list(
    list(ID           = "intOverA",
        Trigger       = "(T >= 5)",
        WhatHappens   = "n_B = n_B * 0.1",
        Repetitions   = 0,
        Periodicity   = Inf
    ))

    interventions <- createInterventions(interventions, afd3)

    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "McFL",
                    mu = 1e-4,
                    #initSize = c(20000, 20000),
                    #initMutant = c("A", "B"),
                    initSize = 2000,
                    sampleEvery = 0.01,
                    finalTime = 5.2,
                    onlyCancer = FALSE,
		            interventions = interventions
                    )

    # we plot the simulation when interventions are specified.
    plot(ep2, show = "genotypes", type = "line")
})

test_that("1.2 Drastically reducing A population (Exp)", {

    fa1 <- data.frame(Genotype = c("WT", "A", "B"),
                    Fitness = c("n_*0",
                                "1",
                                "1"))

    afd3 <- allFitnessEffects(genotFitness = fa1,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

    ep1 <- oncoSimulIndiv(
                    afd3, 
                    model = "Exp",
                    mu = 1e-4,
                    sampleEvery = 0.01,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    finalTime = 5.2,
                    onlyCancer = FALSE
                    )

    # we plot the simulation when no interventions are specified.
    plot(ep1, show = "genotypes", type = "line")

    # now we especify intervention to drastically reduce A population
    interventions <- list(
        list(ID           = "intOverA",
            Trigger       = "(T >= 5)",
            WhatHappens   = "n_A = n_A * 0.1",
            Repetitions   = 0,
            Periodicity   = Inf
        )
    )

    interventions <- createInterventions(interventions, afd3)

    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "Exp",
                    mu = 1e-4,
                    sampleEvery = 0.01,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    finalTime = 5.2,
                    onlyCancer = FALSE,
		            interventions = interventions
                    )

    # we plot the simulation when interventions are specified.
    plot(ep2, show = "genotypes", type = "line")
})

test_that("1.3 Drastically reducing a high-fitness genotype population (McFL)", {
    df3x <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                       Fitness = c("0*n_",
                                    "1.5",
                                    "1.001",
                                    "1.002",
                                    "1.003",
                                    "1.004"))

    afd3 <- allFitnessEffects(genotFitness = df3x,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

    N <- 6
    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "McFL",
                    mu = 1e-4,
                    sampleEvery = 0.01,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    finalTime = 5.2,
                    onlyCancer = FALSE
                    )

    # we plot the simulation when no interventions are specified.
    plot(ep2, show = "genotypes", type = "line")

    interventions <- list(
    list(ID           = "intOverC",
        Trigger       = "(T >= 4)",
        WhatHappens   = "n_C = n_C * 0.1",
        Repetitions   = 0,
        Periodicity   = Inf
    ))

    interventions <- createInterventions(interventions, afd3)

    N <- 6
    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "McFL",
                    mu = 1e-4,
                    sampleEvery = 0.01,
                    initSize = 40000,
                    keepEvery = 0.1,
                    finalTime = 5.2,
                    onlyCancer = FALSE,
                    interventions = interventions
                    )
    # we plot the simulation when interventions are specified.
    plot(ep2, show = "genotypes", type = "line")
})

test_that("1.4 Drastically reducing a high-fitness genotype population (Exp)", {
    df3x <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                       Fitness = c("0*n_",
                                    "1.5",
                                    "1.001",
                                    "1.002",
                                    "1.003",
                                    "1.004"))

    afd3 <- allFitnessEffects(genotFitness = df3x,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

    N <- 6
    ep1 <- oncoSimulIndiv(
                    afd3, 
                    model = "Exp",
                    mu = 1e-4,
                    initMutant = c(),
                    initSize = c(),
                    sampleEvery = 0.01,
                    finalTime = 5.2,
                    )

    # we plot the simulation when no interventions are specified.
    plot(ep1, show = "genotypes", type = "line")

    interventions <- list(
    list(ID           = "intOverB",
        Trigger       = "(n_B >= n_)",
        WhatHappens   = "n_B = n_B * 0.001",
        Repetitions   = 0,
        Periodicity   = 0.001
    ))

    interventions <- createInterventions(interventions, afd3)

    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "Exp",
                    mu = 1e-4,
                    initSize = 40000,
                    sampleEvery = 0.001,
                    finalTime = 5.2,
                    interventions = interventions
                    )
    plot(ep2, show = "genotypes", type = "line")
})

