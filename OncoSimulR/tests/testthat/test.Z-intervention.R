inittime <- Sys.time()
cat(paste("\n Starting interventions tests", date(), "\n"))

test_that("1. Drastically reducing a high-fitness genotype population (McFL) | Trigger depends on T and n_*", {
    set.seed(1)
    df3x <- data.frame(Genotype = c("WT", "B", "R"),
                       Fitness = c("1 + (n_ * 0)",
                                    "1.5",
                                    "1.003 + 0.002 * (n_B > 120)"))

    afd3 <- allFitnessEffects(genotFitness = df3x,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

    interventions <- list(
        list(
            ID          = "intOverBAffectsR",
            Trigger     = "((T >= 20) and (T <= 70)) and (n_B >= 200)",
            WhatHappens = "n_B = 100",
            Repetitions = 50,
            Periodicity = 0.001
        ),
        list(
            ID          = "intOverTotPop",
            Trigger     = "(T >= 80) and (T <= 85) and (n_B >= 40)",
            WhatHappens = "n_B = 20",
            Repetitions = Inf,
            Periodicity = 0.001
        )
    )

    interventions <- createInterventions(interventions, afd3)

    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "McFLD",
                    mu = 1e-4,
                    sampleEvery = 0.001,
                    initSize = c(5000, 10, 300),
                    initMutant = c("WT", "B", "R"),
                    finalTime = 100,
                    onlyCancer = FALSE,
                    interventions = interventions
                    )

    flag <- FALSE
    i <- 20002
    while(i <= 70001){
        if(ep2$pops.by.time[i, 3:3] >= 210){
            flag <- TRUE
        }
        i <- i + 1
    }

    testthat::expect_equal(flag, FALSE)

    # then, between the time intervals, T >= 80 and T<=85
    # we control that the B population
    flag <- FALSE
    i <- 80002
    while(i <= 85000){
        if(ep2$pops.by.time[i, 3:3] > 40){
            flag <- TRUE
        }
        i <- i + 1
    }

    testthat::expect_equal(flag, FALSE)

    # we plot the simulation when no interventions are specified.
    #plot(ep2, show = "genotypes", type = "line")
})





test_that("2. Drastically reducing a high-fitness genotype population (Exp) | Trigger depends on n_*", {
    set.seed(1)
    df3x <- data.frame(Genotype = c("WT", "B", "R"),
                       Fitness = c("1 + (n_ * 0)",
                                    "1.5",
                                    "1.003 + 0.002 * (n_B > 120)"))

    afd3 <- allFitnessEffects(genotFitness = df3x,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

    interventions <- list(
        list(
            ID          = "intOverBAffectsR",
            Trigger     = "((T >= 20) and (T <= 70)) and (n_B >= 200)",
            WhatHappens = "n_B = 100",
            Repetitions = 50,
            Periodicity = 0.001
        ),
        list(
            ID          = "intOverTotPop",
            Trigger     = "(T >= 80) and (T <= 85) and (n_B >= 40)",
            WhatHappens = "n_B = 20",
            Repetitions = Inf,
            Periodicity = 0.001
        )
    )

    interventions <- createInterventions(interventions, afd3)

    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "Exp",
                    mu = 1e-4,
                    sampleEvery = 0.001,
                    initSize = c(5000, 10, 300),
                    initMutant = c("WT", "B", "R"),
                    finalTime = 100,
                    onlyCancer = FALSE,
                    interventions = interventions)

    flag <- FALSE
    i <- 20002
    while(i <= 70001){
        if(ep2$pops.by.time[i, 3:3] >= 210){
            flag <- TRUE
        }
        i <- i + 1
    }

    testthat::expect_equal(flag, FALSE)

    # then, between the time intervals, T >= 80 and T<=85
    # we control that the B population
    flag <- FALSE
    i <- 80002
    while(i <= 85000){
        if(ep2$pops.by.time[i, 3:3] > 40){
            flag <- TRUE
        }
        i <- i + 1
    }

    testthat::expect_equal(flag, FALSE)
    # we plot the simulation when no interventions are specified.
    #plot(ep2, show = "genotypes", type = "line")
})

cat(paste("\n Ending interventions tests", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

