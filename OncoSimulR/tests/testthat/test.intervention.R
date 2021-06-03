inittime <- Sys.time()
cat(paste("\n Starting interventions tests", date(), "\n"))
############################################################################################################################
############################################################################################################################
############################################################################################################################

test_that("A intervention is created correctly",{
    fa1 <- data.frame(Genotype = c("A", "B"),
                    Fitness = c("1+ 0*n_A",
                                "1.002"))

    afd3 <- allFitnessEffects(genotFitness = fa1,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

    interventions <- list(
        list(ID           = "intOverA",
            Trigger       = "n_B >= 5",
            WhatHappens   = "n_A = n_A * 0.1",
            Repetitions   = 0,
            Periodicity   = Inf
        )
    )

    testthat::expect_output(createInterventions(interventions, afd3), "Checking intervention: intOverA" )
    interventions <- createInterventions(interventions, afd3)

    # we check that the transformation of the WhatHappens and Trigger atributte is correct
    testthat::expect_equal(interventions[[1]]$WhatHappens, "n_1 = n_1 * 0.1")
    testthat::expect_equal(interventions[[1]]$Trigger, "n_2 >= 5")
})

test_that("Two interventions cannot have the same ID (check_double_id)", {
    interventions <- list(
        list(ID           = "intOverA",
            Trigger       = "n_B >= 5",
            WhatHappens   = "n_A = n_A * 0.1",
            Repetitions   = 0,
            Periodicity   = Inf
        ),list(
            ID          = "intOverA",
            Trigger     = "(n_B >= 5) and (T > 4)",
            WhatHappens = "n_ = n_ * 0.4",
            Repetitions = Inf,
            Periodicity = Inf
        )   
    )
    testthat::expect_error(createInterventions(interventions, afd3), "Check the interventions, there are 2 or more that have same IDs")    
})

test_that("The attribute WhatHappens is correctly specified (check_what_happens)",{
    interventions <- list(
        list(ID           = "intOverA",
            Trigger       = "n_B >= 5",
            WhatHappens   = "n_A +1 = n_A * 0.1",
            Repetitions   = 0,
            Periodicity   = Inf
        )  
    )

    testthat::expect_error(createInterventions(interventions, afd3), "The specification of WhatHappens is wrong.\n It should be: 
        <genotype_to_apply_some_operation or total_population> = <some_operation>\n Exiting.")

    interventions <- list(
        list(ID           = "intOverA",
            Trigger       = "n_B >= 5",
            WhatHappens   = "n_A = n_A * 0.1 = 32",
            Repetitions   = 0,
            Periodicity   = Inf
        )  
    )

    testthat::expect_error(createInterventions(interventions, afd3), "The specification of WhatHappens is wrong.\n It should be: 
        <genotype_to_apply_some_operation or total_population> = <some_operation>\n Exiting.")

    interventions <- list(
        list(ID           = "intOverA",
            Trigger       = "n_B >= 5",
            WhatHappens   = "= n_A * 0.1 = 32",
            Repetitions   = 0,
            Periodicity   = Inf
        )  
    )

    testthat::expect_error(createInterventions(interventions, afd3), "The specification of WhatHappens is wrong.\n It should be: 
        <genotype_to_apply_some_operation or total_population> = <some_operation>\n Exiting.")
})

test_that("The user cannot create population in an intervention",{
    # in this test, the main goal is to create a scenario where 
    # the whathappens is wrong, and creates population

    list_of_interventions <- list(
        list(ID           = "intOverMultiplicatesA",
            Trigger       = "n_B >= 5",
            WhatHappens   = "n_A = n_A * 2",
            Repetitions   = 0,
            Periodicity   = Inf
        )  
    )

    # we force the A genotype to not have mutationrate of 1 to avoid unexpected messages.
    fa1 <- data.frame(Genotype = c("A", "B"),
                    Fitness = c("1.1 + 0*n_A",
                                "1.002"))

    afd3 <- allFitnessEffects(genotFitness = fa1,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

    interventions <- createInterventions(list_of_interventions, afd3)

    testthat::expect_output(oncoSimulIndiv(
                    afd3, 
                    model = "McFL",
                    mu = 1e-4,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    finalTime = 5.2,
                    sampleEvery = 0.01,
                    onlyCancer = FALSE,
                    interventions = interventions
                    ), , paste0("In intervention:", interventions[[1]]$ID, 
                        " with WhatHappens: ", interventions[[1]]$WhatHappens, 
                        ". You cannot intervene to generate more population."))

})


test_that("1.1 Drastically reducing A population (McFL) | Trigger dependending on T", {

    fa1 <- data.frame(Genotype = c("A", "B"),
                    Fitness = c("1.001 + (0*n_A)",
                                "1.002"))

    afd3 <- allFitnessEffects(genotFitness = fa1,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

    ep1 <- oncoSimulIndiv(
                    afd3, 
                    model = "McFL",
                    mu = 1e-4,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
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
        WhatHappens   = "n_A = n_A * 0.1",
        Repetitions   = 0,
        Periodicity   = Inf
    ))

    interventions <- createInterventions(interventions, afd3)

    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "McFL",
                    mu = 1e-4,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    sampleEvery = 0.01,
                    finalTime = 5.2,
                    onlyCancer = FALSE,
		            interventions = interventions
                    )

    # we plot the simulation when interventions are specified.
    plot(ep2, show = "genotypes", type = "line")

    # to select the interval of data that pops.by.time offers
    ep2$pops.by.time[490:510, 1:3]
})





test_that("1.2 Drastically reducing A population (Exp) | Trigger dependending on T", {

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
            WhatHappens   = "n_A = n_A * 0.01",
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

    # to select the interval of data that pops.by.time offers
    sfd3$pops.by.time[1:5, 1:5]
})






test_that("1.3 Intervening over total population (McFL) | Trigger depends on T", {
    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                    Fitness = c("1",
                    "1 + 0.2 * (n_B > 0)",
                    ".9 + 0.4 * (n_A > 0)"
                    ))
    afd3 <- allFitnessEffects(genotFitness = gffd3,
                            frequencyDependentFitness = TRUE,
                            frequencyType = "abs")

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "McFL",
                            onlyCancer = FALSE,
                            finalTime = 200,
                            mu = 1e-4,
                            initSize = 5000,
                            sampleEvery = 0.01)

    plot(sfd3, show = "genotypes", type = "line")

    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "T > 40",
            WhatHappens   = "N = N * 0.2",
            Repetitions   = 2,   
            Periodicity   = 20
        )
    )

    interventions <- createInterventions(interventions, afd3)

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "McFL",
                            onlyCancer = FALSE,
                            finalTime = 200,
                            mu = 1e-4,
                            initSize = 5000,
                            sampleEvery = 0.01,
                            interventions = interventions)

    plot(sfd3, show = "genotypes", type = "line")

    # to select the interval of data that pops.by.time offers
    # First intervention
    sfd3$pops.by.time[4000:4020, 1:4]

    # Second intervention
    sfd3$pops.by.time[6000:6020, 1:4]

    # Third intervention
    sfd3$pops.by.time[7995:8010, 1:4]
})



test_that("1.4 Intervening over total population (Exp) | Trigger depends on T", {
    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                    Fitness = c("1",
                    "1 + 0.2 * (n_B * 0)",
                    ".9 + 0.4 * (n_B * 0)"
                    ))
    afd3 <- allFitnessEffects(genotFitness = gffd3,
                            frequencyDependentFitness = TRUE,
                            frequencyType = "abs")

    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "T > 40",
            WhatHappens   = "N = N * 0.2",
            Repetitions   = 2,   
            Periodicity   = 20
        )
    )

    interventions <- createInterventions(interventions, afd3)

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "Exp",
                            onlyCancer = FALSE,
                            finalTime = 200,
                            mu = 1e-4,
                            initSize = 5000,
                            sampleEvery = 0.01)

    plot(sfd3, show = "genotypes", type = "line")

    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "T > 50",
            WhatHappens   = "N = N * 0.5",
            Repetitions   = 1,   
            Periodicity   = 0.09
        )
    )

    interventions <- createInterventions(interventions, afd3)

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "Exp",
                            onlyCancer = FALSE,
                            finalTime = 200,
                            mu = 1e-4,
                            initSize = 5000,
                            sampleEvery = 0.01,
                            interventions = interventions)

    plot(sfd3, show = "genotypes", type = "line")

    # to select the interval of data that pops.by.time offers
    sfd3$pops.by.time[1:5, 1:5]
})




test_that("1.5 Intervening over total population (McFL) | Trigger depends on T", {
    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                    Fitness = c("1",
                    "1 + 0.2 * (n_B * 0)",
                    ".9 + 0.4 * (n_B * 0)"
                    ))
    afd3 <- allFitnessEffects(genotFitness = gffd3,
                            frequencyDependentFitness = TRUE,
                            frequencyType = "abs")

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "McFL",
                            onlyCancer = FALSE,
                            finalTime = 200,
                            mu = 1e-4,
                            initSize = 5000,
                            sampleEvery = 0.01)

    plot(sfd3, show = "genotypes", type = "line")

    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "T > 40",
            WhatHappens   = "N = N * 0.4",
            Repetitions   = 1,   
            Periodicity   = 0.09
        )
    )

    interventions <- createInterventions(interventions, afd3)

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "McFL",
                            onlyCancer = FALSE,
                            finalTime = 200,
                            mu = 1e-4,
                            initSize = 5000,
                            sampleEvery = 0.01,
                            interventions = interventions)

    plot(sfd3, show = "genotypes", type = "line")

    # to select the interval of data that pops.by.time offers
    sfd3$pops.by.time[1:5, 1:5]
})





test_that("1.6 Drastically reducing a high-fitness genotype population (Exp) | Trigger depends on T", {
    df3x <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                       Fitness = c("0*n_",
                                    "1.5",
                                    "1.0012",
                                    "1.002",
                                    "1.003",
                                    "1.004"))

    afd3 <- allFitnessEffects(genotFitness = df3x,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "Exp",
                    mu = 1e-4,
                    sampleEvery = 0.01,
                    initSize = c(20000, 20000, 300),
                    initMutant = c("A", "B", "C"),
                    finalTime = 200,
                    onlyCancer = FALSE
                    )

    # we plot the simulation when no interventions are specified.
    plot(ep2, show = "genotypes", type = "line")

    interventions <- list(
    list(ID           = "intOverB",
        Trigger       = "(T >= 10)",
        WhatHappens   = "n_B = n_B * 0.945",
        Repetitions   = Inf,
        Periodicity   = 0.1
    ))

    interventions <- createInterventions(interventions, afd3)

    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "Exp",
                    mu = 1e-4,
                    sampleEvery = 0.01,
                    initSize = c(20000, 20000, 300),
                    initMutant = c("A", "B", "C"),
                    finalTime = 200,
                    onlyCancer = FALSE,
                    interventions = interventions
                    )
    # we plot the simulation when interventions are specified.
    plot(ep2, show = "genotypes", type = "line")

    # to select the interval of data that pops.by.time offers
    sfd3$pops.by.time[1:5, 1:5]
})






test_that("1.7 Drastically reducing a high-fitness genotype population (Exp) | Trigger depends on n_*", {
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

    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "Exp",
                    mu = 1e-4,
                    sampleEvery = 0.01,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    finalTime = 200,
                    onlyCancer = FALSE
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
                    finalTime = 200,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    interventions = interventions
                    )
    plot(ep2, show = "genotypes", type = "line")

    # to select the interval of data that pops.by.time offers
    sfd3$pops.by.time[1:5, 1:5]
})






test_that("1.8 Intervening over 4 genotypes both over specific genotype and total population (McFL) | Trigger depends on N",{
    df3x <- data.frame(Genotype = c("A", "B", "C", "D"),
                      Fitness = c("1",
                                  "1.0001 + (0.0001 * n_A)",
                                  "1.1",
                                  "1.09"))

    afd3 <- allFitnessEffects(genotFitness = df3x,
                            frequencyDependentFitness = TRUE, frequencyType = "abs")

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "McFL",
                            onlyCancer = FALSE,
                            finalTime = 5.2,
                            mu = 1e-4,
                            sampleEvery = 0.01,
                            keepPhylog = FALSE,
                            seed = NULL,
                            errorHitMaxTries = FALSE,
                            errorHitWallTime = FALSE,
                            initMutant = c("A", "B", "C", "D"),
                            initSize = c(20000, 20000, 30, 10))

    plot(sfd3, show = "genotypes", type = "line")

    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "N >= 45000",
            WhatHappens   = "N = N * 0.5",
            Repetitions   = 1,   ## This will be translated to MAX_INT
            Periodicity   = 0.09
        ),
        list(
            ID            = "intOverA",
            Trigger       = "(n_B > n_A) and (n_A > 10) and (T > 30)",
            WhatHappens   = "n_A = 0.97 * n_A",
            Repetitions   = 20,   ## This will be translated to MAX_INT
            Periodicity   = 10
        )
    )

    interventions = createInterventions(interventions, afd3)

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "McFL",
                            onlyCancer = FALSE,
                            finalTime = 200,
                            mu = 1e-4,
                            sampleEvery = 0.01,
                            keepPhylog = FALSE,
                            seed = NULL,
                            errorHitMaxTries = FALSE,
                            errorHitWallTime = FALSE,
                            initMutant = c("A", "B", "C", "D"),
                            initSize = c(20000, 20000, 30, 10),
                            interventions = interventions)

    plot(sfd3, show = "genotypes", type = "line")

    # to select the interval of data that pops.by.time offers
    sfd3$pops.by.time[1:5, 1:5]
})





test_that("1.9 Intervening over 4 genotypes both over specific and total population (Exp) | Trigger depends on T", {

    df3x <- data.frame(Genotype = c("A", "B", "C", "D"),
                      Fitness = c("1",
                                  "0.93 + (0.0001 * n_A)",
                                  "1.4 - (0.000031 * n_B)",
                                  "1.09 + 0.2 * (n_B > n_A)"))

    afd3 <- allFitnessEffects(genotFitness = df3x,
                            frequencyDependentFitness = TRUE, frequencyType = "abs")

    evalAllGenotypes(afd3, spPopSizes = c(A = 2500, B = 2000, C = 5500, D = 700))

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "Exp",
                            onlyCancer = FALSE,
                            finalTime = 5.2,
                            mu = 1e-4,
                            sampleEvery = 0.01,
                            keepPhylog = FALSE,
                            seed = NULL,
                            errorHitMaxTries = FALSE,
                            errorHitWallTime = FALSE,
                            initMutant = c("A", "B", "C", "D"),
                            initSize = c(2500, 2000, 3000, 1000))

    plot(sfd3, show = "genotypes", type = "line")

    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "T > 3",
            WhatHappens   = "N = N * 0.5",
            Repetitions   = 1,   ## This will be translated to MAX_INT
            Periodicity   = 0.09
        ),
        list(
            ID            = "intOverA",
            Trigger       = "(N == 2000) and (n_A > n_B)",
            WhatHappens   = "n_A = 0.92 * n_A",
            Repetitions   = 20,   ## This will be translated to MAX_INT
            Periodicity   = 0.3
        )
    )

    interventions = createInterventions(interventions, afd3)

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "Exp",
                            onlyCancer = FALSE,
                            finalTime = 5.2,
                            mu = 1e-4,
                            sampleEvery = 0.01,
                            keepPhylog = FALSE,
                            seed = NULL,
                            errorHitMaxTries = FALSE,
                            errorHitWallTime = FALSE,
                            initMutant = c("A", "B", "C", "D"),
                            initSize = c(2500, 2000, 3000, 1000),
                            interventions = interventions)

    plot(sfd3, show = "genotypes", type = "line")

    # to select the interval of data that pops.by.time offers
    sfd3$pops.by.time[1:5, 1:5]
})

test_that("1.10 Intervention with Periodicity = Inf, should not execute affect the simulation", {
    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                    Fitness = c("1",
                    "1 + 0.2 * (n_B * 0)",
                    ".9 + 0.4 * (n_B * 0)"
                    ))
    afd3 <- allFitnessEffects(genotFitness = gffd3,
                            frequencyDependentFitness = TRUE,
                            frequencyType = "abs")

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "McFL",
                            onlyCancer = FALSE,
                            finalTime = 200,
                            mu = 1e-4,
                            initSize = 5000,
                            sampleEvery = 0.01)

    plot(sfd3, show = "genotypes", type = "line")

    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "T > 1",
            WhatHappens   = "n_ = n_ * 0.5",
            Repetitions   = Inf,   ## This will be translated to MAX_INT
            Periodicity   = Inf
        )
    )

    interventions <- createInterventions(interventions, afd3)

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "McFL",
                            onlyCancer = FALSE,
                            finalTime = 200,
                            mu = 1e-4,
                            initSize = 5000,
                            sampleEvery = 0.01)

    plot(sfd3, show = "genotypes", type = "line")

    # to select the interval of data that pops.by.time offers
    sfd3$pops.by.time[1:5, 1:5]
})

test_that("1.11 Intervening in the Rock-Paper-Scissors model y bacterial community by Kerr.", {
    crs <- function (a, b, c){
        data.frame(Genotype = c("WT", "C", "R"),
        Fitness = c(paste0("1 + ", a, " * n_R/N - ", b, " * n_C/N"),
        paste0("1 + ", b, " * n_/N - ", c, " * n_R/N"),
        paste0("1 + ", c, " * n_C/N - ", a, " * n_/N")
        ))
    }

    afcrs1 <- allFitnessEffects(genotFitness = crs(1, 1, 1),
    frequencyDependentFitness = TRUE,
    frequencyType = "abs")

    resultscrs1 <- oncoSimulIndiv(afcrs1,
                                model = "McFL",
                                onlyCancer = FALSE,
                                finalTime = 100,
                                mu = 1e-2,
                                initSize = 4000,
                                keepPhylog = TRUE,
                                seed = NULL,
                                errorHitMaxTries = FALSE,
                                errorHitWallTime = FALSE)

    lista_intervenciones = list(
        list(
            ID = "Intervención para Afectar a WT",
            Trigger = "n_C >= 500",
            WhatHappens = "n_C = n_C * 0.1",
            Periodicity = 3,
            Repetitions = Inf
        )
    )

    final_interventions = createInterventions(interventions = lista_intervenciones, genotFitness = afcrs1)

    resultscrs1 <- oncoSimulIndiv(afcrs1,
                                model = "McFL",
                                onlyCancer = FALSE,
                                finalTime = 100,
                                mu = 1e-2,
                                initSize = 4000,
                                keepPhylog = TRUE,
                                seed = NULL,
                                errorHitMaxTries = FALSE,
                                errorHitWallTime = FALSE,
                                interventions = final_interventions)

    plot(resultscrs1, show = "genotypes", type = "line", cex.lab=1.1,
    las = 1)
})

cat(paste("\n Ending interventions tests", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)