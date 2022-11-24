inittime <- Sys.time()
cat(paste("\n Starting interventions tests", date(), "\n"))

test_that("1. A intervention is created correctly",{
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

    testthat::expect_output(createInterventions(interventions, afd3),
                            "Checking intervention: intOverA" )
    interventions <- createInterventions(interventions, afd3)

    ## we check that the transformation of the WhatHappens and Trigger atributte is correct
    testthat::expect_equal(interventions[[1]]$WhatHappens, "n_1 = n_1 * 0.1")
    testthat::expect_equal(interventions[[1]]$Trigger, "n_2 >= 5")
})

test_that("2. Two interventions cannot have the same ID (check_double_id)", {
    interventions <- list(
        list(ID           = "intOverA",
            Trigger       = "n_B >= 5",
            WhatHappens   = "n_A = n_A * 0.1",
            Repetitions   = 0,
            Periodicity   = Inf
        ),list(
            ID          = "intOverA",
            Trigger     = "(n_B >= 5) and (T > 4)",
            WhatHappens = "n_A = n_A * 0.4",
            Repetitions = Inf,
            Periodicity = Inf
        )
    )
    testthat::expect_error(createInterventions(interventions, afd3),
                           "Check the interventions, there are 2 or more that have same IDs")
})

test_that("3. The attribute WhatHappens is correctly specified (check_what_happens)",{
    interventions <- list(
        list(ID           = "intOverA",
            Trigger       = "n_B >= 5",
            WhatHappens   = "n_A +1 = n_A * 0.1",
            Repetitions   = 0,
            Periodicity   = Inf
        )
    )

    testthat::expect_error(createInterventions(interventions, afd3),
                           "The specification of WhatHappens is wrong.\n It should be:")

    interventions <- list(
        list(ID           = "intOverA",
            Trigger       = "n_B >= 5",
            WhatHappens   = "n_A = n_A * 0.1 = 32",
            Repetitions   = 0,
            Periodicity   = Inf
        )
    )

    testthat::expect_error(createInterventions(interventions, afd3),
                           "The specification of WhatHappens is wrong.\n It should be:")

    interventions <- list(
        list(ID           = "intOverA",
            Trigger       = "n_B >= 5",
            WhatHappens   = "= n_A * 0.1 = 32",
            Repetitions   = 0,
            Periodicity   = Inf
        )
    )

    testthat::expect_error(createInterventions(interventions, afd3),
                           "The specification of WhatHappens is wrong.\n It should be:")
})

test_that("4. The user cannot create population in an intervention",{
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


test_that("5. Drastically reducing A-genotype population (McFL) | Trigger dependending on T", {
    seed <- round(runif(1, 1, 1e9))
    
    fa1 <- data.frame(Genotype = c("A", "B"),
                    Fitness = c("1.001 + (0*n_A)",
                                "1.002"))

    afd3 <- allFitnessEffects(genotFitness = fa1,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

    set.seed(seed)
    ## run the simulation without interventions
    ep1 <- oncoSimulIndiv(
        afd3,
        model = "McFL",
        mu = 1e-4,
        initSize = c(20000, 20000),
        initMutant = c("A", "B"),
        finalTime = 5.2,
        sampleEvery = 0.025,
        onlyCancer = FALSE
    )


    ## now we especify intervention to drastically reduce A population

    interv_fract <- runif(1, 0.01, 0.31)
    interventions <- list(
    list(ID           = "intOverA",
         Trigger       = "(T >= 3)",
         WhatHappens   =  paste0("n_A = n_A * ", interv_fract),
         Repetitions   = 0,
         Periodicity   = Inf
    ))

    interventions <- createInterventions(interventions, afd3)

    set.seed(seed)
    ## run the simulation WITH interventions
    ep2 <- oncoSimulIndiv(
        afd3,
        model = "McFL",
        mu = 1e-4,
        initSize = c(20000, 20000),
        initMutant = c("A", "B"),
        sampleEvery = 0.025,
        finalTime = 5.2,
        onlyCancer = FALSE,
        interventions = interventions
    )

    index <- which(ep2$pops.by.time[,1] %in% ep2$other$interventionTimes)

    expect_equal(ep2$pops.by.time[index, 2],
                 floor(interv_fract * ep1$pops.by.time[index, 2]))

})





test_that("6. Drastically reducing A population (Exp) | Trigger dependending on T", {

    seed <- round(runif(1, 1, 1e9))

    fa1 <- data.frame(Genotype = c("WT", "A", "B"),
                      Fitness = c("0 * n_", ## we need an expression
                                  "1.7",
                                  "1.1"))

    afd3 <- allFitnessEffects(genotFitness = fa1,
                              frequencyDependentFitness = TRUE,
                              frequencyType = "abs")

    set.seed(seed)
    ep1 <- oncoSimulIndiv(
        afd3,
        model = "Exp",
        mu = 1e-4,
        sampleEvery = 1, # 0.001,
        initSize = c(20000, 20000),
        initMutant = c("A", "B"),
        finalTime = 10, # 5.2,
        onlyCancer = FALSE
    )


    ## now we especify intervention to drastically reduce A population
    interv_fract <- runif(1, 0.01, 0.31)
    interventions <- list(
        list(ID           = "intOverA",
             Trigger       = "(T >= 5)",
             WhatHappens   = paste0("n_A = n_A * ", interv_fract), ## 0.1
             Repetitions   = 0,
             Periodicity   = Inf
             )
    )

    interventions <- createInterventions(interventions, afd3)

    set.seed(seed)
    ep2 <- oncoSimulIndiv(
        afd3,
        model = "Exp",
        mu = 1e-4,
        sampleEvery = 1, # 0.001,
        initSize = c(20000, 20000),
        initMutant = c("A", "B"),
        finalTime = 10, # 5.2,
        onlyCancer = FALSE,
        interventions = interventions
    )

    index <- which(ep2$pops.by.time[,1] %in% ep2$other$interventionTimes)

    
    last <- nrow(ep1$pops.by.time)
    ## when we do not intervene population of A will be bigger than B,
    ## since it has better fitness
    testthat::expect_gt(ep1$pops.by.time[last, 2],
                        ep1$pops.by.time[last, 3])

    ## once we intervene we test that the value of the population of A
    ## is quite lower once the intervention is made
    testthat::expect_gt(ep2$pops.by.time[index - 1, 2],
                        ep2$pops.by.time[index, 2])

    ## since in the first simulation we do not intervene,
    ## the population should be greater that
    ## when we intervene
    testthat::expect_lt(ep2$pops.by.time[index, 2],
                        ep1$pops.by.time[index, 2])

    ## And this is a much stronger test: the population
    ## of the intervened is exactly the intervention fraction
    ## of the non-intervened

    expect_equal(ep2$pops.by.time[index, 2],
                 floor(interv_fract * ep1$pops.by.time[index, 2]))

})


test_that("7. Intervening over total population (McFL) | Trigger depends on T", {
    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                        Fitness = c("1",
                                    "1 + 0.2 * (n_B > 0)",
                                    ".9 + 0.4 * (n_A > 0)"
                                    ))
    aafd3 <- allFitnessEffects(genotFitness = gffd3,
                               frequencyDependentFitness = TRUE,
                               frequencyType = "abs")


    interv_fract <- round(runif(1, 0.4, 0.8), 2)

    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "T > 5",
            WhatHappens   = paste0("N = N * ", interv_fract),
            Repetitions   = 2,
            Periodicity   = 5
        )
    )

    interventions <- createInterventions(interventions, aafd3)
    
    sfd3 <- oncoSimulIndiv(aafd3,
                           model = "McFL",
                           onlyCancer = FALSE,
                           finalTime = 16,
                           mu = 1e-4,
                           initSize = 2e4, 
                           sampleEvery = 0.025, 
                           interventions = interventions,
                           detectionSize = NA
                           )

    ## it may happen that, in some simulations, the population collapses, in that case,
    ## pops by time is null, and cannot be checked

    if (!is.null(sfd3$pops.by.time)) {
        indexes <- which(sfd3$pops.by.time[,1] %in% sfd3$other$interventionTimes)
        total_before <- rowSums(sfd3$pops.by.time[indexes - 1, -1])
        total_after <-  rowSums(sfd3$pops.by.time[indexes, -1])
        reduction <- round(total_after/total_before, 2)
        print(reduction)
        expect_equal(reduction, rep(interv_fract, 3))
    }
    
})


test_that("8. Intervening over total population (Exp) | Trigger depends on T", {
    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                        Fitness = c("1",
                                    "1 + 0.2 * (n_B > 0)",
                                    ".9 + 0.4 * (n_A > 0)"
                                    ))
    aafd3 <- allFitnessEffects(genotFitness = gffd3,
                               frequencyDependentFitness = TRUE,
                               frequencyType = "abs")


    interv_fract <- round(runif(1, 0.4, 0.8), 2)

    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "T > 10",
            WhatHappens   = paste0("N = N * ", interv_fract),
            Repetitions   = 2,
            Periodicity   = 10
        )
    )

    interventions <- createInterventions(interventions, aafd3)
    
    sfd3 <- oncoSimulIndiv(aafd3,
                           model = "Exp",
                           onlyCancer = FALSE,
                           finalTime = 32,
                           mu = 1e-4,
                           initSize = 2e5, ## 20000,
                           sampleEvery = 1, ## 0.001,
                           interventions = interventions,
                           detectionSize = NA
                           )

    ## it may happen that, in some simulations, the population collapses, in that case,
    ## pops by time is null, and cannot be checked

    if (!is.null(sfd3$pops.by.time)) {
        indexes <- which(sfd3$pops.by.time[,1] %in% sfd3$other$interventionTimes)
        total_before <- rowSums(sfd3$pops.by.time[indexes - 1, -1])
        total_after <-  rowSums(sfd3$pops.by.time[indexes, -1])
        reduction <- round(total_after/total_before, 2)
        expect_equal(reduction, rep(interv_fract, 3))
    }
    
})

## test 9 and 10 found in test.Z-intervention.R

test_that("11. Intervening over 4 genotypes both over specific genotype and total population (Exp) | Trigger depends on N",{
    df3x <- data.frame(Genotype = c("A", "B", "C", "D", "E"),
                       Fitness = c("1",
                                   "1.01 + (0 * n_A)",
                                   "1.1",
                                   "1.19",
                                   "1.17"))

    afd3 <- allFitnessEffects(genotFitness = df3x,
                              frequencyDependentFitness = TRUE, frequencyType = "abs")

    interventions = list(
        list(
            ID            = "intOverA",
            Trigger       = "T >= 1.2",
            WhatHappens   = "n_A = n_A * 0.2",
            Repetitions   = 0,
            Periodicity   = Inf
        ),
        list(
            ID            = "intOverB",
            Trigger       = "T >= 2.2",
            WhatHappens   = "n_B = n_B * 0.2",
            Repetitions   = 0,
            Periodicity   = Inf
        ),
        list(
            ID            = "intOverC",
            Trigger       = "T >= 3.2",
            WhatHappens   = "n_C = n_C * 0.1",
            Repetitions   = 0,
            Periodicity   = Inf
        ),
        list(
            ID            = "intOverD",
            Trigger       = "T >= 4.2",
            WhatHappens   = "n_D = n_D * 0.01",
            Repetitions   = 0,
            Periodicity   = Inf
        )
    )

    interventions = createInterventions(interventions, afd3)
    seed <- round(runif(1, 1, 1e9))
    set.seed(seed)
    sfd3_with_ints <- oncoSimulIndiv( afd3,
                                     model = "Exp",
                                     onlyCancer = FALSE,
                                     finalTime = 5,
                                     mu = 1e-4,
                                     sampleEvery = 1,
                                     keepPhylog = FALSE,
                                     errorHitMaxTries = FALSE,
                                     errorHitWallTime = FALSE,
                                     initMutant = c("A", "B", "C", "D", "E"),
                                     initSize = c(20000, 20000, 3000, 10000, 200),
                                     interventions = interventions)
    set.seed(seed)
    sfd3_without_ints <- oncoSimulIndiv( afd3,
                                        model = "Exp",
                                        onlyCancer = FALSE,
                                        finalTime = 5,
                                        mu = 1e-4,
                                        sampleEvery = 1,
                                        keepPhylog = FALSE,
                                        errorHitMaxTries = FALSE,
                                        errorHitWallTime = FALSE,
                                        initMutant = c("A", "B", "C", "D", "E"),
                                        initSize = c(20000, 20000, 3000, 10000, 200))
    
    indexes <- which(sfd3_with_ints$pops.by.time[,1] %in%
                     sfd3_with_ints$other$interventionTimes)


    ## Same logic as above: compare populations without and without intervention
    ## But after the first difference, trajectories could differ from
    ## random number divergences
    expect_equal(sfd3_with_ints$pops.by.time[indexes[1], 2],
                 floor(.2 * sfd3_without_ints$pops.by.time[indexes[1], 2]))
    expect_lt(sfd3_with_ints$pops.by.time[indexes[2], 3],
              sfd3_without_ints$pops.by.time[indexes[2], 3])
    expect_lt(sfd3_with_ints$pops.by.time[indexes[3], 4],
              sfd3_without_ints$pops.by.time[indexes[3], 4])
    expect_lt(sfd3_with_ints$pops.by.time[indexes[4], 5],
              sfd3_without_ints$pops.by.time[indexes[4], 5])

    
})


test_that("12. Intervening over 4 genotypes both over specific genotype and total population (McFL) | Trigger depends on N",{
    ## Note same logic as in 11, because death rate in McFL decreases when we decrease
    ## population size. So intervene at exact same time on all
    ## and test exact proportions for all
    df3x <- data.frame(Genotype = c("A", "B", "C", "D", "E"),
                       Fitness = c("1",
                                   "1.01 + (0 * n_A)",
                                   "1.1",
                                   "1.09",
                                   "1.07"))

    afd3 <- allFitnessEffects(genotFitness = df3x,
                              frequencyDependentFitness = TRUE, frequencyType = "abs")

    interventions = list(
        list(
            ID            = "intOverA",
            Trigger       = "T >= 1.2",
            WhatHappens   = "n_A = n_A * 0.2",
            Repetitions   = 0,
            Periodicity   = Inf
        ),
        list(
            ID            = "intOverB",
            Trigger       = "T >= 1.2",
            WhatHappens   = "n_B = n_B * 0.42",
            Repetitions   = 0,
            Periodicity   = Inf
        ),
        list(
            ID            = "intOverC",
            Trigger       = "T >= 1.2",
            WhatHappens   = "n_C = n_C * 0.35",
            Repetitions   = 0,
            Periodicity   = Inf
        ),
        list(
            ID            = "intOverD",
            Trigger       = "T >= 1.2",
            WhatHappens   = "n_D = n_D * 0.125",
            Repetitions   = 0,
            Periodicity   = Inf
        )
    )

    interventions = createInterventions(interventions, afd3)
    seed <- round(runif(1, 1, 1e9))
    set.seed(seed)
    sfd3_with_ints <- oncoSimulIndiv( afd3,
                                     model = "McFL",
                                     onlyCancer = FALSE,
                                     finalTime = 1.3,
                                     mu = 1e-4,
                                     sampleEvery = 0.025,
                                     keepPhylog = FALSE,
                                     errorHitMaxTries = FALSE,
                                     errorHitWallTime = FALSE,
                                     initMutant = c("A", "B", "C", "D", "E"),
                                     initSize = c(20000, 20000, 30, 10, 200),
                                     interventions = interventions)
    set.seed(seed)
    sfd3_without_ints <- oncoSimulIndiv( afd3,
                                        model = "McFL",
                                        onlyCancer = FALSE,
                                        finalTime = 1.3,
                                        mu = 1e-4,
                                        sampleEvery = 0.025,
                                        keepPhylog = FALSE,
                                        errorHitMaxTries = FALSE,
                                        errorHitWallTime = FALSE,
                                        initMutant = c("A", "B", "C", "D", "E"),
                                        initSize = c(20000, 20000, 30, 10, 200))
    
    indexes <- which(sfd3_with_ints$pops.by.time[,1] %in%
                     sfd3_with_ints$other$interventionTimes)


    expect_equal(sfd3_with_ints$pops.by.time[indexes, 2],
                 floor(.2 * sfd3_without_ints$pops.by.time[indexes, 2]))
    expect_equal(sfd3_with_ints$pops.by.time[indexes, 3],
                 floor(.42 * sfd3_without_ints$pops.by.time[indexes, 3]))
    expect_equal(sfd3_with_ints$pops.by.time[indexes, 4],
                 floor(.35 * sfd3_without_ints$pops.by.time[indexes, 4]))
    expect_equal(sfd3_with_ints$pops.by.time[indexes, 5],
                 floor(.125 * sfd3_without_ints$pops.by.time[indexes, 5]))
    
})




## ## This is  not testing interventions, but interventions
## ## plus other things in a complex to reason about model
## test_that("13. Intervening in the Rock-Paper-Scissors model for bacterial community by Kerr.", {
##     crs <- function (a, b, c){
##         data.frame(Genotype = c("WT", "C", "R"),
##         Fitness = c(paste0("1 + ", a, " * n_R/N - ", b, " * n_C/N"),
##         paste0("1 + ", b, " * n_/N - ", c, " * n_R/N"),
##         paste0("1 + ", c, " * n_C/N - ", a, " * n_/N")
##         ))
##     }

##     afcrs1 <- allFitnessEffects(genotFitness = crs(1, 1, 1),
##     frequencyDependentFitness = TRUE,
##     frequencyType = "abs")

##     lista_intervenciones = list(
##         list(
##             ID = "Bothering R strain, by reducing C",
##             Trigger = "n_C >= 500",
##             WhatHappens = "n_C = n_C * 0.1",
##             Periodicity = 3,
##             Repetitions = Inf
##         )
##     )

##     final_interventions = createInterventions(interventions = lista_intervenciones, genotFitness = afcrs1)

##     resultscrs1 <- oncoSimulIndiv(afcrs1,
##                                 model = "McFL",
##                                 onlyCancer = FALSE,
##                                 finalTime = 100,
##                                 mu = 1e-2,
##                                 initSize = 4000,
##                                 keepPhylog = FALSE,
##                                 seed = NULL,
##                                 errorHitMaxTries = FALSE,
##                                 errorHitWallTime = FALSE,
##                                 interventions = final_interventions)

##     # by reducing C, R wont spread in the population. This will mean that, with the apropiate
##     # periodicty in the intervention, C will never surpass WT. We try to check this in these tests.
##     i <- 1
##     while(i <= nrow(resultscrs1$pops.by.time)){
##         cur_nWT = resultscrs1$pops.by.time[i, 2:2]
##         cur_nR = resultscrs1$pops.by.time[i, 4:4]
##         testthat::expect_gt(cur_nWT, cur_nR)
##         i <- i + 1
##     }
## })

test_that("14. Intervening over total population (Exp) | Trigger depends on user variables", {
    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                        Fitness = c("1",
                                    "1 + 0.2 * (n_B > 0)",
                                    ".9 + 0.4 * (n_A > 0)"
                                    ))
    afd3 <- allFitnessEffects(genotFitness = gffd3,
                              frequencyDependentFitness = TRUE,
                              frequencyType = "abs")

    userVars <- list(
        list(Name = "user_var_1",
             Value = 0
             )
    )

    userVars <- createUserVars(userVars)

    rules <- list(
        list(ID = "rule_1",
             Condition = "T >= 10",
             Action = "user_var_1 = 1"
             ),list(ID = "rule_2",
                    Condition = "T >= 20",
                    Action = "user_var_1 = 2"
                    ),list(ID = "rule_3",
                           Condition = "T >= 30",
                           Action = "user_var_1 = 3"
                           )
    )

    rules <- createRules(rules, afd3)
    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "user_var_1 = 1",
            WhatHappens   = "N = N * 0.7",
            Repetitions   = 0,
            Periodicity   = Inf
        ),list(
              ID            = "intOverTotPop2",
              Trigger       = "user_var_1 = 2",
              WhatHappens   = "N = N * 0.2",
              Repetitions   = 0,
              Periodicity   = Inf
          ),list(
                ID            = "intOverTotPop3",
                Trigger       = "user_var_1 = 3",
                WhatHappens   = "N = N * 0.5",
                Repetitions   = 0,
                Periodicity   = Inf
            )
    )

    interventions <- createInterventions(interventions, afd3)

    sfd3 <- oncoSimulIndiv( afd3,
                           model = "Exp",
                           onlyCancer = FALSE,
                           finalTime = 35,
                           mu = 1e-4,
                           initSize = 5000,
                           sampleEvery = .025,
                           interventions = interventions,
                           userVars = userVars,
                           rules = rules,
                           ## FIXME
                           ## In Windows sometimes this takes forever
                           max.wall.time = 600)

    if (!is.null(sfd3$pops.by.time)) {
        indexes <- which(sfd3$pops.by.time[,1] %in% sfd3$other$interventionTimes)
        total_before <- rowSums(sfd3$pops.by.time[indexes - 1, -1])
        total_after <-  rowSums(sfd3$pops.by.time[indexes, -1])
        reduction <- round(total_after/total_before, 1)
        expect_equal(reduction, c(0.7, 0.2, 0.5))
    }
})


test_that("15. Intervening over total population (Exp) | WhatHappens uses user variables", {
    ## gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
    ##                 Fitness = c("1",
    ##                             "1 + 0.128 * (n_B > 0)",
    ##                             ".9 + 0.237 * (n_A > 0)"
    ##                             ))
    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                        Fitness = c("1",
                                    "1.21 + 0 * (n_B > 0)",
                                    "1.45 + 0 * (n_A > 0)"
                                    ))
    afd3 <- allFitnessEffects(genotFitness = gffd3,
                              frequencyDependentFitness = TRUE,
                              frequencyType = "abs")

    userVars <- list(
        list(Name = "user_var_1",
             Value = 0
             )
    )

    userVars <- createUserVars(userVars)

    rules <- list(
        list(ID = "rule_1",
             Condition = "T >= 10",
             Action = "user_var_1 = 0.5"
             ),list(ID = "rule_2",
                    Condition = "T >= 20",
                    Action = "user_var_1 = 0.8"
                    ),list(ID = "rule_3",
                           Condition = "T >= 30",
                           Action = "user_var_1 = 0.7"
                           )
    )

    rules <- createRules(rules, afd3)
    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "T >= 15",
            WhatHappens   = "N = N * user_var_1",
            Repetitions   = 2,
            Periodicity   = 10
        )
    )

    interventions <- createInterventions(interventions, afd3)

    sfd3 <- oncoSimulIndiv( afd3,
                           model = "Exp",
                           onlyCancer = FALSE,
                           finalTime = 40,
                           mu = 1e-4,
                           initSize = 5000,
                           sampleEvery = 1,
                           interventions = interventions,
                           userVars = userVars,
                           rules = rules,
                           detectionSize = NA,
                           ## FIXME: again, in Windows this sometimes takes a very long time
                           max.wall.time = 600)


    if (!is.null(sfd3$pops.by.time)) {
        indexes <- which(sfd3$pops.by.time[,1] %in% sfd3$other$interventionTimes)
        total_before <- rowSums(sfd3$pops.by.time[indexes - 1, -1])
        total_after <-  rowSums(sfd3$pops.by.time[indexes, -1])
        reduction <- round(total_after/total_before, 1)
        expect_equal(reduction, c(0.5, 0.8, 0.7))
    }
    
})




## FIXME: write tests like this
## test_that("x1. Total pop intervention cannot increase population size", {

##     gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
##                         Fitness = c("1",
##                                     "1.1 + 0 * (n_B > 0)",
##                                     "1.3 + 0 * (n_A > 0)"
##                                     ))
##     aafd3 <- allFitnessEffects(genotFitness = gffd3,
##                                frequencyDependentFitness = TRUE,
##                                frequencyType = "abs")


##     interventions = list(
##         list(
##             ID            = "intOverTotPop",
##             Trigger       = "T > 5",
##             WhatHappens   = "N = N * 2",
##             Repetitions   = 2,
##             Periodicity   = 5
##         )
##     )
##     interventions <- createInterventions(interventions, aafd3)
##     oncoSimulIndiv(aafd3,
##                    model = "Exp",
##                    onlyCancer = FALSE,
##                    finalTime = 16,
##                    mu = 1e-4,
##                    initSize = 2e4, 
##                    sampleEvery = 0.025, 
##                    interventions = interventions,
##                    detectionSize = NA
##                    )
        
        
## })




## FIXME: additional checks
##  - that when the interventions take place is when they should
##    for example with complex triggers

## Something like this, with a complex trigger and a complex what happens

## test_that("blablba. Drastically reducing a high-fitness genotype population (Exp) | Trigger depends on n_*", {

##     df3x <- data.frame(Genotype = c("WT", "B", "R"),
##                        Fitness = c("1 + (n_ * 0)",
##                                    "1.1",
##                                    "1.3"))

##     afd3 <- allFitnessEffects(genotFitness = df3x,
##                               frequencyDependentFitness = TRUE,
##                               frequencyType = "abs")

##     interventions <- list(
##         list(
##             ID          = "intOverBAffectsR",
##             Trigger     = "((T >= 20) and (T <= 70)) and (n_B >= 200)",
##             WhatHappens = "n_R = 300",
##             Repetitions = Inf,
##             Periodicity = 5
##         ),
##         list(
##             ID          = "intOverTotPop",
##             Trigger     = "(T >= 80) and (T <= 85) and (n_B >= 40)",
##             WhatHappens = "N = 2000",
##             Repetitions = Inf,
##             Periodicity = 2
##         )
##     )

##     interventions <- createInterventions(interventions, afd3)

##     ep2 <- oncoSimulIndiv(
##         afd3,
##         model = "Exp",
##         mu = 1e-4,
##         sampleEvery = 1,
##         initSize = c(5000, 10, 300),
##         initMutant = c("WT", "B", "R"),
##         finalTime = 100,
##         onlyCancer = FALSE,
##         interventions = interventions)

##     ## FIXME: now, verify
## })





cat(paste("\n Ending interventions tests", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)






################################# This used to be in test.Z-intervention.R


## inittime <- Sys.time()
## cat(paste("\n Starting Z-interventions tests", date(), "\n"))

## ## FIXME
## ## These two tests are extremely computationally intensive,
## ## and I think could be better tested otherwise.
## ## And this test has many non-idiomatic R constructs
## ## and it could probably run in a 1/10 of the time

## ## And these are not really good tests: they are deterministic
## ## and focus on irrelevant, accessory details.

## ## Disabled for now

## test_that("1. Drastically reducing a high-fitness genotype population (McFL) | Trigger depends on T and n_*", {
##     set.seed(1)
##     df3x <- data.frame(Genotype = c("WT", "B", "R"),
##                        Fitness = c("1 + (n_ * 0)",
##                                    "1.5",
##                                    "1.003 + 0.002 * (n_B > 120)"))

##     afd3 <- allFitnessEffects(genotFitness = df3x,
##                               frequencyDependentFitness = TRUE,
##                               frequencyType = "abs")

##     ## FIXME: why such periodicity?
##     interventions <- list(
##         list(
##             ID          = "intOverBAffectsR",
##             Trigger     = "((T >= 20) and (T <= 70)) and (n_B >= 200)",
##             WhatHappens = "n_B = 100",
##             Repetitions = 50,
##             Periodicity = 0.001
##         ),
##         list(
##             ID          = "intOverTotPop",
##             Trigger     = "(T >= 80) and (T <= 85) and (n_B >= 40)",
##             WhatHappens = "n_B = 20",
##             Repetitions = Inf,
##             Periodicity = 0.001
##         )
##     )

##     interventions <- createInterventions(interventions, afd3)


##     ep2 <- oncoSimulIndiv(
##         afd3,
##         model = "McFLD",
##         mu = 1e-4,
##         sampleEvery = 0.001,
##         initSize = c(5000, 10, 300),
##         initMutant = c("WT", "B", "R"),
##         finalTime = 100,
##         onlyCancer = FALSE,
##         interventions = interventions,
##         ## FIXME: this test occasionally fails
##         ## as it goes > 200 s in Windows.
##         ## This should not be needed
##         max.wall.time = 600
##     )

##     ## Why the thresholds? 210 here and 40 below.

##     flag <- FALSE
##     i <- 20002
##     while(i <= 70001){
##         if(ep2$pops.by.time[i, 3:3] >= 210){
##             flag <- TRUE
##         }
##         i <- i + 1
##     }

##     testthat::expect_equal(flag, FALSE)

##                                         # then, between the time intervals, T >= 80 and T<=85
##                                         # we control that the B population
##     flag <- FALSE
##     i <- 80002

##     ## FIXME: why simulate to 100 time units if we only look up
##     ## to row 85000?
##     while(i <= 85000){
##         if(ep2$pops.by.time[i, 3:3] > 40){
##             flag <- TRUE
##         }
##         i <- i + 1
##     }

##     testthat::expect_equal(flag, FALSE)

##                                         # we plot the simulation when no interventions are specified.
##                                         #plot(ep2, show = "genotypes", type = "line")
## })





## test_that("2. Drastically reducing a high-fitness genotype population (Exp) | Trigger depends on n_*", {
##     set.seed(1)
##     df3x <- data.frame(Genotype = c("WT", "B", "R"),
##                        Fitness = c("1 + (n_ * 0)",
##                                    "1.5",
##                                    "1.003 + 0.002 * (n_B > 120)"))

##     afd3 <- allFitnessEffects(genotFitness = df3x,
##                               frequencyDependentFitness = TRUE,
##                               frequencyType = "abs")

##     interventions <- list(
##         list(
##             ID          = "intOverBAffectsR",
##             Trigger     = "((T >= 20) and (T <= 70)) and (n_B >= 200)",
##             WhatHappens = "n_B = 100",
##             Repetitions = 50,
##             Periodicity = 0.001
##         ),
##         list(
##             ID          = "intOverTotPop",
##             Trigger     = "(T >= 80) and (T <= 85) and (n_B >= 40)",
##             WhatHappens = "n_B = 20",
##             Repetitions = Inf,
##             Periodicity = 0.001
##         )
##     )

##     interventions <- createInterventions(interventions, afd3)

##     ep2 <- oncoSimulIndiv(
##         afd3,
##         model = "Exp",
##         mu = 1e-4,
##         sampleEvery = 0.001,
##         initSize = c(5000, 10, 300),
##         initMutant = c("WT", "B", "R"),
##         finalTime = 100,
##         onlyCancer = FALSE,
##         interventions = interventions,
##         ## FIXME: This huge wall time should not be necessary.
##         ## See above; this is slow as hell in Windows.
##         max.wall.time = 600)


##     flag <- FALSE
    ##     i <- 20002
    ##     while(i <= 70001) {
    ##         if(ep2$pops.by.time[i, 3:3] >= 210) {
    ##             flag <- TRUE
    ##         }
    ##         i <- i + 1
    ##     }
    ##     testthat::expect_equal(flag, FALSE)


    ##     ## then, between the time intervals, T >= 80 and T<=85
    ##     ## we control that the B population
    ##     flag <- FALSE
    ##     i <- 80002
    ##     while(i <= 85000) {
    ##         if(ep2$pops.by.time[i, 3:3] > 40){
    ##             flag <- TRUE
    ##         }
    ##         i <- i + 1
    ##     }

    ##     testthat::expect_equal(flag, FALSE)
    ## })

    ## cat(paste("\n Ending Z-interventions tests", date(), "\n"))
    ## cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
    ## rm(inittime)









