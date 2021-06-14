inittime <- Sys.time()
cat(paste("\n Starting interventions tests", date(), "\n"))
############################################################################################################################
############################################################################################################################
############################################################################################################################

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

    testthat::expect_output(createInterventions(interventions, afd3), "Checking intervention: intOverA" )
    interventions <- createInterventions(interventions, afd3)

    # we check that the transformation of the WhatHappens and Trigger atributte is correct
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
            WhatHappens = "n_ = n_ * 0.4",
            Repetitions = Inf,
            Periodicity = Inf
        )   
    )
    testthat::expect_error(createInterventions(interventions, afd3), "Check the interventions, there are 2 or more that have same IDs")    
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

    fa1 <- data.frame(Genotype = c("A", "B"),
                    Fitness = c("1.001 + (0*n_A)",
                                "1.002"))

    afd3 <- allFitnessEffects(genotFitness = fa1,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")
    # run the simulation without interventions
    ep1 <- oncoSimulIndiv(
                    afd3, 
                    model = "McFL",
                    mu = 1e-4,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    finalTime = 5.2,
                    sampleEvery = 0.001,
                    onlyCancer = FALSE
                    )

    # now we especify intervention to drastically reduce A population
    interventions <- list(
    list(ID           = "intOverA",
        Trigger       = "(T >= 5)",
        WhatHappens   = "n_A = n_A * 0.1",
        Repetitions   = 0,
        Periodicity   = Inf
    ))

    interventions <- createInterventions(interventions, afd3)

    # run the simulation WITH interventions
    ep2 <- oncoSimulIndiv(
                    afd3, 
                    model = "McFL",
                    mu = 1e-4,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    sampleEvery = 0.001,
                    finalTime = 5.2,
                    onlyCancer = FALSE,
		            interventions = interventions
                    )

    # test that the value of the population of A is quite lower once the intervention is made
    testthat::expect_gt(ep2$pops.by.time[4995:4995, 2:2], ep2$pops.by.time[5005:5005, 2:2])

    # since in the first simulation we do not intervene, the population should be greater that
    # when we intervene
    testthat::expect_lt(ep2$pops.by.time[5005:5005, 2:2], ep1$pops.by.time[5005:5005, 2:2])

})





test_that("6. Drastically reducing A population (Exp) | Trigger dependending on T", {

    fa1 <- data.frame(Genotype = c("WT", "A", "B"),
                    Fitness = c("n_*0",
                                "1.5",
                                "1"))

    afd3 <- allFitnessEffects(genotFitness = fa1,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

    ep1 <- oncoSimulIndiv(
                    afd3, 
                    model = "Exp",
                    mu = 1e-4,
                    sampleEvery = 0.001,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    finalTime = 5.2,
                    onlyCancer = FALSE
                    )

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
                    sampleEvery = 0.001,
                    initSize = c(20000, 20000),
                    initMutant = c("A", "B"),
                    finalTime = 5.2,
                    onlyCancer = FALSE,
		            interventions = interventions
                    )

    # when we do not intervene population of A will be bigger than B, since it has better fitness
    testthat::expect_gt(ep1$pops.by.time[5005:5005, 2:2], ep1$pops.by.time[5005:5005, 3:3])

    # once we intervene we test that the value of the population of A is quite lower once the intervention is made
    testthat::expect_gt(ep2$pops.by.time[4995:4995, 2:2], ep2$pops.by.time[5005:5005, 2:2])

    # since in the first simulation we do not intervene, the population should be greater that
    # when we intervene
    testthat::expect_lt(ep2$pops.by.time[5005:5005, 2:2], ep1$pops.by.time[5005:5005, 2:2])
})


test_that("7. Intervening over total population (McFL) | Trigger depends on T", {
    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                    Fitness = c("1",
                    "1 + 0.25 * (n_B > 0)",
                    ".9 + 0.4 * (n_A > 0)"
                    ))
    afd3 <- allFitnessEffects(genotFitness = gffd3,
                            frequencyDependentFitness = TRUE,
                            frequencyType = "abs")

    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "T > 40",
            WhatHappens   = "N = N * 0.6",
            Repetitions   = 2,   
            Periodicity   = 20
        )
    )

    interventions <- createInterventions(interventions, afd3)

    sfd3 <- oncoSimulIndiv( afd3,
                            model = "McFLD",
                            onlyCancer = FALSE,
                            finalTime = 100,
                            mu = 1e-4,
                            initSize = 5000,
                            sampleEvery = 0.001,
                            interventions = interventions)

    # we can check genotype by genotype that when an intervention ocurs, their population lowers
    # when T = 40:
        #Genotype WT
    if((sfd3$pops.by.time[39997:39997, 2:2] > 0) & (sfd3$pops.by.time[40003:40003, 2:2] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[39997:39997, 2:2]), (sfd3$pops.by.time[40003:40003, 2:2]))
    }
        #Genotype A
    if((sfd3$pops.by.time[39997:39997, 3:3]) > 0 & (sfd3$pops.by.time[40003:40003, 3:3] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[39997:39997, 3:3]), (sfd3$pops.by.time[40003:40003, 3:3]))
    }    #Genotype B
    if((sfd3$pops.by.time[39997:39997, 4:4]) > 0 & (sfd3$pops.by.time[40003:40003, 4:4] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[39997:39997, 4:4]), (sfd3$pops.by.time[40003:40003, 4:4]))
    }
    # when T = 60
        #Genotype WT
    if((sfd3$pops.by.time[59997:59997, 2:2] > 0) & (sfd3$pops.by.time[60003:60003, 2:2] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[59997:59997, 2:2]), (sfd3$pops.by.time[60003:60003, 2:2]))
    }
        #Genotype A
    if((sfd3$pops.by.time[59997:59997, 3:3] > 0) & (sfd3$pops.by.time[60003:60003, 3:3] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[59997:59997, 3:3]), (sfd3$pops.by.time[60003:60003, 3:3]))
    }
        #Genotype B
    if((sfd3$pops.by.time[59997:59997, 4:4] > 0) & (sfd3$pops.by.time[60003:60003, 4:4] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[59997:59997, 4:4]), (sfd3$pops.by.time[60003:60003, 4:4]))
    }
    # when T = 80
        #Genotype WT
    if((sfd3$pops.by.time[79997:79997, 2:2] > 0) & (sfd3$pops.by.time[80003:80003, 2:2] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[79997:79997, 2:2]), (sfd3$pops.by.time[80003:80003, 2:2]))
    }
        #Genotype A
    if((sfd3$pops.by.time[79997:79997, 3:3] > 0) & (sfd3$pops.by.time[80003:80003, 3:3] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[79997:79997, 3:3]), (sfd3$pops.by.time[80003:80003, 3:3]))
    }
        #Genotype B
    if((sfd3$pops.by.time[79997:79997, 4:4] > 0) & (sfd3$pops.by.time[80003:80003, 4:4] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[79997:79997, 4:4]), (sfd3$pops.by.time[80003:80003, 4:4]))
    }

    # it is expected from the interventions to reduce the population in 3 different times:
    # when T = 40
    tot_pop_t_40_int = (sfd3$pops.by.time[40003:40003, 2:2] + sfd3$pops.by.time[40003:40003, 3:3] 
                          + sfd3$pops.by.time[40003:40003, 4:4])
    tot_pop_t_before40_int = (sfd3$pops.by.time[39997:39997, 2:2] + sfd3$pops.by.time[39997:39997, 3:3] 
                                + sfd3$pops.by.time[39997:39997, 4:4])
    testthat::expect_gt(tot_pop_t_before40_int, tot_pop_t_40_int)
    # when T = 60
    tot_pop_t_60_int = (sfd3$pops.by.time[60003:60003, 2:2] + sfd3$pops.by.time[60003:60003, 3:3] 
                          + sfd3$pops.by.time[60003:60003, 4:4])
    tot_pop_t_before60_int = (sfd3$pops.by.time[59997:59997, 2:2] + sfd3$pops.by.time[59997:59997, 3:3] 
                                + sfd3$pops.by.time[59997:59997, 4:4])
    testthat::expect_gt(tot_pop_t_before60_int, tot_pop_t_60_int)
    # when T = 80
    tot_pop_t_80_int = (sfd3$pops.by.time[80003:80003, 2:2] + sfd3$pops.by.time[80003:80003, 3:3] 
                          + sfd3$pops.by.time[80003:80003, 4:4])
    tot_pop_t_before80_int = (sfd3$pops.by.time[79997:79997, 2:2] + sfd3$pops.by.time[79997:79997, 3:3] 
                            + sfd3$pops.by.time[79997:79997, 4:4])
    testthat::expect_gt(tot_pop_t_before80_int, tot_pop_t_80_int)

    # finally we check that the population is in range within the specified in WhatHappens
    # since the intervention reduces all the population to the 60% of its original value, we will check
    # that in fact the population reduces (aprox) in that percentage
    testthat::expect_gt(tot_pop_t_before40_int * 0.61, tot_pop_t_40_int)
    testthat::expect_gt(tot_pop_t_before60_int * 0.61, tot_pop_t_60_int)
    testthat::expect_gt(tot_pop_t_before80_int * 0.61, tot_pop_t_80_int)
    
})



test_that("8. Intervening over total population (Exp) | Trigger depends on T", {
    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                    Fitness = c("1",
                    "1 + 0.2 * (n_B > 0)",
                    ".9 + 0.4 * (n_A > 0)"
                    ))
    afd3 <- allFitnessEffects(genotFitness = gffd3,
                            frequencyDependentFitness = TRUE,
                            frequencyType = "abs")

    interventions = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "T > 10",
            WhatHappens   = "N = N * 0.8",
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
                            sampleEvery = 0.001,
                            interventions = interventions)

    # it may happen that, in some simulations, the population collapses, in that case, 
    # pops by time is null, and cannot be checked

    # we can check genotype by genotype that when an intervention ocurs, their population lowers
    # when T = 10:
        #Genotype WT
    if((sfd3$pops.by.time[9997:9997, 2:2] > 0) & (sfd3$pops.by.time[10005:10005, 2:2] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[9997:9997, 2:2]), (sfd3$pops.by.time[10005:10005, 2:2]))
    }
        #Genotype A
    if((sfd3$pops.by.time[9997:9997, 3:3] > 0) & (sfd3$pops.by.time[10005:10005, 3:3] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[9997:9997, 3:3]), (sfd3$pops.by.time[10005:10005, 3:3]))
    }    #Genotype B
    if((sfd3$pops.by.time[9997:9997, 4:4] > 0) & (sfd3$pops.by.time[10005:10005, 4:4] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[9997:9997, 4:4]), (sfd3$pops.by.time[10005:10005, 4:4]))
    }
    # when T = 20
        #Genotype WT
    if((sfd3$pops.by.time[19997:19997, 2:2] > 0) & (sfd3$pops.by.time[20005:20005, 2:2] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[19997:19997, 2:2]), (sfd3$pops.by.time[20005:20005, 2:2]))
    }
        #Genotype A
    if((sfd3$pops.by.time[19997:19997, 3:3] > 0) & (sfd3$pops.by.time[20005:20005, 3:3] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[19997:19997, 3:3]), (sfd3$pops.by.time[20005:20005, 3:3]))
    }
        #Genotype B
    if((sfd3$pops.by.time[19997:19997, 4:4] > 0) & (sfd3$pops.by.time[20005:20005, 4:4] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[19997:19997, 4:4]), (sfd3$pops.by.time[20005:20005, 4:4]))
    }
    # when T = 30
        #Genotype WT
    if((sfd3$pops.by.time[29997:29997, 2:2] > 0) & (sfd3$pops.by.time[30005:30005, 2:2] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[29997:29997, 2:2]), (sfd3$pops.by.time[30005:30005, 2:2]))
    }
        #Genotype A
    if((sfd3$pops.by.time[29997:29997, 3:3] > 0) & (sfd3$pops.by.time[30005:30005, 3:3] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[29997:29997, 3:3]), (sfd3$pops.by.time[30005:30005, 3:3]))
    }
        #Genotype B
    if((sfd3$pops.by.time[29997:29997, 4:4] > 0) & (sfd3$pops.by.time[30005:30005, 4:4] > 0)){
        testthat::expect_gt((sfd3$pops.by.time[29997:29997, 4:4]), (sfd3$pops.by.time[30005:30005, 4:4]))
    }

        # it is expected from the interventions to reduce the population in 3 different times:
    # when T = 10
    tot_pop_t_10_int = (sfd3$pops.by.time[10005:10005, 2:2] + sfd3$pops.by.time[10005:10005, 3:3] 
                          + sfd3$pops.by.time[10005:10005, 4:4])
    tot_pop_t_before10_int = (sfd3$pops.by.time[9997:9997, 2:2] + sfd3$pops.by.time[9997:9997, 3:3] 
                                + sfd3$pops.by.time[9997:9997, 4:4])
    testthat::expect_gt(tot_pop_t_before10_int, tot_pop_t_10_int)
    # when T = 20
    tot_pop_t_20_int = (sfd3$pops.by.time[20005:20005, 2:2] + sfd3$pops.by.time[20005:20005, 3:3] 
                          + sfd3$pops.by.time[20005:20005, 4:4])
    tot_pop_t_before20_int = (sfd3$pops.by.time[19997:19997, 2:2] + sfd3$pops.by.time[19997:19997, 3:3] 
                            + sfd3$pops.by.time[19997:19997, 4:4])
    testthat::expect_gt(tot_pop_t_before20_int, tot_pop_t_20_int)
    # when T = 30
    tot_pop_t_30_int = (sfd3$pops.by.time[30005:30005, 2:2] + sfd3$pops.by.time[30005:30005, 3:3] 
                          + sfd3$pops.by.time[30005:30005, 4:4])
    tot_pop_t_before30_int = (sfd3$pops.by.time[29997:29997, 2:2] + sfd3$pops.by.time[29997:29997, 3:3] 
                            + sfd3$pops.by.time[29997:29997, 4:4])
    testthat::expect_gt(tot_pop_t_before30_int, tot_pop_t_30_int)

    # finally we check that the population is in range within the specified in WhatHappens
    # since the intervention reduces all the population to thee 80% of its original value, we will check
    # that in fact the population reduces (aprox) in that percentage
    testthat::expect_gt(tot_pop_t_before10_int * 0.81, tot_pop_t_10_int)
    testthat::expect_gt(tot_pop_t_before20_int * 0.81, tot_pop_t_20_int)
    testthat::expect_gt(tot_pop_t_before30_int * 0.81, tot_pop_t_30_int)
})

# test 9 and 10 found in test.Z-intervention.R

test_that("11. Intervening over 4 genotypes both over specific genotype and total population (McFL) | Trigger depends on N",{
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
            WhatHappens   = "n_A = n_A * 0.5",
            Repetitions   = 0,   
            Periodicity   = Inf
        ),
        list(
            ID            = "intOverB",
            Trigger       = "T >= 2.2",
            WhatHappens   = "n_B = n_B * 0.5",
            Repetitions   = 0,   
            Periodicity   = Inf
        ),
        list(
            ID            = "intOverC",
            Trigger       = "T >= 3.2",
            WhatHappens   = "n_C = n_C * 0.5",
            Repetitions   = 0,   
            Periodicity   = Inf
        ),
        list(
            ID            = "intOverD",
            Trigger       = "T >= 4.2",
            WhatHappens   = "n_D = n_D * 0.5",
            Repetitions   = 0,   
            Periodicity   = Inf
        )
    )

    interventions = createInterventions(interventions, afd3)

    sfd3_with_ints <- oncoSimulIndiv( afd3,
                            model = "McFLD",
                            onlyCancer = FALSE,
                            finalTime = 6,
                            mu = 1e-4,
                            sampleEvery = 0.01,
                            keepPhylog = FALSE,
                            seed = NULL,
                            errorHitMaxTries = FALSE,
                            errorHitWallTime = FALSE,
                            initMutant = c("A", "B", "C", "D", "E"),
                            initSize = c(20000, 20000, 30, 10, 200),
                            interventions = interventions)

    
    # when T = 1.2 population of genotype A is intervened
    if((sfd3_with_ints$pops.by.time[118:118, 2:2] > 0) & (sfd3_with_ints$pops.by.time[123:123, 2:2] > 0)){
        testthat::expect_gt((sfd3_with_ints$pops.by.time[118:118, 2:2]), (sfd3_with_ints$pops.by.time[125:125, 2:2]))
    }

    # when T = 2.2 population of genotype B is intervened
    if((sfd3_with_ints$pops.by.time[218:218, 3:3] > 0) & (sfd3_with_ints$pops.by.time[225:225, 3:3] > 0)){
        testthat::expect_gt((sfd3_with_ints$pops.by.time[218:218, 3:3]), (sfd3_with_ints$pops.by.time[225:225, 3:3]))
    }
    # when T = 3.2 population of genotype C is intervened
    if((sfd3_with_ints$pops.by.time[318:318, 4:4] > 0) & (sfd3_with_ints$pops.by.time[323:323, 4:4] > 0)){
        testthat::expect_gt((sfd3_with_ints$pops.by.time[318:318, 4:4]), (sfd3_with_ints$pops.by.time[323:323, 4:4]))
    }
    # when T = 4.2 population of genotype D is intervened
    if((sfd3_with_ints$pops.by.time[418:418, 5:5] > 0) & (sfd3_with_ints$pops.by.time[425:425, 5:5] > 0)){
        testthat::expect_gt((sfd3_with_ints$pops.by.time[418:418, 5:5]), (sfd3_with_ints$pops.by.time[425:425, 5:5]))
    }
})



test_that("12. Intervening over 4 genotypes both over specific and total population (Exp) | Trigger depends on T", {

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
            WhatHappens   = "n_A = n_A * 0.5",
            Repetitions   = 0,   
            Periodicity   = Inf
        ),
        list(
            ID            = "intOverB",
            Trigger       = "T >= 2.2",
            WhatHappens   = "n_B = n_B * 0.5",
            Repetitions   = 0,   
            Periodicity   = Inf
        ),
        list(
            ID            = "intOverC",
            Trigger       = "T >= 3.2",
            WhatHappens   = "n_C = n_C * 0.5",
            Repetitions   = 0,   
            Periodicity   = Inf
        ),
        list(
            ID            = "intOverD",
            Trigger       = "T >= 4.2",
            WhatHappens   = "n_D = n_D * 0.5",
            Repetitions   = 0,   
            Periodicity   = Inf
        )
    )

    interventions = createInterventions(interventions, afd3)

    sfd3_with_ints <- oncoSimulIndiv( afd3,
                            model = "Exp",
                            onlyCancer = FALSE,
                            finalTime = 6,
                            mu = 1e-4,
                            sampleEvery = 0.01,
                            keepPhylog = FALSE,
                            seed = NULL,
                            errorHitMaxTries = FALSE,
                            errorHitWallTime = FALSE,
                            initMutant = c("A", "B", "C", "D", "E"),
                            initSize = c(20000, 20000, 30, 10, 200),
                            interventions = interventions)

    #plot(sfd3, show = "genotypes", type = "line")

    # when T = 1.2 population of genotype A is intervened
    if((sfd3_with_ints$pops.by.time[118:118, 2:2] > 0) & (sfd3_with_ints$pops.by.time[123:123, 2:2] > 0)){
        testthat::expect_gt((sfd3_with_ints$pops.by.time[118:118, 2:2]), (sfd3_with_ints$pops.by.time[125:125, 2:2]))
    }

    # when T = 2.2 population of genotype B is intervened
    if((sfd3_with_ints$pops.by.time[218:218, 3:3] > 0) & (sfd3_with_ints$pops.by.time[225:225, 3:3] > 0)){
        testthat::expect_gt((sfd3_with_ints$pops.by.time[218:218, 3:3]), (sfd3_with_ints$pops.by.time[225:225, 3:3]))
    }
    # when T = 3.2 population of genotype C is intervened
    if((sfd3_with_ints$pops.by.time[318:318, 4:4] > 0) & (sfd3_with_ints$pops.by.time[323:323, 4:4] > 0)){
        testthat::expect_gt((sfd3_with_ints$pops.by.time[318:318, 4:4]), (sfd3_with_ints$pops.by.time[323:323, 4:4]))
    }
    # when T = 4.2 population of genotype D is intervened
    if((sfd3_with_ints$pops.by.time[418:418, 5:5] > 0) & (sfd3_with_ints$pops.by.time[425:425, 5:5] > 0)){
        testthat::expect_gt((sfd3_with_ints$pops.by.time[418:418, 5:5]), (sfd3_with_ints$pops.by.time[425:425, 5:5]))
    }
})




test_that("13. Intervening in the Rock-Paper-Scissors model for bacterial community by Kerr.", {
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

    lista_intervenciones = list(
        list(
            ID = "Bothering R strain, by reducing C",
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

    # by reducing C, R wont spread in the population. This will mean that, with the apropiate 
    # periodicty in the intervention, C will never surpass WT. We try to check this in these tests.
    i <- 1
    while(i <= nrow(resultscrs1$pops.by.time)){
        cur_nWT = resultscrs1$pops.by.time[i, 2:2]
        cur_nR = resultscrs1$pops.by.time[i, 4:4]
        testthat::expect_gt(cur_nWT, cur_nR)
        i <- i + 1
    }
})

cat(paste("\n Ending interventions tests", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)