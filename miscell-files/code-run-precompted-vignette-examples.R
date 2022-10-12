library(OncoSimulR)


## intex2 and  simulationwithinterventionsintex2
local({

    gffd3 <- data.frame(Genotype = c("WT", "A", "B"),
                        Fitness = c("1",
                                    "1 + 0.25 * (n_B > 0)",
                                    ".9 + 0.4 * (n_A > 0)"
                                    ))
    afd3 <- allFitnessEffects(genotFitness = gffd3,
                              frequencyDependentFitness = TRUE,
                              frequencyType = "abs")

    osi <- oncoSimulIndiv( afd3,
                          model = "McFLD",
                          onlyCancer = FALSE,
                          finalTime = 200,
                          mu = 1e-4,
                          initSize = 5000,
                          sampleEvery = 0.001,
                          keepEvery = 1)
    osi$other$userVarValues <- NULL
    save(file = "osi_intex2.RData", osi)

    
    intervention_tot_pop = list(
        list(
            ID            = "intOverTotPop",
            Trigger       = "T > 40",
            WhatHappens   = "N = N * 0.2",
            Repetitions   = 2,   
            Periodicity   = 20
        )
    )

    intervention_tot_pop <-
        createInterventions(intervention_tot_pop, afd3)

    osi_with_ints <- oncoSimulIndiv( afd3,
                                    model = "McFLD",
                                    onlyCancer = FALSE,
                                    finalTime = 200,
                                    mu = 1e-4,
                                    initSize = 5000,
                                    sampleEvery = 0.001,
                                    interventions = intervention_tot_pop,
                                    keepEvery = 1)
    
    osi_with_ints$other$userVarValues <- NULL
    save(file = "osi_with_ints.RData", osi_with_ints)
})

## HansenExample3

local({
    dfat4 <- data.frame(Genotype = c("WT", "A", "B"), 
                        Fitness = c("n_/n_",
                                    "1.005",
                                    "1.1"
                                    ))
    afat4 <- allFitnessEffects(genotFitness = dfat4, 
                               frequencyDependentFitness = TRUE)
    interventions <- list(
        list(ID           = "i1",
             Trigger       = "T > 10",
             WhatHappens   = "n_B = n_B*0.8",
             Periodicity   = 1,
             Repetitions   = Inf
             )
    )

    interventions <- createInterventions(interventions, afat4) 
    set.seed(1) ## for reproducibility
    atex4 <- oncoSimulIndiv(afat4,
                            model = "McFLD", 
                            onlyCancer = FALSE, 
                            finalTime = 2000,
                            mu = 1e-4,
                            initSize = c(10000, 50, 1000), 
                            initMutant = c("WT", "A", "B"),
                            keepPhylog = FALSE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE,
                            interventions = interventions,
                            keepEvery = 1)

    atex4$other$userVarValues <- NULL
    save(file = "HansenExample3.RData", atex4)

})



## HansenExampleAT4
local({
    dfat5 <- data.frame(Genotype = c("WT", "A", "B"), 
                        Fitness = c("n_/n_",
                                    "1.005",
                                    "1.1"
                                    ))
    afat5 <- allFitnessEffects(genotFitness = dfat5, 
                               frequencyDependentFitness = TRUE)

    userVars <- list(
        list(Name           = "measure",
             Value       = 0
             ),
        list(Name           = "lastTime",
             Value       = 0
             ),
        list(Name           = "treatment",
             Value       = 0
             ),
        list(Name           = "totalPopMeasured",
             Value       = 0
             )
    )

    userVars <- createUserVars(userVars)

    rules <- list(
        list(ID = "rule_1",
             Condition = "T - lastTime < 10",
             Action = "measure = 0"
             ),
        list(ID = "rule_2",
             Condition = "T - lastTime >= 10",
             Action = "measure = 1;lastTime = T"
             ),
        list(ID = "rule_3",
             Condition = "measure == 1",
             Action = "totalPopMeasured = n_A + n_B"
             ),
        list(ID = "rule_4",
             Condition = "totalPopMeasured < 2000",
             Action = "treatment = 0"
             ),
        list(ID = "rule_5",
             Condition = "totalPopMeasured >= 2000",
             Action = "treatment = 1"
             )
    )

    rules <- createRules(rules, afat5)

    interventions <- list(
        list(ID           = "i1",
             Trigger       = "treatment == 1",
             WhatHappens   = "n_B = n_B*0.8",
             Periodicity   = 1,
             Repetitions   = Inf
             )
    )

    interventions <- createInterventions(interventions, afat5)

    set.seed(1) ## for reproducibility
    atex5 <- oncoSimulIndiv(afat5,
                            model = "McFLD", 
                            onlyCancer = FALSE, 
                            finalTime = 1500,
                            mu = 1e-4,
                            initSize = c(10000, 50, 1000), 
                            initMutant = c("WT", "A", "B"),
                            keepPhylog = FALSE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE,
                            userVars = userVars,
                            rules = rules,
                            interventions = interventions,
                            keepEvery = 1)
    ## Not needed for the plot
    atex5$other$userVarValues <- NULL
    save(file = "HansenExampleAT4.RData", atex5)
})





## AdaptiveTherapyExample3

local({
    dfat2 <- data.frame(Genotype = c("WT", "A", "B"), 
                        Fitness = c("1",
                                    "1 + 0.2 * (n_B > 10)",
                                    ".9 + 0.4 * (n_A > 10)"
                                    ))
    afat2 <- allFitnessEffects(genotFitness = dfat2, 
                               frequencyDependentFitness = TRUE)
    

    userVars <- list(
        list(Name           = "genADiff",
             Value       = 0
             ),
        list(Name           = "genBDiff",
             Value       = 0
             ),
        list(Name           = "genWTDiff",
             Value       = 0
             )
    )

    userVars <- createUserVars(userVars)

    rules <- list(
        list(ID = "rule_1",
             Condition = "TRUE",
             Action = "genADiff = b_1-d_1"
             ),
        list(ID = "rule_2",
             Condition = "TRUE",
             Action = "genBDiff = b_2-d_2"
             ),
        list(ID = "rule_3",
             Condition = "TRUE",
             Action = "genWTDiff = b_-d_"
             ),
        list(ID = "rule_4",
             Condition = "n_A < 1000",
             Action = "genADiff = -1"
             ),
        list(ID = "rule_5",
             Condition = "n_B < 1000",
             Action = "genBDiff = -1"
             ),
        list(ID = "rule_6",
             Condition = "n_ < 1000",
             Action = "genWTDiff = -1"
             )
    )

    rules <- createRules(rules, afat2)

    interventions <- list(
        list(ID           = "i1",
             Trigger       = "genADiff > genBDiff and genADiff > genWTDiff",
             WhatHappens   = "n_A = n_A*0.5",
             Periodicity   = 10,
             Repetitions   = Inf
             ),
        list(ID           = "i2",
             Trigger       = "genBDiff > genADiff and genBDiff > genWTDiff",
             WhatHappens   = "n_B = n_B*0.5",
             Periodicity   = 10,
             Repetitions   = Inf
             ),
        list(ID           = "i3",
             Trigger       = "genWTDiff > genADiff and genWTDiff > genBDiff",
             WhatHappens   = "n_ = n_*0.5",
             Periodicity   = 10,
             Repetitions   = Inf
             )
    )

    interventions <- createInterventions(interventions, afat2)
    set.seed(1) ## for reproducibility
    atex2 <- oncoSimulIndiv(afat2,
                            model = "McFLD", 
                            onlyCancer = FALSE, 
                            finalTime = 200,
                            mu = 1e-4,
                            initSize = 5000, 
                            keepPhylog = FALSE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE,
                            userVars = userVars,
                            rules = rules,
                            interventions = interventions,
                            keepEvery = 1)

    atex2$other$userVarValues <- NULL
    save(file = "AdaptiveTherapyExample3.RData", atex2)
})



## AdaptiveTherapyComplexExample

local({
    dfat3 <- data.frame(Genotype = c("WT", "A", "B"), 
                        Fitness = c("1",
                                    "0.8 + 0.2 * (n_B > 10) + 0.1 (n_A > 10)",
                                    "0.8 + 0.25 * (n_B > 10)"
                                    ))
    afat3 <- allFitnessEffects(genotFitness = dfat3, 
                               frequencyDependentFitness = TRUE)
    

    userVars <- list(
        list(Name           = "lastMeasuredA",
             Value       = 0
             ),
        list(Name           = "lastMeasuredB",
             Value       = 0
             ),
        list(Name           = "previousA",
             Value       = 0
             ),
        list(Name           = "previousB",
             Value       = 0
             ),
        list(Name           = "lastTime",
             Value       = 0
             ),
        list(Name           = "measure",
             Value       = 0
             ),
        list(Name           = "treatment",
             Value       = 0
             )
    )

    userVars <- createUserVars(userVars)

    rules <- list(
        list(ID = "rule_1",
             Condition = "T - lastTime < 10",
             Action = "measure = 0"
             ),
        list(ID = "rule_2",
             Condition = "T - lastTime >= 10",
             Action = "measure = 1;lastTime = T"
             ),
        list(ID = "rule_3",
             Condition = "measure == 1",
             Action = "previousA = lastMeasuredA;previousB = lastMeasuredB;lastMeasuredA = n_A;lastMeasuredB = n_B"
             ),
        list(ID = "rule_4",
             Condition = "TRUE",
             Action = "treatment = 0"
             ),
        list(ID = "rule_5",
             Condition = "lastMeasuredA + lastMeasuredB > 100",
             Action = "treatment = 1"
             ),
        list(ID = "rule_6",
             Condition = "lastMeasuredA - PreviousA > 500",
             Action = "treatment = 2"
             ),
        list(ID = "rule_7",
             Condition = "lastMeasuredB - PreviousB > 500",
             Action = "treatment = 3"
             ),
        list(ID = "rule_8",
             Condition = "lastMeasuredA - PreviousA > 500 and lastMeasuredB - PreviousB > 500",
             Action = "treatment = 4"
             )
    )

    rules <- createRules(rules, afat3)


    interventions <- list(
        list(ID           = "basicTreatment",
             Trigger       = "treatment == 1",
             WhatHappens   = "N = 0.8*N",
             Periodicity   = 10,
             Repetitions   = Inf
             ),
        list(ID           = "treatmentOverA",
             Trigger       = "treatment == 2 or treatment == 4",
             WhatHappens   = "n_B = n_B*0.3",
             Periodicity   = 20,
             Repetitions   = Inf
             ),
        list(ID           = "treatmentOverB",
             Trigger       = "treatment == 3 or treatment == 4",
             WhatHappens   = "n_B = n_B*0.3",
             Periodicity   = 20,
             Repetitions   = Inf
             ),
        list(ID           = "intervention",
             Trigger       = "lastMeasuredA+lastMeasuredB > 5000",
             WhatHappens   = "N = 0.1*N",
             Periodicity   = 70,
             Repetitions   = Inf
             )
    )

    interventions <- createInterventions(interventions, afat3)

    set.seed(1) ## for reproducibility
    atex2 <- oncoSimulIndiv(afat3,
                            model = "McFLD", 
                            onlyCancer = FALSE, 
                            finalTime = 200,
                            mu = 1e-4,
                            initSize = 5000, 
                            keepPhylog = FALSE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE,
                            userVars = userVars,
                            rules = rules,
                            interventions = interventions,
                            keepEvery = 1)

    atex2$other$userVarValues <- NULL
    save(file = "AdaptiveTherapyComplexExample.RData", atex2)
})

## ## userVarsBasicExample
## local({
##     dfuv2 <- data.frame(Genotype = c("WT", "B", "A", "B, A", "C, A"),
##                         Fitness = c("0*n_",
##                                     "1.5",
##                                     "1.002",
##                                     "1.003",
##                                     "1.004"))

##     afuv2 <- allFitnessEffects(genotFitness = dfuv2,
##                                frequencyDependentFitness = TRUE,
##                                frequencyType = "abs")

##     userVars <- list(
##         list(Name           = "genAProp",
##              Value       = 0.5
##              ),
##         list(Name           = "genBProp",
##              Value       = 0.5
##              ),
##         list(Name           = "genABProp",
##              Value       = 0.0
##              ),
##         list(Name           = "genACProp",
##              Value       = 0.0
##              )
##     )

##     userVars <- createUserVars(userVars)

##     rules <- list(
##         list(ID = "rule_1",
##              Condition = "TRUE",
##              Action = "genBProp = n_B/N"
##              ),
##         list(ID = "rule_2",
##              Condition = "TRUE",
##              Action = "genAProp = n_A/N"
##              ),
##         list(ID = "rule_3",
##              Condition = "TRUE",
##              Action = "genABProp = n_A_B/N"
##              ),
##         list(ID = "rule_4",
##              Condition = "TRUE",
##              Action = "genACProp = n_A_C/N"
##              )
##     )

##     rules <- createRules(rules, afuv2)

##     set.seed(1)
##     uvex2 <- oncoSimulIndiv(
##         afuv2, 
##         model = "McFLD",
##         mu = 1e-4,
##         sampleEvery = 0.01,
##         initSize = c(20000, 20000),
##         initMutant = c("A", "B"),
##         finalTime = 10,
##         onlyCancer = FALSE,
##         userVars = userVars,
##         rules = rules,
##         keepEvery = 0.1
##     )
##     uvex2$other$userVarValues <- NULL
##     save(file = "usersVarsBasicExample.RData", uvex2)

## })


## userVarsBasicExample2
local({
    dfuv3 <- data.frame(Genotype = c("WT", "A", "B"), 
                        Fitness = c("1",
                                    "1 + 0.2 * (n_B > 10)",
                                    ".9 + 0.4 * (n_A > 10)"
                                    ))
    afuv3 <- allFitnessEffects(genotFitness = dfuv3, 
                               frequencyDependentFitness = TRUE)

    userVars <- list(
        list(Name           = "genWTRateDiff",
             Value       = 0.5
             ),
        list(Name           = "genARateDiff",
             Value       = 0.5
             ),
        list(Name           = "genBRateDiff",
             Value       = 0.0
             )
    )

    userVars <- createUserVars(userVars)

    rules <- list(
        list(ID = "rule_1",
             Condition = "TRUE",
             Action = "genWTRateDiff = b_-d_"
             ),
        list(ID = "rule_2",
             Condition = "TRUE",
             Action = "genARateDiff = b_1-d_1"
             ),
        list(ID = "rule_3",
             Condition = "TRUE",
             Action = "genBRateDiff = b_2-d_2"
             )
    )

    rules <- createRules(rules, afuv3)

    set.seed(1)

    uvex3 <- oncoSimulIndiv(afuv3,
                            model = "McFLD", 
                            onlyCancer = FALSE, 
                            finalTime = 150,
                            mu = 1e-4,
                            initSize = 5000, 
                            keepPhylog = FALSE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE,
                            userVars = userVars,
                            rules = rules,
                            keepEvery = 1)

    uvex3$other$userVarValues <- NULL
    save(file = "usersVarsBasicExample2.RData", uvex3)

})



## userVarsOncoSimulIndivExample
local({
    dfuv <- data.frame(Genotype = c("WT", "A", "B"),
                       Fitness = c("1",
                                   "1 + 0.2 * (n_B > 0)",
                                   ".9 + 0.4 * (n_A > 0)"
                                   ))
    afuv <- allFitnessEffects(genotFitness = dfuv,
                              frequencyDependentFitness = TRUE,
                              frequencyType = "abs")

    userVars <- createUserVars(userVars)
    rules <- createRules(rules, afuv)

    uvex <- oncoSimulIndiv(
        afuv, 
        model = "McFLD",
        mu = 1e-4,
        sampleEvery = 0.001,
        initSize = c(20000, 20000),
        initMutant = c("A", "B"),
        finalTime = 5.2,
        onlyCancer = FALSE,
        userVars = userVars,
        rules = rules
    )

    uvex$other$userVarValues <- NULL
    save(file = "usersVarsOncoSimulIndivExample.RData", uvex)
})


## interventionsOncoSimulIndivExample

local({
    fa1 <- data.frame(Genotype = c("WT", "A", "B"),
                      Fitness = c("n_*0",
                                  "1.5",
                                  "1"))

    afd3 <- allFitnessEffects(genotFitness = fa1,
                              frequencyDependentFitness = TRUE,
                              frequencyType = "abs")

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

    ep2$other$userVarValues <- NULL
    save(file = "interventionsOncoSimulIndivExample.RData", ep2)
})
