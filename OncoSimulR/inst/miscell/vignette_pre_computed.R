### Precompute some examples for the vignette


rm(list = ls())
library(OncoSimulR)

local({
    f_cells <- function(c1, c2, c3, r11, r12, r13, 
                        r21, r22, r23, r31, r32, r33, M, awt = 1e-4,
                        gt = c("WT", "OC", "OB", "MM")) {
        data.frame(Genotype = gt,
                   Fitness = c(
                       paste0("max(0.1, 1 - ", awt, " * (f_OC + f_OB+f_MM)*N)"),
                       paste0("1", "+(((f_OC * (", M, "-1)+1)*", c1, ")/", M, ")*",r11,
                              "+((((1-f_MM) * (", M, "-1)-f_OC*(", M, "-1)-1)*", c2, ")/", M, ")*", r12,
                              "+(((", M, "-(1-f_MM)*(", M, "-1))*", c3, ")/", M, ")*", r13,
                              "-", c1
                              ),
                       paste0("1", "+(((f_OB*(", M, "-1)+1)*", c2, ")/", M, ")*", r22,
                              "+((((1-f_OC)*(", M, "-1)-f_OB*(", M, "-1)-1)*", c3, ")/", M, ")*", r23,
                              "+(((", M, "-(1-f_OC)*(", M, "-1))*", c1, ")/", M, ")*", r21,
                              "-", c2
                              ),
                       paste0("1", "+(((f_MM*(", M, "-1)+1)*", c3, ")/", M, ")*", r33,
                              "+((((1-f_OB)*(", M, "-1)-f_MM*(", M, "-1)-1)*", c1, ")/", M, ")*", r31,
                              "+(((", M, "-(1-f_OB)*(", M, "-1))*", c2, ")/", M, ")*", r32,
                              "-", c3
                              )
                   )
                  ,stringsAsFactors = FALSE
                   )
    }


    
    N <- 40000
    M <- 10
    c1 <- 1
    c2 <- 1.2
    c3 <- 1.4
    r11 <- 0
    r12 <- 1
    r13 <- 2.5
    r21 <- 1
    r22 <- 0
    r23 <- -0.3
    r31 <- 2.5
    r32 <- 0
    r33 <- 0

    fe_cells <-
        allFitnessEffects(
            genotFitness =
                f_cells(c1, c2, c3, r11, r12, r13, 
                        r21, r22, r23, r31, r32, r33, M,
                        gt = c("WT", "OC", "OB", "MM")),
            frequencyDependentFitness = TRUE,
            frequencyType = "rel")

    ## Simulated trajectories

    data(smyelo3v57)

    set.seed(2)
    smyelo3v57 <- oncoSimulIndiv(fe_cells,
                                 model = "McFL",
                                 onlyCancer = FALSE,
                                 finalTime = 20, 
                                 mu = c("OC"=1e-1, "OB"=1e-1, "MM"=1e-4),
                                 initSize = N,
                                 keepPhylog = FALSE,
                                 seed = NULL,
                                 errorHitMaxTries = FALSE,
                                 errorHitWallTime = FALSE,
                                 keepEvery = 0.1)

    smyelo3v57$other$userVarValues <- NULL
    smyelo3v57$PerSampleStats <- NULL
    smyelo3v57$interventionTimes <- NULL
    save(file = "../../data/smyelo3v57.RData", smyelo3v57)  
})


local({
    create_fe <- function(a, b, c, d, e, f, g,
                          gt = c("WT", "A", "P", "C")) {
        data.frame(Genotype = gt,
                   Fitness = c(
                       paste0("1 + ",
                              d, " * f_A ",
                              "- ", c, " * f_C"),
                       paste0("1 - ", a, 
                              " + ", d, " + ",
                              f, " * f_A ",
                              "- ", c, " * f_C"),
                       paste0("1 + ", g, " + ",
                              d, " * f_A ",
                              "- ", c, " * (1 + ",
                              g, ") * f_C"),
                       paste0("1 - ", b, " + ",
                              e, " * f_ + ",
                              "(", d, " + ", 
                              e, ") * f_A + ",
                              e , " * f_P")))
    }

    afe_3_b <- allFitnessEffects(
        genotFitness =
            create_fe(0.02, 0.04, 0.08, 0.1,
                      0.15, 0.1, 0.05),
        frequencyDependentFitness = TRUE)
    set.seed(2)
    ## Use a short finalTime, for speed of vignette execution
    s_3_b <- oncoSimulIndiv(afe_3_b,
                            model = "McFL", 
                            onlyCancer = FALSE, 
                            finalTime = 100, 
                            mu = 1e-4,
                            initSize = 5000, 
                            keepPhylog = FALSE, 
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE,
                            keepEvery = 1)

    s_3_b$other$userVarValues <- NULL
    s_3_b$PerSampleStats <- NULL
    s_3_b$interventionTimes <- NULL
    save(file = "../../data/s_3_b.RData", s_3_b)  
})


## userVarsBasicExample
local({
    dfuv2 <- data.frame(Genotype = c("WT", "B", "A", "B, A", "C, A"),
                        Fitness = c("0*n_",
                                    "1.5",
                                    "1.002",
                                    "1.003",
                                    "1.004"))

    afuv2 <- allFitnessEffects(genotFitness = dfuv2,
                               frequencyDependentFitness = TRUE,
                               frequencyType = "abs")

    userVars <- list(
        list(Name           = "genAProp",
             Value       = 0.5
             ),
        list(Name           = "genBProp",
             Value       = 0.5
             ),
        list(Name           = "genABProp",
             Value       = 0.0
             ),
        list(Name           = "genACProp",
             Value       = 0.0
             )
    )

    userVars <- createUserVars(userVars)

    rules <- list(
        list(ID = "rule_1",
             Condition = "TRUE",
             Action = "genBProp = n_B/N"
             ),
        list(ID = "rule_2",
             Condition = "TRUE",
             Action = "genAProp = n_A/N"
             ),
        list(ID = "rule_3",
             Condition = "TRUE",
             Action = "genABProp = n_A_B/N"
             ),
        list(ID = "rule_4",
             Condition = "TRUE",
             Action = "genACProp = n_A_C/N"
             )
    )

    rules <- createRules(rules, afuv2)

    set.seed(1)
    uvex2 <- oncoSimulIndiv(
        afuv2, 
        model = "McFLD",
        mu = 1e-4,
        sampleEvery = 0.01,
        initSize = c(20000, 20000),
        initMutant = c("A", "B"),
        finalTime = 10,
        onlyCancer = FALSE,
        userVars = userVars,
        rules = rules,
        keepEvery = 0.1
    )
    uvex2$PerSampleStats <- NULL
    uvex2$interventionTimes <- NULL
    ## uvex2$other$userVarValues <- NULL
    save(file = "../../data/uvex2.RData", uvex2)

})



local({
    dfT2 <- data.frame(Genotype = c("WT", "A", "B"),
                       Fitness = c(
                           "1",
                           "if (T>0 and T<50) 0; else if (T>100 and T<150) 0.05; else 1.2 + 0.35*f_;",
                           "0.8 + 0.45*(f_)"
                       ),
                       stringsAsFactors = FALSE)
    
    afeT2 <- allFitnessEffects(genotFitness = dfT2,
                               frequencyDependentFitness = TRUE,
                               frequencyType = "rel")
    
    set.seed(1)
    simT2 <- oncoSimulIndiv(afeT2,
                            model = "McFL",
                            mu = 1e-5,
                            initSize = 10000,
                            finalTime = 225,
                            onlyCancer = FALSE,
                            keepEvery = 1,
                            seed = NULL)

    simT2$other$userVarValues <- NULL
    simT2$PerSampleStats <- NULL
    simT2$interventionTimes <- NULL
    save(file = "../../data/simT2.RData", simT2)
})


local({

    ## Healthy     Sensitive   Resistant
    a <- 1;       b <- 0.5;   c <- 0.5    # Healthy
    d <- 1;       e <- 1.25;  f <- 0.7    # Sensitive
    g <- 0.975;   h <- -0.5;  i <- 0.75   # Resistant

    drug_eff <- 0.01

    wt_fitness <- paste0(a, "*f_+", b, "*f_S+", c, "*f_S_R")
    sens_fitness <- paste0(d, "*f_+", e, "*f_S+", f, "*f_S_R")
    res_fitness <- paste0(g, "*f_+", h, "*f_S+", i, "*f_S_R")

    players_1 <- data.frame(Genotype = c("WT", "S", "R", "S, R"),
                            Fitness = c(wt_fitness,                                               #WT
                                        paste0("if (T>50) ", drug_eff, "*(",sens_fitness, ")",";  
                                           else ", sens_fitness, ";"),                        #S
                                        "0",                                                      #R
                                        res_fitness),                                             #S,R
                            stringsAsFactors = FALSE)

    period_1 <- allFitnessEffects(genotFitness = players_1,
                                  frequencyDependentFitness = TRUE,
                                  frequencyType = "rel")

    set.seed(2)

    final_time <- 170 ## for speed
    simul_period_1 <- oncoSimulIndiv(period_1,
                                     model = "McFL",
                                     onlyCancer = FALSE,
                                     finalTime = final_time,
                                     mu = 0.01,
                                     initSize = 5000,
                                     keepPhylog = FALSE,
                                     seed = NULL,
                                     keepEvery = .2)

    simul_period_1$other$userVarValues <- NULL
    simul_period_1$PerSampleStats <- NULL
    simul_period_1$interventionTimes <- NULL
    save(file = "../../data/simul_period_1.RData", simul_period_1)

})



local({
    dfT3 <- data.frame(Genotype = c("WT", "A", "B"),
                       Fitness = c(
                           "1",
                           "1 + 0.2 * (n_2 > 10)",
                           "if (T>50 and T<80) 0.80; else 0.9 + 0.4 * (n_1 > 10)"),
                       stringsAsFactors = FALSE)

    afeT3 <- allFitnessEffects(genotFitness = dfT3,
                               frequencyDependentFitness = TRUE,
                               frequencyType = "abs")
    set.seed(2)
    simT3 <- oncoSimulIndiv(afeT3,
                            model = "McFLD",
                            mu = 1e-4,
                            initSize = 5000,
                            finalTime = 500,
                            onlyCancer = FALSE,
                            seed = NULL,
                            errorHitWallTime = FALSE,
                            errorHitMaxTries = FALSE,
                            keepEvery = 1)
    
    simT3$other$userVarValues <- NULL
    simT3$PerSampleStats <- NULL
    simT3$interventionTimes <- NULL
    save(file = "../../data/simT3.RData", simT3)
})

local({
    create_fe <- function(a, b, c, d, e, f, g,
                          gt = c("WT", "A", "P", "C")) {
        data.frame(Genotype = gt,
                   Fitness = c(
                       paste0("1 + ",
                              d, " * f_A ",
                              "- ", c, " * f_C"),
                       paste0("1 - ", a, 
                              " + ", d, " + ",
                              f, " * f_A ",
                              "- ", c, " * f_C"),
                       paste0("1 + ", g, " + ",
                              d, " * f_A ",
                              "- ", c, " * (1 + ",
                              g, ") * f_C"),
                       paste0("1 - ", b, " + ",
                              e, " * f_ + ",
                              "(", d, " + ", 
                              e, ") * f_A + ",
                              e , " * f_P")))
    }
    afe_3_a <- allFitnessEffects(
        genotFitness =
            create_fe(0.02, 0.04, 0.08, 0.06,
                      0.15, 0.1, 0.06),
        frequencyDependentFitness = TRUE)
    s_3_a <- oncoSimulIndiv(afe_3_a,
                            model = "McFL", 
                            onlyCancer = FALSE, 
                            finalTime = 160,
                            mu = 1e-4,
                            initSize = 5000, 
                            keepPhylog = FALSE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE,
                            keepEvery = 1)
    s_3_a$other$userVarValues <- NULL
    s_3_a$PerSampleStats <- NULL
    s_3_a$interventionTimes <- NULL
    save(file = "../../data/s_3_a.RData", s_3_a)
})


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
    osi$PerSampleStats <- NULL
    osi$interventionTimes <- NULL
    save(file = "../../data/osi.RData", osi)

    
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
    ## osi_with_ints$PerSampleStats <- NULL
    osi_with_ints$other$interventionTimes <- NULL
    save(file = "../../data/osi_with_ints.RData",
         osi_with_ints)
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
                            keepEvery = 2)

    atex4$other$userVarValues <- NULL
    atex4$PerSampleStats <- NULL
    atex4$other$interventionTimes <- NULL
    save(file = "../../data/atex4.RData", atex4)

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
                            keepEvery = 2)
    ## Not needed for the plot
    atex5$other$userVarValues <- NULL
    atex5$PerSampleStats <- NULL
    atex5$other$interventionTimes <- NULL
    save(file = "../../data/atex5.RData", atex5)
})





## ## AdaptiveTherapyExample3

## local({
##     dfat2 <- data.frame(Genotype = c("WT", "A", "B"), 
##                         Fitness = c("1",
##                                     "1 + 0.2 * (n_B > 10)",
##                                     ".9 + 0.4 * (n_A > 10)"
##                                     ))
##     afat2 <- allFitnessEffects(genotFitness = dfat2, 
##                                frequencyDependentFitness = TRUE)


##     userVars <- list(
##         list(Name           = "genADiff",
##              Value       = 0
##              ),
##         list(Name           = "genBDiff",
##              Value       = 0
##              ),
##         list(Name           = "genWTDiff",
##              Value       = 0
##              )
##     )

##     userVars <- createUserVars(userVars)

##     rules <- list(
##         list(ID = "rule_1",
##              Condition = "TRUE",
##              Action = "genADiff = b_1-d_1"
##              ),
##         list(ID = "rule_2",
##              Condition = "TRUE",
##              Action = "genBDiff = b_2-d_2"
##              ),
##         list(ID = "rule_3",
##              Condition = "TRUE",
##              Action = "genWTDiff = b_-d_"
##              ),
##         list(ID = "rule_4",
##              Condition = "n_A < 1000",
##              Action = "genADiff = -1"
##              ),
##         list(ID = "rule_5",
##              Condition = "n_B < 1000",
##              Action = "genBDiff = -1"
##              ),
##         list(ID = "rule_6",
##              Condition = "n_ < 1000",
##              Action = "genWTDiff = -1"
##              )
##     )

##     rules <- createRules(rules, afat2)

##     interventions <- list(
##         list(ID           = "i1",
##              Trigger       = "genADiff > genBDiff and genADiff > genWTDiff",
##              WhatHappens   = "n_A = n_A*0.5",
##              Periodicity   = 10,
##              Repetitions   = Inf
##              ),
##         list(ID           = "i2",
##              Trigger       = "genBDiff > genADiff and genBDiff > genWTDiff",
##              WhatHappens   = "n_B = n_B*0.5",
##              Periodicity   = 10,
##              Repetitions   = Inf
##              ),
##         list(ID           = "i3",
##              Trigger       = "genWTDiff > genADiff and genWTDiff > genBDiff",
##              WhatHappens   = "n_ = n_*0.5",
##              Periodicity   = 10,
##              Repetitions   = Inf
##              )
##     )

##     interventions <- createInterventions(interventions, afat2)
##     set.seed(1) ## for reproducibility
##     atex2 <- oncoSimulIndiv(afat2,
##                             model = "McFLD", 
##                             onlyCancer = FALSE, 
##                             finalTime = 200,
##                             mu = 1e-4,
##                             initSize = 5000, 
##                             keepPhylog = FALSE,
##                             seed = NULL, 
##                             errorHitMaxTries = FALSE, 
##                             errorHitWallTime = FALSE,
##                             userVars = userVars,
##                             rules = rules,
##                             interventions = interventions,
##                             keepEvery = 1)

##     atex2$PerSampleStats <- NULL
##     atex2$other$interventionTimes <- NULL
##     save(file = "../../data/AdaptiveTherapyExample3.RData", atex2)
## })



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
    atex2b <- oncoSimulIndiv(afat3,
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

    atex2b$other$userVarValues <- NULL
    atex2b$PerSampleStats <- NULL
    atex2b$other$interventionTimes <- NULL
    save(file = "../../data/atex2b.RData", atex2b)
})


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
                            finalTime = 105,
                            mu = 1e-4,
                            initSize = 5000, 
                            keepPhylog = FALSE,
                            seed = NULL, 
                            errorHitMaxTries = FALSE, 
                            errorHitWallTime = FALSE,
                            userVars = userVars,
                            rules = rules,
                            keepEvery = 1)

    uvex3$PerSampleStats <- NULL
    uvex3$other$interventionTimes <- NULL
    save(file = "../../data/uvex3.RData", uvex3)

})



## ## userVarsOncoSimulIndivExample
## local({

##     userVars <- list(
##         list(Name           = "user_var1",
##              Value       = 0
##              ),
##         list(Name           = "user_var2",
##              Value       = 3
##              ),
##         list(Name           = "user_var3",
##              Value       = 2.5
##              )
##     )

##     rules <- list(
##         list(ID = "rule_1",
##              Condition = "T > 20",
##              Action = "user_var_1 = 1"
##              ),
##         list(ID = "rule_2",
##              Condition = "T > 30",
##              Action = "user_var_2 = 2; user_var3 = 2*N"
##              ),
##         list(ID = "rule_3",
##              Condition = "T > 40",
##              Action = "user_var_3 = 3;user_var_2 = n_A*n_B"
##              )
##     )
##     dfuv <- data.frame(Genotype = c("WT", "A", "B"),
##                        Fitness = c("1",
##                                    "1 + 0.2 * (n_B > 0)",
##                                    ".9 + 0.4 * (n_A > 0)"
##                                    ))
##     afuv <- allFitnessEffects(genotFitness = dfuv,
##                               frequencyDependentFitness = TRUE,
##                               frequencyType = "abs")

##     userVars <- createUserVars(userVars)
##     rules <- createRules(rules, afuv)

##     uvex <- oncoSimulIndiv(
##         afuv, 
##         model = "McFLD",
##         mu = 1e-4,
##         sampleEvery = 0.001,
##         initSize = c(20000, 20000),
##         initMutant = c("A", "B"),
##         finalTime = 5.2,
##         onlyCancer = FALSE,
##         userVars = userVars,
##         rules = rules
##     )

##     uvex$other$userVarValues <- NULL
##     uvex$PerSampleStats <- NULL
##     uvex$other$interventionTimes <- NULL
##     save(file = "../../data/userVarsOncoSimulIndivExample.RData", uvex)
## })


## ## interventionsOncoSimulIndivExample

## local({

##     interventions <- list(
##         list(ID           = "i2",
##              Trigger       = "(N > 1e6) & (T > 100)",
##              WhatHappens   = "N = 0.001 * N",
##              Repetitions   = 7,  
##              Periodicity    = Inf
##              ),
##         list(ID           = "i1",
##              Trigger       = "(T > 10)",
##              WhatHappens   = "N = 0.3 * N",
##              Periodicity   = 10,
##              Repetitions   = 0
##              ),
##         list(ID           = "i3", 
##              Trigger       = "(T > 1) & (T < 200)",
##              WhatHappens   = "n_A = n_A * 0,3 / n_C",
##              Repetitions   = Inf,  
##              Periodicity    = 10
##              ),
##         list(ID           = "i5",
##              Trigger       = "(N > 1e8) & (T> 1.2)",
##              WhatHappens   = "n_A_B = n_B * 0,3 / n_SRL",
##              Repetitions   = 0,   
##              Periodicity    = Inf
##              )
##     )


##     fa1 <- data.frame(Genotype = c("WT", "A", "B"),
##                       Fitness = c("n_*0",
##                                   "1.5",
##                                   "1"))

##     afd3 <- allFitnessEffects(genotFitness = fa1,
##                               frequencyDependentFitness = TRUE,
##                               frequencyType = "abs")

##     interventions <- createInterventions(interventions, afd3)
##     ep2 <- oncoSimulIndiv(
##         afd3, 
##         model = "Exp",
##         mu = 1e-4,
##         sampleEvery = 0.001,
##         initSize = c(20000, 20000),
##         initMutant = c("A", "B"),
##         finalTime = 5.2,
##         onlyCancer = FALSE,
##         interventions = interventions
##     )

##     ep2$other$userVarValues <- NULL
##     ep2$PerSampleStats <- NULL
##     ep2$other$interventionTimes <- NULL
##     save(file = "../../data/interventionsOncoSimulIndivExample.RData", ep2)
## })




