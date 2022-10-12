
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
    save("HansenExampleAT4.RData", file = atex5)
})


atex5 <- OncoSimulR:::thin.pop.data(atex5o, keep = 0.01, min.keep = 1)
print(object.size(atex5), units = "MB")
