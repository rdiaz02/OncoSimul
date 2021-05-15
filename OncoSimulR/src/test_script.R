interventions <- list(
    list(ID           = "i2",
        Trigger       = "(N > 1e6) & (T > 100)",
        WhatHappens   = "N = 0.001 * N",
        Repetitions   = 7,  
        Periodicity    = Inf
    ),
    list(ID           = "i1",
        Trigger       = "(N > 1e9)",
        WhatHappens   = "N = 0.3 * N",
        Periodicity   = 10,
        Repetitions   = 0
    ),
    list(ID           = "i3",
        Trigger       = "(N > 1e9)",
        WhatHappens   = "n_A = n_A * 0,3 / n_C",
        Repetitions   = Inf,   ## This will be translated to MAX_INT
        Periodicity    = Inf
    ),
    list(ID           = "i5",
        Trigger       = "(N > 1e9)",
        WhatHappens   = "n_A = n_B * 0,3 / n_C",
        Repetitions   = 0,   
        Periodicity    = Inf
    ),
    list(ID           = "i4",
        Trigger       = "(N > 1e9)",
        WhatHappens   = "n_B   = n_A * 0,3 / n_C",
        Repetitions   = Inf,   ## This will be translated to MAX_INT
        Periodicity    = Inf
    )
)

df3x <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                      Fitness = c("n_A + 0.1",
                                  "log(n_A_B)",
                                  "sqrt(n_A) + 3 * exp(n_C)",
                                  "-1 * n_A + 7 * (n_C_A > 2)",
                                  "2 - 0.1/n_B",
                                  "min(n_C, n_B) - 1 * (n_B_A > 2)"))

adf3x <- allFitnessEffects(genotFitness = df3x,frequencyDependentFitness = TRUE)

interventions <- create_interventions(interventions, adf3x)

ep1 <- oncoSimulIndiv(adf3x, model = "McFL",
                      mu = 5e-6,
                      sampleEvery = 0.025,
                      keepEvery = 0.5,
                      initSize = 2000,
                      finalTime = 4,
                      onlyCancer = FALSE,
		              interventions = interventions)