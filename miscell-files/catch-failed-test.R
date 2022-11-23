

library(OncoSimulR)
library(testthat)

i <- round(runif(1, 1, 1e9))


while (TRUE) {
    i <- i + 1
    set.seed(i)
    cat("\n\n\n\n  ##################   seed = ", i, "\n")
    
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
        indexes <- vector()
                                        # We get the indexes that match the intervention times in pops.by.time
        for(time in sfd3$other$interventionTimes){
            indexes <- append(indexes, which(sfd3$pops.by.time[,1] == time))
        }

                                        # For each intervention time (T = 10, 20, 30)
        for(index in indexes){
            line <- sfd3$pops.by.time[index, ]
            prev_line <- sfd3$pops.by.time[index-1, ]
                                        #Total
            total <- line[2] + line[3] + line[4]
            prev_total <- prev_line[2] + prev_line[3] + prev_line[4]
            testthat::expect_gt(total, prev_total*0.8 - 0.2*prev_total)
            testthat::expect_lt(total, prev_total*0.8 + 0.2*prev_total)
                                        #Genotype WT
            if((prev_line[2] > 0) & (line[2] > 0)){
                testthat::expect_gte(prev_line[2], line[2])
            }
                                        #Genotype A
            ## FIXME this failed once in Linux
            if((prev_line[3] > 0) & (line[3] > 0)){
                testthat::expect_gte(prev_line[3], line[3])
            }
                                        #Genotype B
            if((prev_line[4] > 0) & (line[4] > 0)){
                testthat::expect_gte(prev_line[4], line[4])
            }
        }
    })
}
