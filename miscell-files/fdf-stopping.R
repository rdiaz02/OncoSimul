## Hurlbut et al., 2018, example

options(stringsAsFactors = FALSE) ## Get rid of the messages

create_fe <- function(a, b, c, d, e, f, g,
                        gt = c("WT", "A", "P", "C")) {
  data.frame(Genotype = gt,
             Fitness = c(
                 paste0("1 + ",
                       as.character(d), " * f_1 ",
                       "- ", as.character(c), " * f_3"),
                 paste0("1 - ", as.character(a), " + ", as.character(d), " + ",
                       as.character(f), " * f_1 ",
                      "- ", as.character(c), " * f_3"),
                 paste0("1 + ", as.character(g), " + ",
                       as.character(d), " * f_1 ",
                       "- ", as.character(c), " * (1 + ",
                       as.character(g), ") * f_3"),
                 paste0("1 - ", as.character(b), " + ",
                       as.character(e), " * f_ + ",
                       "(", as.character(d), " + ", as.character(e), ") * f_1 + ",
                       as.character(e) , " * f_2")),
             stringsAsFactors = FALSE)
}


afe_3_a <- allFitnessEffects(genotFitness =
                                 create_fe(0.02, 0.04, 0.08, 0.06,
                                           0.15, 0.1, 0.06),
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

set.seed(2)
s_3_a <- oncoSimulIndiv(afe_3_a,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 200,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)

set.seed(2)
s_3_a_s <- oncoSimulIndiv(afe_3_a,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     detectionProb = c(p2 = 0.2, 
                                       n2 = 5000, 
                                       PDBaseline = 2000,
                                       checkSizePEvery = 2),
                     ## finalTime = 200,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)

set.seed(2)
s_3_a_s_2 <- oncoSimulIndiv(afe_3_a,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     detectionProb = c(p2 = 0.05, 
                                       n2 = 5000, 
                                       PDBaseline = 2000,
                                       checkSizePEvery = 2),
                     ## finalTime = 200,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)


set.seed(2)
s_3_a_s_3 <- oncoSimulIndiv(afe_3_a,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     detectionProb = c(p2 = 0.2, 
                                       n2 = 5000, 
                                       PDBaseline = 4950,
                                       checkSizePEvery = 2),
                     ## finalTime = 200,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)

set.seed(2)
s_3_a_s_4 <- oncoSimulIndiv(afe_3_a,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     detectionProb = c(p2 = 0.1, 
                                       n2 = 5000, 
                                       PDBaseline = 4950,
                                       checkSizePEvery = 2),
                     ## finalTime = 200,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)

s_3_a_s
s_3_a_s_2
s_3_a_s_3
s_3_a_s_4


## 2a

a <- b <- c <- d <- e <- f <- g <- 0.1

afe_2_a <- allFitnessEffects(genotFitness = create_fe(a, b, c, d, e, f, g),
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

set.seed(1)
s_2_a <- oncoSimulIndiv(afe_2_a,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 80,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)

set.seed(1)
s_2_a_2 <- oncoSimulIndiv(afe_2_a,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 200,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)

set.seed(1)
s_2_a_3 <- oncoSimulIndiv(afe_2_a,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 2000,
                     fixation= "P",
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)

set.seed(1)
s_2_a_4 <- oncoSimulIndiv(afe_2_a,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 2000,
                     fixation= c("P", fixation_tolerance = 0.2),
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)

s_2_a
s_2_a_2
s_2_a_3
s_2_a_4
