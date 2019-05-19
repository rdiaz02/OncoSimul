
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


## Different assumption about origins from mutation:
## WT -> P; P -> A,P; P -> C,P

create_fe2 <- function(a, b, c, d, e, f, g,
                        gt = c("WT", "A", "P", "C", "A, P", "A, C", "P, C")) {
  data.frame(Genotype = gt,
             Fitness = c(
                 paste0("1 + ",
                       as.character(d), " * f_1_2 ",
                       "- ", as.character(c), " * f_2_3"),
                 "0",
                 paste0("1 + ", as.character(g), " + ",
                       as.character(d), " * f_1_2 ",
                       "- ", as.character(c), " * (1 + ",
                       as.character(g), ") * f_2_3"),
                 "0",
                 paste0("1 - ", as.character(a), " + ", as.character(d), " + ",
                       as.character(f), " * f_1_2 ",
                       "- ", as.character(c), " * f_2_3"),
                 "0",
                 paste0("1 - ", as.character(b), " + ",
                       as.character(e), " * f_ + ",
                       "(", as.character(d), " + ", as.character(e), ") * f_1_2 + ",
                       as.character(e) , " * f_2")),
             stringsAsFactors = FALSE)
}


## Show the expressions, as such
create_fe("a", "b", "c", "d", "e", "f", "g")
create_fe2("a", "b", "c", "d", "e", "f", "g")








## Figure 2a in paper.
a <- b <- c <- d <- e <- f <- g <- 0.1

afe_2_a <- allFitnessEffects(genotFitness = create_fe(a, b, c, d, e, f, g),
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")


s_2_a <- oncoSimulIndiv(afe_2_a,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 100,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s_2_a, show = "genotypes",
     xlim = c(40, 100))


sp_2_a <- oncoSimulPop(4,
                   afe_2_a,
                   model = "McFL", 
                   onlyCancer = FALSE, 
                   finalTime = 100,
                   mu = 1e-4,
                   initSize = 5000, 
                   keepPhylog = TRUE,
                   seed = NULL, 
                   errorHitMaxTries = FALSE, 
                   errorHitWallTime = FALSE,
                   mc.cores = 4)

plot(sp_2_a, show = "genotypes",
     xlim = c(40, 100))

## ## Very similar to 2a. Commented out
## ## Figure 2b in paper.
## a <- 0.02
## b <- 0.02
## c <- 0.11
## d <- 0.1
## e <- 0.1
## f <- 0.1
## g <- 0.15

## afe_2_b <- allFitnessEffects(genotFitness = create_fe(a, b, c, d, e, f, g),
##                          frequencyDependentFitness = TRUE,
##                          frequencyType = "rel")

## sp_2_b <- oncoSimulPop(4,
##                    afe_2_b,
##                    model = "McFL", 
##                    onlyCancer = FALSE, 
##                    finalTime = 100,
##                    mu = 1e-4,
##                    initSize = 5000, 
##                    keepPhylog = TRUE,
##                    seed = NULL, 
##                    errorHitMaxTries = FALSE, 
##                    errorHitWallTime = FALSE,
##                    mc.cores = 4)

## plot(sp_2_b, show = "genotypes",
##      xlim = c(1, 100))


## Figure 3a
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
## plot(s_3_a, show = "genotypes",
##      xlim = c(40, 200))
plot(s_3_a, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue"))


set.seed(3)
sp_3_a <- oncoSimulPop(10,
                   afe_3_a,
                   model = "McFL", 
                   onlyCancer = FALSE, 
                   finalTime = 200,
                   mu = 1e-4,
                   initSize = 10000, 
                   keepPhylog = TRUE,
                   seed = NULL, 
                   errorHitMaxTries = FALSE, 
                   errorHitWallTime = FALSE,
                   mc.cores = 4)
plot(sp_3_a, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue"))

plot(sp_3_a, show = "genotypes")



## Figure 3b
afe_3_b <- allFitnessEffects(genotFitness =
                                 create_fe(0.02, 0.04, 0.08, 0.1,
                                           0.15, 0.1, 0.05),
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")
set.seed(2)
s_3_b <- oncoSimulIndiv(afe_3_b,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 200,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
## plot(s_3_b, show = "genotypes", col = c("black", "green", "red", "blue"))
plot(s_3_b, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue"))


sp_3_b <- oncoSimulPop(10,
                   afe_3_b,
                   model = "McFL", 
                   onlyCancer = FALSE, 
                   finalTime = 200,
                   mu = 1e-4,
                   initSize = 5000, 
                   keepPhylog = TRUE,
                   seed = NULL, 
                   errorHitMaxTries = FALSE, 
                   errorHitWallTime = FALSE,
                   mc.cores = 8)
plot(sp_3_b, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue"))

plot(sp_3_b, show = "genotypes")




## Figure 3c

afe_3_c <- allFitnessEffects(genotFitness =
                                 create_fe(0.02, 0.04,
                                           0.01, 0.1,
                                           0.06, 0.05, 0.06),
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

s_3_c <- oncoSimulIndiv(afe_3_c,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 200,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s_3_c, show = "genotypes",
     xlim = c(40, 200))

sp_3_c <- oncoSimulPop(10,
                   afe_3_c,
                   model = "McFL", 
                   onlyCancer = FALSE, 
                   finalTime = 400,
                   mu = 1e-4,
                   initSize = 5000, 
                   keepPhylog = TRUE,
                   seed = NULL, 
                   errorHitMaxTries = FALSE, 
                   errorHitWallTime = FALSE,
                   mc.cores = 8)

plot(sp_3_c, show = "genotypes")



## Figure 3d

afe3_d <- allFitnessEffects(genotFitness =
                                 create_fe(0.02, 0.04,
                                           0.01, 0.1,
                                           0.12, 0.05, 0.05),
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

s3_d <- oncoSimulIndiv(afe3_d,
                     model = "McFL", 
                     onlyCancer = FALSE, 
                     finalTime = 200,
                     mu = 1e-4,
                     initSize = 5000, 
                     keepPhylog = TRUE,
                     seed = NULL, 
                     errorHitMaxTries = FALSE, 
                     errorHitWallTime = FALSE)
plot(s3_d, show = "genotypes",
     xlim = c(40, 200))

sp3_d <- oncoSimulPop(10,
                   afe3_d,
                   model = "McFL", 
                   onlyCancer = FALSE, 
                   finalTime = 400,
                   mu = 1e-4,
                   initSize = 5000, 
                   keepPhylog = TRUE,
                   seed = NULL, 
                   errorHitMaxTries = FALSE, 
                   errorHitWallTime = FALSE,
                   mc.cores = 8)

plot(sp3_d, show = "genotypes")



### Compare changing inheritance.



## Figure 3b. This one is interesting.
afe_3_b_2 <- allFitnessEffects(genotFitness =
                                 create_fe2(0.02, 0.04, 0.08, 0.1,
                                           0.15, 0.1, 0.05),
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")
set.seed(2)
s_3_b_2 <- oncoSimulIndiv(afe_3_b_2,
                   model = "McFL", 
                   onlyCancer = FALSE, 
                   finalTime = 300,
                   mu = 1e-4,
                   initSize = 5000, 
                   keepPhylog = TRUE,
                   seed = NULL, 
                   errorHitMaxTries = FALSE, 
                   errorHitWallTime = FALSE)
plot(s_3_b_2, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue"))

sp_3_b_2 <- oncoSimulPop(10,
                   afe_3_b_2,
                   model = "McFL", 
                   onlyCancer = FALSE, 
                   finalTime = 500,
                   mu = 1e-4,
                   initSize = 5000, 
                   keepPhylog = TRUE,
                   seed = NULL, 
                   errorHitMaxTries = FALSE, 
                   errorHitWallTime = FALSE,
                   mc.cores = 8)

plot(sp_3_b_2, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue"))

plot(sp_3_b_2, show = "genotypes")



## Figure 3 c

afe_3_c_2 <- allFitnessEffects(genotFitness =
                                 create_fe2(0.02, 0.04,
                                           0.01, 0.1,
                                           0.06, 0.05, 0.06),
                         frequencyDependentFitness = TRUE,
                         frequencyType = "rel")

sp_3_c_2 <- oncoSimulPop(10,
                   afe_3_c_2,
                   model = "McFL", 
                   onlyCancer = FALSE, 
                   finalTime = 400,
                   mu = 1e-4,
                   initSize = 5000, 
                   keepPhylog = TRUE,
                   seed = NULL, 
                   errorHitMaxTries = FALSE, 
                   errorHitWallTime = FALSE,
                   mc.cores = 8)

plot(sp_3_c_2, show = "genotypes")



## For vignette, just show 3a and 3b and 3b2
