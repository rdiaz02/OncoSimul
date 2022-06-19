# Ecology examples

## Predator-prey with carrying capacity

G_fe_LVm <- function(r1, r2, K1, K2, a_12, a_21, gt = c("S1", "S2")) {
    data.frame(Genotype = gt,
               Birth = c(paste0(r1, "-", r1, "*(", a_12, "*n_", gt[2], ")/", K1), r2),
               Death = c(paste0(r1, "*(n_", gt[1], ")/", K1),
                         paste0(r2, "*(n_", gt[2], "+", a_21, "*n_", gt[1], ")/", K2)))
}
fe_pred_prey <- allFitnessEffects(
		genotFitness = G_fe_LVm(1.4, 1.5, 4000, 10000, -0.5, 1.1, gt = c("Predator", "Prey")),
		frequencyDependentBirth = TRUE, 
		frequencyDependentDeath = TRUE,
		deathSpec = TRUE)

s_pred_preym <- oncoSimulIndiv(fe_pred_prey, model = "Arb",
                                initMutant = c("Predator", "Prey"),
                                initSize = c(1000, 1000), 
                                onlyCancer = FALSE, 
                                finalTime = 200, mu = 1e-3,
                                keepPhylog = TRUE, seed = NULL,
                                errorHitMaxTries = FALSE,
                                errorHitWallTime = FALSE)
								
## Standard predator-prey

C_fe_pred_prey2 <- function(r, a, c, e, d, gt = c("s1", "s2")) {
    data.frame(Genotype = gt,
               Birth = c(r, paste0(e, "*", a, "*", c, "*n_1")),
               Death = c(paste0(a, "*", c, "*n_2"), d))
}

fe_pred_prey2 <- allFitnessEffects(
		genotFitness = C_fe_pred_prey2(r = .8, a = 1, c = 0.007,
									   e = 0.1, d = 0.4,
									   gt = c("Fly", "Lizard")),
		frequencyDependentBirth = TRUE,
		frequencyDependentDeath = TRUE,
		deathSpec = TRUE)

set.seed(1)
pred_prey2 <- oncoSimulIndiv(fe_pred_prey2,
                             model = "Arb",
                             initMutant = c("Fly", "Lizard"),
                             initSize = c(500, 200),
                             mu = 1e-8,
                             onlyCancer = FALSE, 
                             finalTime = 100,
                             keepPhylog = TRUE,
                             seed = NULL, 
                             errorHitMaxTries = FALSE, 
                             errorHitWallTime = FALSE)
							 
plot(pred_prey2, show="genotypes", type="line", log = "")
# Hawks and doves

H_D_fitness <- function(c, v,
                    gt = c("H", "D")) {
  data.frame(Genotype = gt,
             Birth = c(
			   paste0("max(1e-5, f_H *", (v-c)/2, "+ f_D *", v, ")"),
               paste0("f_D *", v/2)))
}

## H = D

HD_eq <-allFitnessEffects(
	genotFitness = H_D_fitness(10, 4, gt = c("H", "D")),
	frequencyDependentBirth = TRUE,
	frequencyType = "rel")
	
osi_eq <- oncoSimulIndiv(HD_eq,
                           model = "Const",
                           onlyCancer = FALSE,
                           finalTime = 50,
                           mu = 1e-6,
                           initSize = c(2000, 2000),
						   initMutant = c("H", "D"),
                           keepPhylog = TRUE,
                           seed = NULL,
                           errorHitMaxTries = FALSE,
                           errorHitWallTime = FALSE)
plot(osi_eq, show="genotypes", ylim=c(1, 5000))


## H > D

HD_hawks <-allFitnessEffects(
	genotFitness = H_D_fitness(3, 4, gt = c("H", "D")),
	frequencyDependentBirth = TRUE,
	frequencyType = "rel")
	
osi_hawks <- oncoSimulIndiv(HD_hawks,
                           model = "Const",
                           onlyCancer = FALSE,
                           finalTime = 50,
                           mu = 1e-6,
                           initSize = c(2000, 2000),
						   initMutant = c("H", "D"),
                           keepPhylog = TRUE,
                           seed = NULL,
                           errorHitMaxTries = FALSE,
                           errorHitWallTime = FALSE)
						   
plot(osi_hawks, show="genotypes", ylim=c(1, 5000))

## H < D

HD_doves <-allFitnessEffects(
	genotFitness = H_D_fitness(30, 2, gt = c("H", "D")),
	frequencyDependentBirth = TRUE,
	frequencyType = "rel")
	
osi_doves <- oncoSimulIndiv(HD_doves,
                           model = "Const",
                           onlyCancer = FALSE,
                           finalTime = 50,
                           mu = 1e-6,
                           initSize = c(2000, 2000),
						   initMutant = c("H", "D"),
                           keepPhylog = TRUE,
                           seed = NULL,
                           errorHitMaxTries = FALSE,
                           errorHitWallTime = FALSE)
						   
plot(osi_doves, show="genotypes", ylim=c(1, 5000))

 
# Breast cancer example


create_fe <- function(cG, bG, cS, cMMph, cMMTC, bR, cD, initSize,
                      gt = c("WT", "Mph", "BTC", "MTC")) {
                      
    K = initSize/(exp(1)-1)
    data.frame(Genotype = gt,
               Birth = c(paste0(bG, "*(f_ + f_Mph)"),
                         paste0(bG, "*(f_ + f_Mph)"),
                         paste0(bR, " + ", bG, "*(f_ + f_Mph)"),
                         paste0(bR, " + ", bG, " *(f_ + f_Mph)")),
               Death = c(paste0(cG, " + ", cS, "*(f_ + f_BTC)+log(1+N/", K, ")"),
                         paste0(cG, " + ", cMMph, "+log(1+N/", K, ")"),
                         paste0(cS, "* (f_ + f_BTC) +", cD , " * f_Mph+log(1+N/", K, ")"),
                         paste0(cMMTC, " + ", cD , " * f_Mph+log(1+N/", K, ")")),
               stringsAsFactors = FALSE)
}

evalAllGenotypes(allFitnessEffects(genotFitness = create_fe(2, 5, 1, 0.8, 1, 1, 9), 
                                   frequencyDependentBirth = TRUE,
                                   frequencyDependentDeath = TRUE,
                                   deathSpec = TRUE,
                                   frequencyType = "rel"),
                 spPopSizes = c(WT = 1000, Mph = 0, BTC = 0, MTC = 0))
				 
evalAllGenotypes(allFitnessEffects(genotFitness = create_fe(2, 5, 1, 0.8, 1, 1, 9), 
                                   frequencyDependentBirth = TRUE,
                                   frequencyDependentDeath = TRUE,
                                   deathSpec = TRUE,
                                   frequencyType = "rel"),
                 spPopSizes = c(WT = 1000, Mph = 1000, BTC = 0, MTC = 0))
				 
evalAllGenotypes(allFitnessEffects(genotFitness = create_fe(2, 5, 1, 0.8, 1, 1, 9), 
                                   frequencyDependentBirth = TRUE,
                                   frequencyDependentDeath = TRUE,
                                   deathSpec = TRUE,
                                   frequencyType = "rel"),
                 spPopSizes = c(WT = 1000, Mph = 1000, BTC = 100, MTC = 100))


## Controlled cancer
afe_control <- allFitnessEffects(genotFitness = create_fe(0.5, 4, 1, 0.2, 1, 0.5, 4, 10000),
                                 frequencyDependentBirth = TRUE,
                                 frequencyDependentDeath = TRUE,
                                 deathSpec = TRUE,
                                 frequencyType = "rel")
								 
osi_control <- oncoSimulIndiv(afe_control,
                              model = "Arb",
                              onlyCancer = FALSE,
                              finalTime = 50,
                              mu = 1e-4,
                              initSize = 10000,
                              keepPhylog = TRUE,
                              seed = NULL,
                              errorHitMaxTries = FALSE,
                              errorHitWallTime = FALSE)

plot(osi_control, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue", "yellow"), xlim=c(0, 50))

## Non metastasic cancer
afe_no_met <- allFitnessEffects(genotFitness = create_fe(1, 4, 0.5, 1, 1.5, 1, 4, 10000),
                                frequencyDependentBirth = TRUE,
                                frequencyDependentDeath = TRUE,
                                deathSpec = TRUE,
                                frequencyType = "rel")
osi_no_met <- oncoSimulIndiv(afe_no_met,
                             model = "Arb", 
                             onlyCancer = FALSE, 
                             finalTime = 50,
                             mu = 1e-4,
                             initSize = 10000, 
                             keepPhylog = TRUE,
                             seed = NULL, 
                             errorHitMaxTries = FALSE, 
                             errorHitWallTime = FALSE)
plot(osi_no_met, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue", "yellow"))


## Metastasic cancer
afe_met <- allFitnessEffects(genotFitness = create_fe(0.5, 4, 2, 0.5, 0.5, 1, 4, 10000),
                             frequencyDependentBirth = TRUE,
                             frequencyDependentDeath = TRUE,
                             deathSpec = TRUE,
                             frequencyType = "rel")
							 
osi_met <- oncoSimulIndiv(afe_met,
                          model = "Arb",
                          onlyCancer = FALSE,
                          finalTime = 50,
                          mu = 1e-4,
                          initSize = 10000,
                          keepPhylog = TRUE,
                          seed = NULL,
                          errorHitMaxTries = FALSE,
                          errorHitWallTime = FALSE)
plot(osi_met, show = "genotypes", type = "line",
     col = c("black", "green", "red", "blue", "yellow"))