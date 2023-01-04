# Create an allFitnesEffects object with a WT and a resistant population.
library(OncoSimulR)
genotFitness <- data.frame("genotype"= c("WT","R","R2"),
                           "fitness"= c("1.11", "1.05", "1.05+0*n_"),
                           stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = genotFitness, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs")
intervenciones <- list(
  list(ID="ANTIBIOTICO SOBRE WT",
       Trigger       = "TRUE",
       WhatHappens   = "n_ = n_*(1+(-x-y+e*x*y))",
       Periodicity   = 1,
       Repetitions   = Inf),
  list(ID="ANTIBIOTICO SOBRE R",
       Trigger       = "TRUE",
       WhatHappens   = "n_R = n_R*(1+(-x*a-y*b+d*a*b+e*x*y))",
       Periodicity   = 1,
       Repetitions   = Inf),
  list(ID="ANTIBIOTICO SOBRE R2",
       Trigger       = "TRUE",
       WhatHappens   = "n_R2 = n_R2*(1+(-x*a2-y*b2+d*a2*b2+e*x*y))",
       Periodicity   = 1,
       Repetitions   = Inf))

a <- 0.5; b <- 0.2; a2 <- 0.1 ; b2 <- 0.6;  d <- -0.1; 
e <- -0.5 ;x <- 0.05; y <- 0.05; tiempo <- 3000
-x*a-y*b+d*a*b+e*x*y
stopifnot((-x*a-y*b+d*a*b+e*x*y)<0 & (-x*a2-y*b2+d*a2*b2+e*x*y)<0)
variables <- list(
  list(Name = "a",
       Value = a
  ),
  list(Name = "a2",
       Value = a2
  ),
  list(Name = "b",
       Value = b
  ),
  list(Name = "b2",
       Value = b2
  ),
  list(Name = "x",
       Value = x
  ),
  list(Name = "y",
       Value = y
  ),
  list(Name= "d",
       Value= d
  ),
  list(Name= "e",
       Value= e
  )
)

RULES <- list(list(ID="CAMBIO ANTIBIOTICO ANTI R2",
                   Condition="n_R2 > 2*n_R",
                   Action="x = 0.04 ; y = 0.07"),
              list(ID="CAMBIO ANTIBIOTICO ANTI R",
                   Condition="n_R > 2*n_R2",
                   Action="x = 0.07 ; y = 0.04"))
vaRiables <- createUserVars(variables)
inteRvenciones <- createInterventions(intervenciones,afe)
rules <- createRules(RULES,afe)

#SIMULACIONES
sim<- oncoSimulIndiv(afe, 
                     keepEvery = 1,
                     mu= 0.001,
                     initMutant= c("WT","R","R2"),
                     initSize = c("WT"=9900,"R"=100,"R2"=100),
                     interventions = inteRvenciones,
                     userVars = vaRiables,
                     finalTime = tiempo,
                     rules = rules)
par(mfrow=c(1,2))
plot(sim,show = "genotypes")
N_tot <- rowSums(unlist(sim$pops.by.time)[,2:4])
tiempo <- unlist(sim$pops.by.time)[,1]
plot(tiempo,N_tot,
     col="orange",
     main="Evolución población total",
     type="line",
     ylab="Células",
     xlab="TIEMPO")