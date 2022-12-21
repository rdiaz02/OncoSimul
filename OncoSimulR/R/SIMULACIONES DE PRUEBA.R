
dfat <- data.frame(Genotype = c("WT", "B", "A", "B, A"),
                   Fitness = c("0*n_+1.9",
                               "1.77",
                               "1.7",
                               "1.5"))

afe2 <- allFitnessEffects(genotFitness = dfat,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

evalAllGenotypes(afe2,spPopSizes = c(900,33,33,33))


intervenciones <- list(
  list(ID="ANTIBIOTICO SOBRE WT",
       Trigger       = "TRUE",
       WhatHappens   = "n_ = n_ -x*n_ -y*n_ ",
       Periodicity   = 1,
       Repetitions   = Inf),
 list(ID="ANTIBIOTICO SOBRE A",
      Trigger       = "TRUE",
      WhatHappens   = "n_A = n_A -x*n_A*a-y*n_A",
      Periodicity   = 1,
      Repetitions   = Inf),
 list(ID="ANTIBIOTICO SOBRE B",
      Trigger       = "TRUE",
      WhatHappens   = "n_B = n_B - x*n_B-y*n_B*b",
      Periodicity   = 1,
      Repetitions   = Inf),
 list(ID="ANTIBIOTICO SOBRE A,B",
      Trigger       = "TRUE",
      WhatHappens   = "n_A_B = n_A_B - x*n_A_B*a-y*n_A_B*b",
      Periodicity   = 1,
      Repetitions   = Inf))



variables_modelo <- list(
  list(Name = "a",
       Value = 0.4
  ),
  list(Name = "b",
       Value       = 0.5
  ),
  list(Name = "x",
       Value = 0.2
  ),
  list(Name = "y",
       Value = 0.3),
  list(Name = "Ma",
       Value = 1),
  list(Name = "nWT"
       Value= 9100),
  list(Name = "nA"
       Value= 300),
  list(Name = "nB"
       Value= 300),
  list(Name = "nAB"
       Value= 300),
  list(Name = "gWT"
       Value= 0),
  list(Name = "gA"
       Value= 0),
  list(Name = "gB"
       Value= 0),
  list(Name = "gAB"
       Value= 0))

rules <- list(
  list(ID="Calcular media de alfa Ma",
       Condition= "TRUE",
       Action= ),
  )
v_Model <- createUserVars(variables_modelo)  
                        
inter2 <- createInterventions(intervenciones,afe2)

simu2 <- oncoSimulIndiv(afe2,
               initMutant = c("WT", "B", "A", "B, A"),
               initSize = c(9100,300,300,300),
               finalTime = 20,
               interventions = inter2,
               userVars = v_Model,
               keepEvery = 1)
plot(simu2,show="genotypes")

pobs <- unlist(simu2$pops.by.time)[,2:5]
totalpob <- rowSums(unlist(simu2$pops.by.time)[,2:5])
freqs <- pobs/totalpob
time <- unlist(simu2$pops.by.time)[,1]
max <- totalpob/totalpob
plot(max~time,type="l", col="black",ylim=c(0,1.1),ylab="FREQUENCY")
lines(time,freqs[,1],type="l",col="#A6761D")
lines(time,freqs[,2],type="l",col="#666666")
lines(time,freqs[,3],type="l",col="#1B9E89")
lines(time,freqs[,4],type="l",col="red")

unlist(simu2$other)

