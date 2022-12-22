#SE TIENEN 4 POBLACIONES WT, A, B, C 
dfat <- data.frame(Genotype = c("WT","A", "B", "C"),
                   Fitness = c("0*n_+1.9",
                               "1.8",
                               "1.8",
                               "1.8"))

afe2 <- allFitnessEffects(genotFitness = dfat,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

evalAllGenotypes(afe2,spPopSizes = c(900,33,33,33))


intervenciones <- list(
  list(ID="ANTIBIOTICO SOBRE WT",
       Trigger       = "T>1",
       WhatHappens   = "n_ = n_ -x*n_*a0 -y*n_*b0+d0*a1*b1+e*x*y ",
       Periodicity   = 1,
       Repetitions   = Inf),
  list(ID="ANTIBIOTICO SOBRE A",
       Trigger       = "T>1",
       WhatHappens   = "n_A = n_A-x*n_A*a1-y*n_A*b1+d1*a1*b1+e*x*y",
       Periodicity   = 1,
       Repetitions   = Inf),
  list(ID="ANTIBIOTICO SOBRE B",
       Trigger       = "T>1",
       WhatHappens   = "n_B = n_B-x*n_B*a1-y*n_B*b1+d0*a1*b1+e*x*y",
       Periodicity   = 1,
       Repetitions   = Inf),
  list(ID="ANTIBIOTICO SOBRE A,B",
       Trigger       = "T>1",
       WhatHappens   = "n_C = n_C-x*n_C*a1-y*n_C*b1+d2*a1*b1+e*x*y",
       Periodicity   = 1,
       Repetitions   = Inf))



variables_modelo <- list(
  list(Name = "a0",
       Value = 1
  ),
  list(Name = "a1",
       Value = 0.6
  ),
  list(Name = "b0",
       Value  = 1
  ),
  list(Name = "b1",
       Value = 0.5
  ),
  list(Name = "x",
       Value = 0.5
  ),
  list(Name = "y",
       Value = 0.5
  ),
  list(Name= "d0",
       Value= 0
  ),
  list(Name= "d1",
       Value= 1
  ),
  list(Name= "d2",
       Value= -1
  ),
  list(Name= "e",
       Value= 0
  )
)




v_Model <- createUserVars(variables_modelo)  

inter2 <- createInterventions(intervenciones,afe2)

simu2 <- oncoSimulIndiv(afe2,
                        initMutant = c("WT", "B", "A", "C"),
                        initSize = c(970,10,10,10),
                        finalTime = 40,
                        mu=0.0000000001,
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



