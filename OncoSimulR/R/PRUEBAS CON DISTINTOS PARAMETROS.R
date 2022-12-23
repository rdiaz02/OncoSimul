#SE TIENEN 4 POBLACIONES WT, A, B, C 
# SEED 18,14,6,88,78,70,60,44,142,114,111,110,107,103,94,93
seed <- c(6,14,18,88,78,70,60,44,142,114,111,110,107,103,94,93)



{i=18
set.seed(i)
# ALFA Y BETA
a0 <- 1
a1 <- round(runif(1,0,1),2)
a2 <- round(runif(1,0,1),2)
b0 <- 1
b1 <- round(runif(1,0,1),2)
b2 <- round(runif(1,0,1),2)

# DELTA Y EPSILON

d0 <- 0
d1 <- round(runif(1,-0.5,0.5),2)
d2 <- round(runif(1,-0.5,0.5),2)
d3 <- round(runif(1,-0.5,0.5),2)
e <- round(runif(1,-0.5,0.5),2)

# CONCENTRACIÓN ANTIBIÓTICOS X, Y
x <- round(runif(1,0,1),2)
y <- round(runif(1,0,1),2)



dfat <- data.frame(Genotype = c("WT","A", "B"),
                   Fitness = c("0*n_+2.2",
                               "2",
                               "2"))

afe2 <- allFitnessEffects(genotFitness = dfat,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "abs")

evalAllGenotypes(afe2,spPopSizes = c(990,5,5))


intervenciones <- list(
  list(ID="ANTIBIOTICO SOBRE WT",
       Trigger       = "T>0",
       WhatHappens   = paste("n_ = n_-",x,"*",a0,"*n_-",y,"*n_*",b0,"+(",
                             d0,"*",a0,"*",b0,")+(",e,"*",x,"*",y,")"),
       Periodicity   = 1,
       Repetitions   = Inf),
  list(ID="ANTIBIOTICO SOBRE A",
       Trigger       = "T>0",
       WhatHappens   = paste("n_A = n_A-",x,"*",a1,"*n_A-",y,"*n_B*",b1,
                             "+",d1,"*",a1,"*",b1,"+",e,"*",x,"*",y),
       Periodicity   = 1,
       Repetitions   = Inf),
  list(ID="ANTIBIOTICO SOBRE B",
       Trigger       = "T>0",
       WhatHappens   = paste("n_B = n_B-",x,"*",a2,"*n_B-",y,"*n_B*",
                             b2,"+",d2,"*",a2,"*",b2,"+",e,"*",x,"*",y),
       Periodicity   = 1,
       Repetitions   = Inf)
)
  


inter2 <- createInterventions(intervenciones,afe2)
t <- 30 # TIEMPO DE SIMULACIÓN

simu2 <- oncoSimulIndiv(afe2,
                        initMutant = c("WT", "B", "A"),
                        initSize = c(970,10,10),
                        finalTime = t,
                        mu=0.0000000001,
                        interventions = inter2,
                        keepEvery = 1)
plot(simu2,show="genotypes",main=paste("SEED",i))

pobs <- unlist(simu2$pops.by.time)[,2:4]
totalpob <- rowSums(unlist(simu2$pops.by.time)[,2:4])
freqs <- pobs/totalpob
time <- unlist(simu2$pops.by.time)[,1]
max <- totalpob/totalpob
par(mfrow = c(1, 2))
plot(max~time,
     main=paste("SEED",i),
     type="n",
     ylim=c(0,1.1),xlim=c(0,t),
     ylab="FREQUENCY",
     xlab="TIME")
legend(18,1.15,
       legend=c("WT","A","B"),
       col=c("RED","BLUE","ORANGE"))
lines(time,freqs[,1],type="l",col="BLUE")
lines(time,freqs[,2],type="l",col="GREEN")
lines(time,freqs[,3],type="l",col="RED")

alfa0= (pobs[,1]*a0)/totalpob
alfa1= (pobs[,2]*a1)/totalpob
alfa2= (pobs[,3]*a2)/totalpob
alfa_medios= (alfa0+alfa1+alfa2)
plot(time,alfa_medios,type="l",
     main=paste("ALFA MEDIA - SEED",i),
     xlim=c(0,t),ylim=c(0,1),col="red",
     ylab="ALFA BETA MEDIOS")
beta0= (pobs[,1]*b0)/totalpob
beta1= (pobs[,2]*b1)/totalpob
beta2= (pobs[,3]*b2)/totalpob
beta_medios= (beta0+beta1+beta2)
lines(time,beta_medios,type="l",col="blue")
}


