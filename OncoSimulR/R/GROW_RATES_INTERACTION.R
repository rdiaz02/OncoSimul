# Create an allFitnesEffects object with a WT and a resistant population.

genotFitness <- data.frame("genotype"= c("WT","R","NS"),
                           "fitness"= c("1.1", "1.05", "0*n_"),
                           stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = genotFitness, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs")

# INTERVENCIONES Y VARIABLES EN EL MODELO, EFECTO DEL ANTIBIÓTICO 
vd <- c(-0.1,0,0.1)
GR <- data.frame(Grow_rate=1,INTERACTION=1, RESISTENCE=1)
for (d_int in vd){

ve <- c(2,0,-2)


for (it in ve){
intervenciones <- list(
  list(ID="ANTIBIOTICO SOBRE WT",
       Trigger       = "T>50",
       WhatHappens   = "n_ = n_*(1+(-x-y+e*x*y))",
       Periodicity   = 1,
       Repetitions   = Inf),
  list(ID="ANTIBIOTICO SOBRE R",
       Trigger       = "T>50",
       WhatHappens   = "n_R = n_R*(1+(-x*a-y*b+d*a*b+e*x*y))",
       Periodicity   = 1,
       Repetitions   = Inf))

a <- 0.3; b <- 0.3; d <-d_int; e <- it; x <- 0.05; y <- 0.05; tiempo <- 1000
-x*a-y*b+d*a*b+e*x*y
stopifnot((-x*a-y*b+d*a*b+e*x*y)<0)
variables <- list(
  list(Name = "a",
       Value = a
  ),
  list(Name = "b",
       Value = b
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

vaRiables <- createUserVars(variables)
inteRvenciones <- createInterventions(intervenciones,afe)

#SIMULACIONES
for (iteracion in 1:60){
sim<- oncoSimulIndiv(afe, 
                     keepEvery = 1,
                     mu= 0.000000001,
                     initMutant= c("WT","R"),
                     initSize = c("WT"=9900,"R"=100),
                     interventions = inteRvenciones,
                     userVars = vaRiables,
                     finalTime = tiempo)
  
# PLOTS

par(mfrow=c(1,2))

plot(sim, show="genotypes")

pobs <-  unlist(sim$pops.by.time)[,2:3]
totalpob <- rowSums(pobs)
grow_rate <- totalpob
time <- unlist(sim$pops.by.time)[,1]
init_pob <- 10000

for (i in 1:length(totalpob)){
actual_pob <- totalpob[i]
grow_rate[i] <- actual_pob /init_pob
init_pob <- actual_pob} 

plot(grow_rate~time,lty=3)
if (grow_rate[length(totalpob)]>1){
GR <- rbind(GR,c(grow_rate[length(totalpob)],it,d_int))
print(paste(it,d_int))}


}
}
}

par(mfrow=c(3,1))
GR <- GR[2:length(GR[,1]),]
colnames(GR) <- c("TASA_CRECIMIENTO","INTERACCIÓN","EFECTO_MUTACIÓN")
boxplot(TASA_CRECIMIENTO~INTERACCIÓN,
        data=GR[GR[,3]==-0.1,],
        main="SENSIBILIDAD CRUZADA (d = -0.1)")

boxplot(TASA_CRECIMIENTO~INTERACCIÓN,
        data=GR[GR[,3]==0,],
        main="EFECTO ADITIVO (d = 0)")

boxplot(TASA_CRECIMIENTO~INTERACCIÓN,
        data=GR[GR[,3]==0.1,],
        main="RESISTENCIA CRUZADA (d = 0.1)")


