# Create an allFitnesEffects object with a WT and a resistant population.

genotFitness <- data.frame("genotype"= c("WT","R","NS"),
                           "fitness"= c("1.1", "1.05", "0*n_"),
                           stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = genotFitness, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs")

# INTERVENCIONES Y VARIABLES EN EL MODELO, EFECTO DEL ANTIBIÓTICO 
par(mfrow=c(3,3))
vd <- c(-0.1,0,0.1)
ve <- c(-2,0,2)
for (tipo_grafico in 1:3){
for (d_int in vd){
  
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
      
      a <- 0.3; b <- 0.3; d <-d_int; e <- it
      x <- 0.05; y <- 0.05; tiempo <- 1500
      
      # Comprobar que el efecto del anribiótico no aumenta la población
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
      
      sim<- oncoSimulIndiv(afe, 
                           keepEvery = 1,
                           mu= 0.000000001,
                           initMutant= c("WT","R"),
                           initSize = c("WT"=9900,"R"=100),
                           interventions = inteRvenciones,
                           userVars = vaRiables,
                           finalTime = tiempo)
      
      # PLOTS
      
      
      if (tipo_grafico==1)plot(sim, show="genotypes")
      
      
      pobs <-  unlist(sim$pops.by.time)[,2:3]
      totalpob <- rowSums(pobs)
      TASA_CRECIMIENTO <- totalpob
      TIEMPO <- unlist(sim$pops.by.time)[,1]
      init_pob <- 10000
      FRECUENCIA<- pobs/totalpob
      
      if (tipo_grafico==2){
      plot(TIEMPO,
           ylim=c(0,max(pobs)),
           type="n",
           ylab="FRECUENCIA",
           xlab="TIEMPO",
           xlim=c(0,tiempo))
      
      lines(TIEMPO,pobs[,1],col="red")
      lines(TIEMPO,pobs[,2],col="blue")
      legend("topleft",
             legend = c("WT","R"),
             col= c("red", "blue" ),lwd=1, lty=c(1,1),
            title="Genotypes" )
      }
      par(cex=1)
      if (tipo_grafico==3){
      for (i in 1:length(totalpob)){
        actual_pob <- totalpob[i]
        TASA_CRECIMIENTO[i] <- actual_pob /init_pob
        init_pob <- actual_pob} 
      plot(totalpob~TIEMPO,xlim=c(0,1500))}
      
      #plot(TASA_CRECIMIENTO~TIEMPO,xlim=c(0,1000))
    }
  }
}

  

