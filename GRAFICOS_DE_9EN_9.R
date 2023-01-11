# Create an allFitnesEffects object with a WT and a resistant population.
library(OncoSimulR)
par(mfrow=c(3,3),cex.axis=0.8)
genotFitness <- data.frame("genotype"= c("WT","R","NS"),
                           "fitness"= c("1.1", "1.05", "0*n_"),
                           stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = genotFitness, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs")
plotFitnessLandscape(evalAllGenotypes(afe,spPopSizes = c(1,1,1)))

# INTERVENCIONES Y VARIABLES EN EL MODELO, EFECTO DEL ANTIBIÓTICO 
vd <- c(0.1,0,-0.1)
ve <- c(2,0,-2)
for (tipo_grafico in 4){
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
      
      a <- 0.35; b <- 0.25; d <-d_int; e <- it
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
      
      
      pobs <-  unlist(sim$pops.by.time)[,2:3]
      totalpob <- rowSums(pobs)
      TASA_CRECIMIENTO <- totalpob
      TIEMPO <- unlist(sim$pops.by.time)[,1]
      init_pob <- 10000
      FRECUENCIA<- pobs/totalpob
    
      
      if (tipo_grafico==1){
      plot(TIEMPO,
           ylim=c(0,max(pobs)),
           type="n",
           ylab="Frecuencia",
           xlab="Tiempo",
           xlim=c(0,tiempo))
      
      lines(TIEMPO,pobs[,1],col="red")
      lines(TIEMPO,pobs[,2],col="blue")
      
      }
      
      if (tipo_grafico==2){
      
      plot(totalpob~TIEMPO,xlim=c(0,1500),
           ylab ="Población total",
           xlab = "Tiempo")}
      
      if (tipo_grafico==3){
        
        tasa_crecimiento <- rep(1,(length(totalpob)-length(totalpob)%%5)/5)
        time <-  rep(1,(length(TIEMPO)-length(TIEMPO)%%5)/5)
        length(tiempo1)
        totalpob1 <- totalpob[1:(length(totalpob)-(length(totalpob)%%5))]
        tiempo1 <-  TIEMPO[1:(length(TIEMPO)-(length(TIEMPO)%%5))]
        last_pob<-10000
        
        for (z in seq(from=1,to=length(totalpob1),by=5)){
        actual_pob <- mean(totalpob1[z:(5+z)])
        tasa_crecimiento[z] <-actual_pob/last_pob 
        time[z] <- mean(tiempo1[z:(z+5)])
        last_pob <- actual_pob}
        
        plot(tasa_crecimiento~time,xlim=c(0,1500), 
             xlab="Tiempo", ylab="Tasa crecimiento")
        
      }
      
      if (tipo_grafico==4){
        alfa0= (pobs[,1]*1)/totalpob
        alfa1= (pobs[,2]*a)/totalpob
        alfa_medios= (alfa0+alfa1)
        beta0= (pobs[,1]*1)/totalpob
      beta1= (pobs[,2]*b)/totalpob
      beta_medios=beta0+beta1
      plot(TIEMPO,alfa_medios,col="orange",type="l",
           ylab="Alfa y Beta",
           xlab="Tiempo",
           xlim=c(0,tiempo),
           ylim=c(0,1))
      lines(TIEMPO,beta_medios,col="green",lty=3)
        
      }
      
      
    }
  }
}

  
