library(OncoSimulR)

# GENERAR EL OBJETO ALLFITNESSEFECTS


genotFitness <- data.frame("genotype"= c("WT","R","NS"),
                           "fitness"= c("1.1", "1.05", "0*n_"),
                           stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = genotFitness, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs")



# ITERAR SOBRE DISTINTOS SUPUESTOS DE INTERACCIÓN UNA VEZ POR TIPO DE GRÁFICA

par(mfrow=c(3,3),cex.axis=0.8)

# VALORES EFECTOS COLATERALES DE MUTACIÓN
vd <- c(0.1,0,-0.1)

# VALORES INTERACCIÓN ANTIBIÓTICOS
ve <- c(2,0,-2)

for (tipo_grafico in 1:4){
for (d_int in vd){
for (it in ve){
  
# ASIGNAR VALORES A LAS VARIABLES  
a <- 0.35; b <- 0.25; d <-d_int; e <- it
x <- 0.05; y <- 0.05; tiempo <- 1500

# Comprobar que el efecto del anribiótico no aumenta la población
stopifnot( ( -x * a - y * b + d * a * b + e * x * y ) < 0 )  

# CREAR OBJETOS DE userVariables Y intervencitions

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
      
vaRiables <- createUserVars(variables)
inteRvenciones <- createInterventions(intervenciones,afe)

# CORRER SIMULACIÓN

sim<- oncoSimulIndiv(afe, 
                     keepEvery = 1,
                     mu= 0.000000001,
                     initMutant= c("WT","R"),
                     initSize = c("WT"=9900,"R"=100),
                     interventions = inteRvenciones,
                     userVars = vaRiables,
                     finalTime = tiempo)

# EXTRAER LOS DATOS DE LA SIMULACIÓN

pobs <-  unlist(sim$pops.by.time)[,2:3]
totalpob <- rowSums(pobs)

TASA_CRECIMIENTO <- totalpob # SE INICIALIZA LA TASA DE CRECIMIENTO

TIEMPO <- unlist(sim$pops.by.time)[,1]


FRECUENCIA<- pobs/totalpob


      
# MOSTRAR LA FRECUENCIA DE LAS SUBPOBLACIONES

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

# MOSTRAR LA POBLACIÓN TOTAL
      
if (tipo_grafico==2){

plot(totalpob~TIEMPO,xlim=c(0,1500),
     ylab ="Población total",
     xlab = "Tiempo")}

# MOSTRAR LA TASA DE CRECIMIENTO

if (tipo_grafico == 3){

#CALCULAR LA TASA DE CRECIMIENTO AGRUPANDO LOS VALORES DE 5 EN 5

# INICIALIZAR TASA DE CRECIMIENTO Y EL TIEMPO ASOCIADO
  
tasa_crecimiento <- rep(1,(length(totalpob) - length(totalpob) %% 5)/ 5)

time <-  rep(1,(length(TIEMPO) - length(TIEMPO) %% 5)/ 5)

# HACER QUE LA LONGITUD DE LOS DATOS SEA DIVISIBLE POR 5

totalpob1 <- totalpob[1:(length(totalpob)-(length(totalpob)%%5))]

tiempo1 <-  TIEMPO[1:(length(TIEMPO)-(length(TIEMPO)%%5))]

last_pob<-10000
        
for (z in seq(from=1,to=length(totalpob1),by=5)){
  
actual_pob <- mean(totalpob1[z:(5 + z)])

tasa_crecimiento[z] <-actual_pob/last_pob 

time[z] <- mean(tiempo1[z:(z + 5)])

last_pob <- actual_pob}

plot(tasa_crecimiento~time,xlim=c(0,1500), 
     xlab="Tiempo", ylab="Tasa crecimiento")
}

# MOSTRAR LOS FACTORES ALFA Y BETA

if (tipo_grafico==4){

# CALCULAR ALFA
  
alfa0= (pobs[,1] * 1)/totalpob
alfa1= (pobs[,2] * a)/totalpob
alfa_medios= alfa0 + alfa1

# CALCULAR BETA

beta0= (pobs[,1] * 1)/totalpob
beta1= (pobs[,2] * b)/totalpob
beta_medios = beta0 + beta1

# PLOT ALFA

plot(TIEMPO,alfa_medios,col="orange",
   type="l",
   ylab="Alfa y Beta",
   xlab="Tiempo",
   xlim=c(0,tiempo),
   ylim=c(0,1))

# AÑADIR BETA
lines(TIEMPO,beta_medios,col="green", lty=3)

}
}
}
}

  
