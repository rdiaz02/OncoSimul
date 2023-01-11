

library(OncoSimulR)
library(HH)



# DEFINIR EL NUMERO DE ITERACIONES PARA CADA CASO Y SI MOSTRAR O NO LOS PLOTS
iteraciones = 30
plot_while_runing = FALSE

# DEFINIR VALORES PARA LOS TIPOS DE INTERACCIÓN SOBRE LOS QUE ITERAR
vd <- c(0.1, 0.05, 0, -0.05, -0.1)
ve <- c(2, 1, 0, -1, -2)

# INICIALIZAR UN DATAFRAME DONDE SE ALMACENAN LOS RESULTADOS

GR <- data.frame(Grow_rate = 1,
                 INTERACTION = 1,
                 RESISTENCE = 1)


# CREAR UN OBJETO allFitnessEffects

genotFitness <- data.frame("genotype"= c("WT","R","NS"),
                           "fitness"= c("1.1", "1.05", "0*n_"),
                           stringsAsFactors = FALSE)

afe <- allFitnessEffects(genotFitness = genotFitness, 
                         frequencyDependentFitness = TRUE,
                         frequencyType = "abs")

# INICIALIZAR DOS BUCLES PARA ITERAR SOBRE vd Y ve

for (d_int in vd){

for (e_int in ve){
  
  
# DEFINICIÓN DE VARIABLES DE LAS INTERVENCIONES
  
# FACTOR DE RESISTENCIA BACTERIANA
a <- 0.3
b <- 0.3

# EFECTO INTERACCIÓN
d <-d_int
e <- e_int

# CONCENTRACIÓN ANTIBIOTICOS
x <- 0.05
y <- 0.05

# DURACIÓN DE CADA SIMULACIÓN
tiempo <- 1100 

# COMPROBAR QUE LAS INTERVENCIONES NO HACEN QUE AUMENTE LA POBLACIÓN

stopifnot( ( -x * a - y * b + d * a * b + e * x * y ) < 0 )  
  
# CREAR LAS INTERVENCIONES (EFECTO DE LOS ANTIBIOTICOS)

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

# CREAR LAS VARIABLES 

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
inteRvenciones <- createInterventions(intervenciones, afe)

# REPETIR CADA SIMULACIÓN PARA TENER MAYOR CANTIDAD DE DATOS

for (iteracion in 1:iteraciones){
  
# CORRER SIMULACIÓN
  
sim <- oncoSimulIndiv(afe, 
                     keepEvery = 1,
                     mu= 0.000000001,
                     initMutant= c("WT","R"),
                     initSize = c("WT"=90000,"R"=1000),
                     interventions = inteRvenciones,
                     userVars = vaRiables,
                     finalTime = tiempo)
  
# EXTRAER VALORES POBLACIONALES DE LA SIMULACIÓN

pobs <-  unlist(sim$pops.by.time)[,2:3]
totalpob <- rowSums(pobs)
time <- unlist(sim$pops.by.time)[,1]

# CALCULAR LA EVOLUCIÓN LOS VALORES POBLACIONALES DE ALFA Y BETA

alfa0 = (pobs[,1] * 1)/totalpob
alfa1 = (pobs[,2] * a)/totalpob

alfa_medios = alfa0 + alfa1

beta0 = (pobs[,1] * 1)/totalpob
beta1 = (pobs[,2] * a)/totalpob

beta_medios = beta0 + beta1

# CALCULAR EL VALOR DE TASA DE CRECIMIENTO A CADA MOMENTO

grow_rate <- totalpob
init_pob <- 10000

for (i in 1:length(totalpob)){
  
actual_pob <- totalpob[i]
grow_rate[i] <- actual_pob /init_pob
init_pob <- actual_pob
} 

# MOSTRAR LOS RESULTADOS DE CADA SIMULACIÓN A LO LARGO DE LAS ITERACIONES

if (plot_while_runing){
  
par( mfrow = c(2, 2) )
  
plot(sim, show="genotypes")

plot( grow_rate[52:length(grow_rate)]~time[52:length(grow_rate)], lty=3)

plot( beta_medios[52:length(alfa_medios)]~time[52:length(grow_rate)] )

plot( alfa_medios[52:length(beta_medios)]~time[52:length(grow_rate)] )
}

# GUARDAR LA TASA DE CRECIMIENTO FINAL SI LA POBLACIÓN NO HA COLAPSADO

if (grow_rate[length(totalpob)]>1){
GR <- rbind(GR, c(grow_rate[length(totalpob)],
                  e_int,
                  d_int))}

}
}
}

# FIN DE ITERACIONES 

# NO SEGIR SI NO SE HAN ALMACENADO DATOS

stopifnot( length(GR[,1]) > 2 )


# GUARDAR TODOS LOS VALORES MENOS EL DE INICIALIZACIÓN
GR0 <- GR[2:length(GR[,1]),]

colnames(GR0) <- c("TASA_CRECIMIENTO","INTERACCIÓN","EFECTO_MUTACIÓN")

# ENCONTRAR OUTGROUPS Y ELIMINARLOS MANUALMENTE
par(mfrow=c(2,2))

plot(aov(TASA_CRECIMIENTO~INTERACCIÓN+EFECTO_MUTACIÓN, data=GR0))

GR1 <- GR0[-c(724,742),]

plot(aov(TASA_CRECIMIENTO~INTERACCIÓN+EFECTO_MUTACIÓN, data=GR1))



# RESUMEN DE LA VARIACIÓN DE LOS DATOS EN FUNCIÓN DE e Y d

colnames(GR1) <- c("TC","e","d")
interaction2wt(TC~e+d,
               data=GR1,main = "RESUMEN TASA DE CRECIMIENTO")



# MOSTRAR COMPARACIONES EN BOXPLOT
par( mfrow = c( 1, 3 ) )

colnames(GR1) <- c("TASA_CRECIMIENTO",
                   "INTERACCIÓN_ANTIBIÓTICO",
                   "EFECTO_MUTACIÓN")

boxplot(TASA_CRECIMIENTO~INTERACCIÓN_ANTIBIÓTICO,
        data=GR1[GR1[,3]==-0.1,],
        main="SENSIBILIDAD CRUZADA (d = -0.1)")

boxplot(TASA_CRECIMIENTO~INTERACCIÓN_ANTIBIÓTICO,
        data=GR1[GR1[,3]==0,],
        main="SIN CORRELACIÓN (d = 0)")

boxplot(TASA_CRECIMIENTO~INTERACCIÓN_ANTIBIÓTICO,
        data=GR1[GR1[,3]==0.1,],
        main="RESISTENCIA CRUZADA (d = 0.1)")


