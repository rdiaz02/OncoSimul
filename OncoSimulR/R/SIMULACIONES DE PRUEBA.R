# PRUEBAS SIMULACIÓN CON OncoSimulR
# Autor: Germán Vallejo Palma.
# Simulación interacción fago bacteria.
library(OncoSimulR)
# 1º Crear un objeto allFitnessEffects frequencyDependentFitnes=TRUE

# a: tasa crecimiento bacterias
a <- 10
# b: tasa de "muerte" fagos
b <- 50
# c: tasa de infección por fagos
c <- 0.01
# d: tasa de infección por fago y replicación
d <- 500


genotipes_and_formulas <- data.frame("GENOTYPE"=c("B","F"),"FITNESS"=c("f_B*10-f_B*f_F*0.01","f_B*f_F*500-f_F*50"))

fitnessbac_fag <- allFitnessEffects(genotFitness = genotipes_and_formulas, frequencyDependentFitness = TRUE,frequencyType = "rel")
evalAllGenotypes(fitnessbac_fag,spPopSizes = c(100000,1000),model="exp",addwt=TRUE)
