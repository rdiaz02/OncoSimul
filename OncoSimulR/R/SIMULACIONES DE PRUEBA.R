# PRUEBAS SIMULACIÓN CON OncoSimulR
# Autor: Germán Vallejo Palma.
# Simulación interacción fago bacteria.
library(OncoSimulR)
# 1º Crear un objeto allFitnessEffects frequencyDependentFitnes=TRUE

# a: tasa crecimiento bacterias
a <- 1.5
# b: tasa de "muerte" fagos
b <- 5
# c: tasa de infección por fagos
c <- 0.00001
# d: tasa de infección por fago y replicación
d <- 500000*c


genotipes_and_formulas <- data.frame("GENOTYPE"=c("B","F"),
                                     "FITNESS"=c(paste("f_B*",a,"-f_B*f_F*",c),
                                                 paste("f_B*f_F*",d,"-f_F*",b)))

fitnessbac_fag <- allFitnessEffects(genotFitness = genotipes_and_formulas,
                                    frequencyDependentFitness = TRUE,
                                    frequencyType = "rel")
evalGenotype("B",fitnessbac_fag,spPopSizes = c("B"=1,"F"=1))
