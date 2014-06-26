
## ----echo=FALSE,results='hide',error=FALSE-------------------------------
require(knitr, quietly = TRUE)
opts_knit$set(concordance = TRUE)
## opts_knit$set(stop_on_error = 2L)


## ------------------------------------------------------------------------
library(OncoSimulR)
library(graph)


## ----fig.height=3--------------------------------------------------------
## Node 2 and 3 depend on 1, and 4 depends on no one
p1 <- cbind(c(1, 1, 0), c(2, 3, 4))
plotPoset(p1, addroot = TRUE)


## ----fig.height=3--------------------------------------------------------
## A simple way to create a poset where no gene (in a set of 15) depends
## on any other.
p4 <- cbind(0, 15)
plotPoset(p4, addroot = TRUE)


## ------------------------------------------------------------------------
## use poset p1101
data(examplePosets)
p1101 <- examplePosets[["p1101"]]

## Bozic Model
b1 <- oncoSimulIndiv(p1101, keepEvery = 15)
summary(b1)


## ----fig.height=5, fig.width=5-------------------------------------------
b2 <- oncoSimulIndiv(p1101, keepEvery = 1)
                    
summary(b2)
plot(b2)


## ------------------------------------------------------------------------

m2 <- oncoSimulIndiv(examplePosets[["p1101"]], model = "McFL", 
                     numPassengers = 0, detectionDrivers = 10, 
                     mu = 5e-7, initSize = 4000, 
                     sampleEvery = 0.025,
                     finalTime = 25000, keepEvery = 5, 
                     detectionSize = 1e6) 
plot(m2, addtot = TRUE, log = "")



## ------------------------------------------------------------------------
b3 <- oncoSimulIndiv(p1101, onlyCancer = FALSE)
summary(b3)

b4 <- oncoSimulIndiv(p1101, onlyCancer = FALSE)
summary(b4)


## ----fig.width=8, fig.height=4-------------------------------------------
par(mfrow = c(1, 2))
par(cex = 0.8) ## smaller font
plot(b3)
plot(b4)


## ------------------------------------------------------------------------
p1 <- oncoSimulPop(4, p1101)
par(mfrow = c(2, 2))
plot(p1)


## ------------------------------------------------------------------------

m1 <- oncoSimulPop(100, examplePosets[["p1101"]], 
                   numPassengers = 0)



## ------------------------------------------------------------------------
genotypes <- samplePop(m1)


## ----fig.width=4, fig.height=4-------------------------------------------
colSums(genotypes)/nrow(genotypes)

require(Oncotree)
ot1 <- oncotree.fit(genotypes)
plot(ot1)


## ----fig.width=4, fig.height=4-------------------------------------------
genotypesSC <- samplePop(m1, typeSample = "single")
colSums(genotypesSC)/nrow(genotypesSC)

ot2 <- oncotree.fit(genotypesSC)
plot(ot2)


