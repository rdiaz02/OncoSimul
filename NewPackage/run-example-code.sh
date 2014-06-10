#!/bin/bash
echo '

library(OncoSimulR)
data(examplePosets)
p701 <- examplePosets[["p701"]]
 set.seed(1)
 b1 <- oncoSimulIndiv(p701, keepEvery = 5,
                      sampleEvery = 5,
                      detectionSize = 1e5)
 print("########### changing seed")
 set.seed(2)
 b1 <- oncoSimulIndiv(p701, keepEvery = 5,
                      sampleEvery = 5,
                      detectionSize = 1e5)

 print("########### seed one again")
 set.seed(1)
 b1 <- oncoSimulIndiv(p701, keepEvery = 5,
                      sampleEvery = 5,
                      detectionSize = 1e5)
set.seed(1)
 b1 <- oncoSimulIndiv(p701, keepEvery = 5,
                      numGenes = 60,
                      sampleEvery = 5,
                      detectionSize = 1e8)
 b1

m2 <- oncoSimulIndiv(examplePosets[["p1101"]],
                     model = "McFL",
                     numGenes = 60,
                     detectionDrivers = 3,
                     mu = 5e-7,
                     initSize = 4000,
                     sampleEvery = 0.025,
                     finalTime = 25000,
                     keepEvery = 5,
                     detectionSize = 1e6) 
m2

' > i1.R

R-3.1.0 -q --vanilla < i1.R

