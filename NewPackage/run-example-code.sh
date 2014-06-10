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
                     sampleEvery = 5,
                     detectionSize = 1e7)
' > i1.R

R-3.1.0 -q --vanilla < i1.R

