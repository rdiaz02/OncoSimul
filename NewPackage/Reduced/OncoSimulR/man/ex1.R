library(OncoSimulR)

data(example_trees)
set.seed(40) ## algo2
r1 <- oncoSimulIndiv(posetToAdjmat(p705))


## maxmemory does work in Linux
r2 <- oncoSimulIndiv(posetToAdjmat(p705), max.memory = 0.5)

set.seed(1)
r2 <- oncoSimulIndiv(posetToAdjmat(p705), max.memory = 0.7)

set.seed(1)
r2 <- oncoSimulIndiv(posetToAdjmat(p705), max.memory = 2)



posetToAdjmat(p705)




set.seed(40) ## algo2
r1 <- oncoSimulIndiv(posetToAdjmat(p705))

set.seed(35) ## ti_DBL_MIN
r1 <- oncoSimulIndiv(posetToAdjmat(p705))


set.seed(2)
r2 <- oncoSimulIndiv(posetToAdjmat(p705), "McF", sampleEvery = 0.05,
                     detectionDrivers = 3,
                     finalTime = 25*365)


set.seed(2)
r3 <- oncoSimulIndiv(posetToAdjmat(p705), "Exp",
                     detectionDrivers = 3,
                     finalTime = 0.25 * 25*365)

