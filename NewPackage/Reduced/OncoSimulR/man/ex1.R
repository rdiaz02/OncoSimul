### Looks like not needed
### Suggests: igraph, RBGL, Rgraphviz, Rtreemix, Oncotree


library(OncoSimulR)
data(examplePosets)
plotPoset(examplePosets[["p1101"]])
plotAdjMat(poset2AdjMat(examplePosets[["p1101"]]))

p705 <- examplePosets[["p705"]]

r1 <- oncoSimulIndiv(p705, detectionSize = 1e6)

p1 <- oncoSimulPop(20, p705, detectionSize = 1e6,
                   sampleEvery = 10, mc.cores = 4)



set.seed(40) ## algo2
r1 <- oncoSimulIndiv(posetToAdj(p705), detectionSize = 1e9)

m1 <- oncoSimulIndiv(posetToAdj(p705), model = "McF", detectionDrivers = 3, sampleEvery = 0.05)


set.seed(40) ## OK
r1 <- oncoSimulIndiv(posetToAdj(p705), detectionSize = 1e9, sampleEvery = 14)



set.seed(40) ## algo2

p1 <- oncoSimulPop(20, posetToAdj(p705), detectionSize = 1e9,
                   sampleEvery = 10, mc.cores = 4)


p2 <- oncoSimulPop(20, posetToAdj(p705), detectionSize = 1e9,
                   sampleEvery = 10, keep.every = 20, mc.cores = 4)


## do not restrict to cancer reaching
set.seed(1)
RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()
set.seed(1)
unix.time(
    p0 <- oncoSimulPop(5, posetToAdj(p705), onlyCancer = FALSE,
                   sampleEvery = 5, keepEvery = 20, mc.cores = 1)
    )

## ugly or break
plot(p0[[1]], thinPops = FALSE, log = "") ## unlog, to show it reaches 0
plot(p0[[2]], thinPops = FALSE)

set.seed(1)
RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()
set.seed(1)
unix.time(
    p1 <- oncoSimulPop(50, posetToAdj(p705), detectionSize = 1e6,
                   sampleEvery = 5, keepEvery = 20, mc.cores = 4)
    )

plot(p1[[1]], thinPops = FALSE)



set.seed(1)
RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()
set.seed(1)
unix.time(
    p2 <- oncoSimulPop(50, posetToAdj(p705), detectionSize = 1e6,
                   sampleEvery = 1, mc.cores = 4)
    )

plot(p2[[1]], thinPops = FALSE)





set.seed(1)
RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()
set.seed(1)
unix.time(
    p2 <- oncoSimulPop(500, posetToAdj(p705), detectionSize = 1e8,
                   sampleEvery = 5, keepEvery = 20, mc.cores = 4)
    )

Rprof()
p2 <- oncoSimulPop(500, posetToAdj(p705), detectionSize = 1e8,
                       sampleEvery = 5, keepEvery = 20, mc.cores = 4)
Rprof(NULL)


set.seed(1)
RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()
set.seed(1)
unix.time(
    p3 <- oncoSimulPop(50, posetToAdj(p705), detectionSize = 1e7,
                   sampleEvery = 10, keepEvery = 35, mc.cores = 4)
    )




posetToAdj(p705)




set.seed(40) ## algo2
r1 <- oncoSimulIndiv(posetToAdj(p705))

set.seed(35) ## ti_DBL_MIN
r1 <- oncoSimulIndiv(posetToAdj(p705))


set.seed(2)
r2 <- oncoSimulIndiv(posetToAdj(p705), "McF", sampleEvery = 0.05,
                     detectionDrivers = 3,
                     finalTime = 25*365)


set.seed(2)
r3 <- oncoSimulIndiv(posetToAdj(p705), "Exp",
                     detectionDrivers = 3,
                     finalTime = 0.25 * 25*365)




## ## maxmemory does work in Linux
## r2 <- oncoSimulIndiv(posetToAdj(p705), max.memory = 0.5)

## set.seed(1)
## r2 <- oncoSimulIndiv(posetToAdj(p705), max.memory = 0.7)

## set.seed(1)
## r2 <- oncoSimulIndiv(posetToAdj(p705), max.memory = 2)
