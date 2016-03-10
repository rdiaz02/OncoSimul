### This is code to make sure that, even if tests pass, there are no
### changes in the output of versions. No plotting here.

### Run with the version you want, and diff the output. seeds are fixed at
### values to make sure no changes. Of course, expect changes if there are
### changes in random number generation, such as if using randutils.


rm(list = ls())
library(OncoSimulR)
packageVersion("OncoSimulR")
library(help = OncoSimulR)


set.seed(1)
s1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = 0.1,
                 sh = c(rep(0, 4), c(-.9, -.9), rep(-.95, 3)),
                 typeDep = "SM")
(allFitnessEffects(s1))


set.seed(123)
tmp <- oncoSimulIndiv(allFitnessEffects(s1), model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("b, a")
                      )
summary(tmp)

data(examplesFitnessEffects)
set.seed(456)
tmp1 <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.015,
                       max.num.tries = 10,
                       keepEvery = 5,
                       initSize = 2000,
                       finalTime = 20000,
                       onlyCancer = FALSE,
                       extraTime = 1500,
                       keepPhylog = TRUE,
                       initMutant = "d > m")
summary(tmp1)


set.seed(789)
s1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = 0.1,
                 sh = c(rep(0, 4), c(-.9, -.9), rep(-.95, 3)),
                 typeDep = "SM")
(allFitnessEffects(s1))



set.seed(123)

tmp2 <- oncoSimulIndiv(allFitnessEffects(s1), model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000, keepPhylog = TRUE)
summary(tmp2)


for(i in 11:15){
    cat(paste("\n i = ", i))
    set.seed(i)
    tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                           model = "McFL", 
                           mu = 5e-5,
                           detectionSize = 1e8, 
                           detectionDrivers = 3,
                           sampleEvery = 0.015,
                           max.num.tries = 10,
                           keepEvery = 5,
                           initSize = 2000,
                           finalTime = 7000,
                           onlyCancer = FALSE,
                           extraTime = 1500,
                           keepPhylog = TRUE);
    print(summary(tmp))
}


set.seed(1010)
tmp6 <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.025,
                       max.num.tries = 10,
                       keepEvery = 5,
                       initSize = 2000,
                       finalTime = 3000,
                       onlyCancer = FALSE,
                       seed = 4);
summary(tmp6)


set.seed(99999)
tmp7 <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
                       model = "McFL", 
                       mu = 5e-5,
                       detectionSize = 1e8, 
                       detectionDrivers = 3,
                       sampleEvery = 0.025,
                       max.num.tries = 10,
                       keepEvery = 5,
                       initSize = 2000,
                       finalTime = 8000,
                       onlyCancer = FALSE,
                       seed = 4,
                       keepPhylog = TRUE);
summary(tmp7)


set.seed(987654)
nin <- 100
ne <- 10
s <- 0.1
sj <- -0.05
nn <- sapply(c("A", "B"), function(x) paste0(seq.int(ne), x))
int <- apply(nn, 1, function(x) paste(x, collapse = " : "))
single <- as.vector(nn)
epi <- c(rep(s, length(single)), rep(sj, length(int)))
names(epi) <- c(single, int)
ee <- allFitnessEffects(epistasis = epi,
                        noIntGenes = rexp(nin, 20))
see <- oncoSimulIndiv(ee, model = "Exp",
                      detectionSize = 1e7,
                      detectionDrivers = 20,
                      sampleEvery = 10,
                      keepEvery = -9)
summary(see)




set.seed(1234456)
data(examplePosets)
p705 <- examplePosets[["p705"]]
p1 <- oncoSimulIndiv(p705)
class(p1)
lp1 <- OncoSimulWide2Long(p1)
head(lp1)
summary(p1)
summary(lp1)

set.seed(33)
data(examplesFitnessEffects)
sm <-  oncoSimulIndiv(examplesFitnessEffects$cbn1,
                       model = "McFL", 
                       mu = 5e-7,
                       detectionSize = 1e8, 
                       detectionDrivers = 2,
                       sampleEvery = 0.025,
                       keepEvery = 5,
                       initSize = 2000,
                       onlyCancer = FALSE)
class(sm)
lsm <- OncoSimulWide2Long(sm)
head(lsm)
summary(lsm)
summary(sm)



set.seed(98123)
cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = 0.1,
                 sh = -0.9,
                 typeDep = "MN")
cbn1 <- allFitnessEffects(cs)
cbn1


p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                 child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                 s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                 sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                 typeDep = c(rep("--", 4), 
                     "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
oe <- c("C > F" = -0.1, "H > I" = 0.12)
sm <- c("I:J"  = -1)
sv <- c("-K:M" = -.5, "K:-M" = -.5)
epist <- c(sm, sv)
modules <- c("Root" = "Root", "A" = "a1",
             "B" = "b1, b2", "C" = "c1",
             "D" = "d1, d2", "E" = "e1",
             "F" = "f1, f2", "G" = "g1",
             "H" = "h1, h2", "I" = "i1",
             "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")
set.seed(1) ## for repeatability
noint <- rexp(5, 10)
names(noint) <- paste0("n", 1:5)
fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                         noIntGenes = noint, geneToModule = modules)
fea




# A three-gene epistasis example
sa <- 0.1
sb <- 0.15
sc <- 0.2
sab <- 0.3
sbc <- -0.25
sabc <- 0.4
sac <- (1 + sa) * (1 + sc) - 1
E3A <- allFitnessEffects(epistasis =
                            c("A:-B:-C" = sa,
                              "-A:B:-C" = sb,
                              "-A:-B:C" = sc,
                              "A:B:-C" = sab,
                              "-A:B:C" = sbc,
                              "A:-B:C" = sac,
                              "A : B : C" = sabc)
                                                )
evalAllGenotypes(E3A, order = FALSE, addwt = FALSE)
evalAllGenotypes(E3A, order = FALSE, addwt = TRUE,  model = "Bozic")
evalGenotype("B, C", E3A, verbose = TRUE)


## Order effects and modules
ofe2 <- allFitnessEffects(orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
evalAllGenotypes(ofe2, max = 325)[1:15, ]
## Next two are identical
evalGenotype("d1 > d2 > f3", ofe2, verbose = TRUE)
evalGenotype("d1 , d2 , f3", ofe2, verbose = TRUE)
## This is different
evalGenotype("f3 , d1 , d2", ofe2, verbose = TRUE)
## but identical to this one
evalGenotype("f3 > d1 > d2", ofe2, verbose = TRUE)



## Restrictions in mutations as a graph. Modules present.
p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                  child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                  s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                  sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                  typeDep = c(rep("--", 4), 
                      "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
fp4m <- allFitnessEffects(p4,
                          geneToModule = c("Root" = "Root", "A" = "a1",
                              "B" = "b1, b2", "C" = "c1",
                              "D" = "d1, d2", "E" = "e1",
                              "F" = "f1, f2", "G" = "g1"))
evalAllGenotypes(fp4m, order = FALSE, max = 1024, addwt = TRUE)[1:15, ]
evalGenotype("b1, b2, e1, f2, a1", fp4m, verbose = TRUE)
## Of course, this is identical; b1 and b2 are same module
## and order is not present here
evalGenotype("a1, b2, e1, f2", fp4m, verbose = TRUE)
evalGenotype("a1 > b2 > e1 > f2", fp4m, verbose = TRUE)
## We can use the exact same integer numeric id codes as in the
##   fitnessEffects geneModule component:
evalGenotype(c(1L, 3L, 7L, 9L), fp4m, verbose = TRUE)




set.seed(99999)
data(examplePosets)
p701 <- examplePosets[["p701"]]
## Bozic Model
b1 <- oncoSimulIndiv(p701)
summary(b1)

## McFarland; use a small sampleEvery, but also a reasonable
##   keepEvery.
## We also modify mutation rate to values similar to those in the
##   original paper.
## Note that detectionSize will play no role
## finalTime is large, since this is a slower process
## initSize is set to 4000 so the default K is larger and we are likely
## to reach cancer. Alternatively, set K = 2000.
set.seed(999)
m1 <- oncoSimulIndiv(p701,
                     model = "McFL",
                     mu = 5e-7,
                     initSize = 4000,
                     sampleEvery = 0.025,
                     finalTime = 15000,
                     keepEvery = 10,
                     onlyCancer = FALSE)
summary(m1)

## Simulating 4 individual trajectories
## (I set mc.cores = 2 to comply with --as-cran checks, but you
##  should either use a reasonable number for your hardware or
##  leave it at its default value).

set.seed(8888)
p1 <- oncoSimulPop(4, p701,
                   keepEvery = 10,
                   mc.cores = 2)
summary(p1)

set.seed(777)
p2 <- oncoSimulSample(4, p701)
p2




#### A model similar to the one in McFarland. We use 2070 genes.
set.seed(45612321)
nd <- 70  
np <- 2000 
s <- 0.1  
sp <- 1e-3 
spp <- -sp/(1 + sp)
mcf1 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                          drv = seq.int(nd))
mcf1s <-  oncoSimulIndiv(mcf1,
                         model = "McFL", 
                         mu = 1e-7,
                         detectionSize = 1e8, 
                         detectionDrivers = 100,
                         sampleEvery = 0.02,
                         keepEvery = 2,
                         initSize = 2000,
                         finalTime = 1000,
                         onlyCancer = FALSE)
summary(mcf1s)


#### Order effects with modules, and 5 genes without interactions
#### with fitness effects from an exponential distribution
set.seed(123214)
oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
oiI1 <- oncoSimulIndiv(oi, model = "Exp")
oiI1$GenotypesLabels
oiI1   ## note the order and separation by "_"
oiP1 <- oncoSimulPop(2, oi,
                     keepEvery = 10,
                     mc.cores = 2)
summary(oiP1) 
## Even if order exists, this cannot reflect it;
## G1 to G10 are d1, d2, f1..,f3, and the 5 genes without
## interaction
samplePop(oiP1)
oiS1 <- oncoSimulSample(2, oi)
## The output contains only the summary of the runs AND
## the sample:
oiS1 
## And their sizes do differ
object.size(oiS1)
object.size(oiP1)



######## Using a poset for pancreatic cancer from Gerstung et al.
###      (s and sh are made up for the example; only the structure
###       and names come from Gerstung et al.)

set.seed(876678)
pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                          "TP53", "TP53", "MLL3"),
                                      child = c("KRAS","SMAD4", "CDNK2A", 
                                          "TP53", "MLL3",
                                          rep("PXDN", 3), rep("TGFBR2", 2)),
                                      s = 0.05,
                                      sh = -0.3,
                                      typeDep = "MN"))
### Use an exponential growth model
pancr1 <- oncoSimulIndiv(pancr, model = "Exp")
pancr1
summary(pancr1)
pancr1$GenotypesLabels


## Pop and Sample
set.seed(89862)
pancrPop <- oncoSimulPop(4, pancr,
                        keepEvery = 10,
                       mc.cores = 2)
summary(pancrPop)
set.seed(89862)
pancrSPop <- samplePop(pancrPop)
pancrSPop
set.seed(89862)
pancrSamp <- oncoSimulSample(2, pancr)
pancrSamp



set.seed(341208)
sa <- 0.1
sb <- -0.2
sab <- 0.25
sac <- -0.1
sbc <- 0.25
sv2 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                                       "A : -B" = sa,
                                       "A : C" = sac,
                                       "A:B" = sab,
                                       "-A:B:C" = sbc),
                         geneToModule = c(
                             "Root" = "Root",
                             "A" = "a1, a2",
                             "B" = "b",
                             "C" = "c"))
evalAllGenotypes(sv2, order = FALSE, addwt = TRUE)
e1 <- oncoSimulIndiv(sv2, model = "McFL",
                     mu = 5e-6,
                     sampleEvery = 0.02,
                     keepEvery = 1,
                     initSize = 2000,
                     finalTime = 3000,
                     onlyCancer = FALSE)
summary(e1)
e1


set.seed(11)
data(examplePosets)
p705 <- examplePosets[["p705"]]
## (I set mc.cores = 2 to comply with --as-cran checks, but you
##  should either use a reasonable number for your hardware or
##  leave it at its default value).
p1 <- oncoSimulPop(4, p705, mc.cores = 2)
samplePop(p1)


## Now single cell sampling
set.seed(12)
r1 <- oncoSimulIndiv(p705)
samplePop(r1, typeSample = "single")

set.seed(21)
r1 <- oncoSimulIndiv(p705)
samplePop(r1, typeSample = "single")
