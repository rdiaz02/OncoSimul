### This is code to make sure that, even if tests pass, there are no
### changes in the output of versions. No plotting here.

### Run with the version you want, and diff the output. seeds are fixed at
### values to make sure no changes. Of course, expect changes if there are
### changes in random number generation, such as if using randutils.


rm(list = ls())
## for reproducibility with mclapply
date()
RNGkind("L'Ecuyer-CMRG")
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
set.seed(123)
oiP1 <- oncoSimulPop(2, oi,
                     keepEvery = 10,
                     mc.cores = 2,
                     seed = NULL)
summary(oiP1) 
## Even if order exists, this cannot reflect it;
## G1 to G10 are d1, d2, f1..,f3, and the 5 genes without
## interaction
samplePop(oiP1)

set.seed(33)
oiS1 <- oncoSimulSample(2, oi, seed = NULL)
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
                         mc.cores = 2, seed = NULL)
summary(pancrPop)
set.seed(89862)
pancrSPop <- samplePop(pancrPop)
pancrSPop
set.seed(89862)
pancrSamp <- oncoSimulSample(2, pancr, seed = NULL)
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
p1 <- oncoSimulPop(4, p705, mc.cores = 2, seed = NULL)
samplePop(p1)


## Now single cell sampling
set.seed(12)
r1 <- oncoSimulIndiv(p705)
samplePop(r1, typeSample = "single")

set.seed(21)
r1 <- oncoSimulIndiv(p705)
samplePop(r1, typeSample = "single")








############################################################

sd <- 0.1
sdp <- 0.15
sp <- 0.05
bauer <- data.frame(parent = c("Root", rep("p", 5)),
                    child = c("p", paste0("s", 1:5)),
                    s = c(sd, rep(sdp, 5)),
                    sh = c(0, rep(sp, 5)),
                    typeDep = "MN")
b1 <- evalAllGenotypes(allFitnessEffects(bauer), order = FALSE)
b2 <- evalAllGenotypes(allFitnessEffects(bauer), order = TRUE, max = 2000)
b1
b2


sd <- 0.1
sdp <- 0.15
sp <- 0.05
bauer <- data.frame(parent = c("Root", rep("p", 5)),
                    child = c("p", paste0("s", 1:5)),
                    s = c(sd, rep(sdp, 5)),
                    sh = c(0, rep(sp, 5)),
                    typeDep = "MN")
b1 <- evalAllGenotypes(allFitnessEffects(bauer), order = FALSE)
b2 <- evalAllGenotypes(allFitnessEffects(bauer), order = TRUE, max = 2000)
b1
b2


sd <- 0.1
sdp <- 0.15
sp <- 0.05
bauer <- data.frame(parent = c("Root", rep("p", 5)),
                    child = c("p", paste0("s", 1:5)),
                    s = c(sd, rep(sdp, 5)),
                    sh = c(0, rep(sp, 5)),
                    typeDep = "MN")
b1 <- evalAllGenotypes(allFitnessEffects(bauer), order = FALSE)
bauer3 <- data.frame(parent = c("Root", rep("u", 5)),
                     child = c("u", paste0("s", 1:5)),
                     s = c(sd, rep(sdp, 5)),
                     sh = c(0, rep(sp, 5)),
                     typeDep = "MN")
b3 <- evalAllGenotypes(allFitnessEffects(bauer), order = TRUE, max = 2000)
b1
b3


sd <- 0.1
sdp <- 0.15
sp <- 0.05
bauer <- data.frame(parent = c("Root", rep("p", 5)),
                    child = c("p", paste0("s", 1:5)),
                    s = c(sd, rep(sdp, 5)),
                    sh = c(0, rep(sp, 5)),
                    typeDep = "MN")
b1 <- evalAllGenotypes(allFitnessEffects(bauer), order = FALSE)
bauer3 <- data.frame(parent = c(rep("u", 5), "Root"),
                     child = c(paste0("s", 1:5), "u"),
                     s = c(sd, rep(sdp, 5)),
                     sh = c(0, rep(sp, 5)),
                     typeDep = "MN")
b3 <- evalAllGenotypes(allFitnessEffects(bauer), order = TRUE, max = 2000)
b1
b3


o1 <- evalAllGenotypes(allFitnessEffects(
    orderEffects = c("d>f" = 0.4, "f > d" = -0.3) ))
o2 <- evalAllGenotypes(allFitnessEffects(
    orderEffects = c("f > d" = -0.3, "d > f" = 0.4) ))
o1
o2


ofe1 <- allFitnessEffects(orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2",
                                "D" = "d1, d2") )
ag <- evalAllGenotypes(ofe1)
ag


ofe2 <- allFitnessEffects(orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
ag2 <- evalAllGenotypes(ofe2, max = 326)
ag2


o1 <- evalAllGenotypes(allFitnessEffects(
    orderEffects = c("F > D" = -0.3, "D > F" = 0.4),  
    geneToModule = c("Root" = "Root", "F" = "d", "D" = "f")))
o2 <- evalAllGenotypes(allFitnessEffects(
    orderEffects = c("D>F" = 0.4, "F >D" = -0.3),  
    geneToModule = c("Root" = "Root", "F" = "d", "D" = "f")))
o1
o2

o3 <- allFitnessEffects(orderEffects = c(
                            "F > D > M" = -0.3,
                            "D > F > M" = 0.4,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.5
                        ),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d") )
ag <- evalAllGenotypes(o3)
ag


ai1 <- evalAllGenotypes(allFitnessEffects(
    noIntGenes = c(0.05, -.2, .1)), order = FALSE)

ai2 <- evalAllGenotypes(allFitnessEffects(
    noIntGenes = c(0.05, -.2, .1)), order = TRUE)
ai1
ai2


ai3 <- evalAllGenotypes(allFitnessEffects(
    noIntGenes = c("a" = 0.05, "b" = -.2, "c" = .1)), order = FALSE)

ai4 <- evalAllGenotypes(allFitnessEffects(
    noIntGenes = c("a" = 0.05, "b" = -.2, "c" = .1)), order = TRUE)
ai3
ai4


ai3 <- evalAllGenotypes(allFitnessEffects(
    noIntGenes = c("m" = 0.05, "b" = -.2, "f" = .1)), order = FALSE)

ai4 <- evalAllGenotypes(allFitnessEffects(
    noIntGenes = c("m" = 0.05, "b" = -.2, "f" = .1)), order = TRUE)
ai3
ai4

foi1 <- allFitnessEffects(
    orderEffects = c("D>B" = -0.3, "B > D" = 0.3),
    noIntGenes = c("A" = 0.05, "C" = -.2, "E" = .1))
agoi1 <- evalAllGenotypes(foi1,  max = 325)
agoi1


foi1 <- allFitnessEffects(
    orderEffects = c("D>B" = -0.3, "B > D" = 0.3),
    noIntGenes = c("M" = 0.05, "A" = -.2, "J" = .1))
agoi1 <- evalAllGenotypes(foi1,  max = 325)
agoi1

s <- 0.2
sv <- allFitnessEffects(epistasis = c("-A : B" = -1,
                                      "A : -B" = -1,
                                      "A:B" = s))
evalAllGenotypes(sv, order = TRUE)
evalAllGenotypes(sv, order = FALSE)


sa <- -0.1
sb <- -0.2
sab <- 0.25
sv2 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                                       "A : -B" = sa,
                                       "A:B" = sab),
                         geneToModule = c(
                             "Root" = "Root",
                             "A" = "a1, a2",
                             "B" = "b"))
evalAllGenotypes(sv2, order =TRUE)
evalAllGenotypes(sv2, order = FALSE)

sa <- 0.1
sb <- 0.2
sab <- -0.8
sm1 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                                       "A : -B" = sa,
                                       "A:B" = sab))
evalAllGenotypes(sm1, order = TRUE)
evalAllGenotypes(sm1, order = FALSE)


sa <- 0.2
sb <- 0.3
sab <- 0.7
e2 <- allFitnessEffects(epistasis =
                            c("A: -B" = sa,
                              "-A:B" = sb,
                              "A : B" = sab))
evalAllGenotypes(e2, order = TRUE)
evalAllGenotypes(e2, order = FALSE)

sa <- 0.2
sb <- 0.3
sab <- 0.7
s2 <- ((1 + sab)/((1 + sa) * (1 + sb))) - 1
e3 <- allFitnessEffects(epistasis =
                            c("A" = sa,
                              "B" = sb,
                              "A : B" = s2))
evalAllGenotypes(e3, order = TRUE)
evalAllGenotypes(e3, order = FALSE)



sa <- 0.1
sb <- 0.15
sc <- 0.2
sab <- 0.3
sbc <- -0.25
sabc <- 0.4
sac <- (1 + sa) * (1 + sc) - 1
E3 <- allFitnessEffects(epistasis =
                            c("A:-B:-C" = sa,
                              "-A:B:-C" = sb,
                              "-A:-B:C" = sc,
                              "A:B:-C" = sab,
                              "-A:B:C" = sbc,
                              "A:-B:C" = sac,
                              "A : B : C" = sabc)
                        )
evalAllGenotypes(E3, order = TRUE)
evalAllGenotypes(E3, order = FALSE)



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
Sab <- ( (1 + sab)/((1 + sa) * (1 + sb))) - 1
Sbc <- ( (1 + sbc)/((1 + sb) * (1 + sc))) - 1
Sabc <- ( (1 + sabc)/( (1 + sa) * (1 + sb) * (1 + sc) * (1 + Sab) * (1 + Sbc) ) ) - 1
E3B <- allFitnessEffects(epistasis =
                             c("A" = sa,
                               "B" = sb,
                               "C" = sc,
                               "A:B" = Sab,
                               "B:C" = Sbc,
                               ## "A:C" = sac,
                               "A : B : C" = Sabc)
                         )

evalAllGenotypes(E3A, order = TRUE)
evalAllGenotypes(E3A, order = FALSE)

evalAllGenotypes(E3B, order = TRUE)
evalAllGenotypes(E3B, order = FALSE)

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
ge3a <- evalAllGenotypes(E3A, order = FALSE, addwt = FALSE)
ge3ao <- evalAllGenotypes(E3A, order = TRUE, addwt = FALSE)
Sab <- ( (1 + sab)/((1 + sa) * (1 + sb))) - 1
Sbc <- ( (1 + sbc)/((1 + sb) * (1 + sc))) - 1
Sabc <- ( (1 + sabc)/( (1 + sa) * (1 + sb) * (1 + sc) * (1 + Sab) * (1 + Sbc) ) ) - 1
E3B <- allFitnessEffects(epistasis =
                             c("A" = sa,
                               "B" = sb,
                               "C" = sc,
                               "A:B" = Sab,
                               "B:C" = Sbc,
                               ## "A:C" = sac,
                               "A : B : C" = Sabc)
                         )
ge3b <- evalAllGenotypes(E3A, order = FALSE, addwt = FALSE)
ge3bo <- evalAllGenotypes(E3A, order = TRUE, addwt = FALSE)

ge3a
ge3ao
ge3b
ge3bo


c1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                 sh = c(rep(0, 4), c(-.1, -.1), rep(-.05, 3)),
                 typeDep = "MN")
fc1 <- allFitnessEffects(c1)
(gfc1 <- evalAllGenotypes(fc1, order = FALSE))
(gfc1o <- evalAllGenotypes(fc1, order = TRUE, max = 1956))


c1b <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                  child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                  s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                  sh = c(rep(0, 4), c(-.1, -.1), rep(-.05, 3)),
                  typeDep = c("-", "--", "SM", "XMPN", rep("MN",5)))
fc1b <- allFitnessEffects(c1b)
(gfc1b <- evalAllGenotypes(fc1b, order = FALSE))


s1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                 sh = c(rep(0, 4), c(-.1, -.1), rep(-.05, 3)),
                 typeDep = "SM")
fs1 <- allFitnessEffects(s1)
(gfs1 <- evalAllGenotypes(fs1, order = FALSE))
(gfs1o <- evalAllGenotypes(fs1, order = TRUE, max = 1956))

zzz <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                  child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                  s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                  sh = c(rep(0, 4), c(-.1, -.1), rep(-.05, 3)),
                  typeDep = c("-", "--", "SM", "XMPN", rep("SM", 5)))
zzz <- allFitnessEffects(zzz)
(zzz <- evalAllGenotypes(zzz, order = FALSE))



x1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                 sh = c(rep(0, 4), c(-.9, -.9), rep(-.95, 3)),
                 typeDep = "XMPN")
fx1 <- allFitnessEffects(x1)
(gfx1 <- evalAllGenotypes(fx1, order = FALSE))
(gfx1o <- evalAllGenotypes(fx1, order = TRUE, max = 1956))

zzz <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                  child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                  s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                  sh = c(rep(0, 4), c(-.9, -.9), rep(-.95, 3)),
                  typeDep = c("SM", "-", "--", "XMPN", rep("XMPN", 5))) 
zzz <- allFitnessEffects(zzz)
(zzz <- evalAllGenotypes(zzz, order = FALSE))

p3 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c", "f"),
                 child = c("a", "b", "d", "e", "c", "c", "f", "f", "g", "g"),
                 s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                 sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                 typeDep = c(rep("--", 4), 
                             "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
fp3 <- allFitnessEffects(p3)
(gfp3 <- evalAllGenotypes(fp3, order = FALSE))

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
(gfp4 <- evalAllGenotypes(fp4m, order = FALSE, max = 1024))



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
set.seed(1111) ## for repeatability
## These are seeds in R; no problems with different compilers, etc.
noint <- rexp(5, 10)
names(noint) <- paste0("n", 1:5)
(fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                          noIntGenes = noint, geneToModule = modules))





###############################################


set.seed(53)
              oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.3, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2") )
              out <- oncoSimulPop(4,
                                  oi, 
                                  detectionSize = 1e4,
                                  onlyCancer = FALSE, seed = NULL)
             out







######################################################################

set.seed(26)
            o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c(0.01, 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              tmp <- oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("d > m")
                      )
tmp


set.seed(27)
              o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c("u" = 0.01, "z" = 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              tmp <- oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("d > m > z")
                      )


set.seed(97)
              o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c("u" = 0.01, "z" = 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              tmp <- oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("m > d")
                      )

tmp



set.seed(98)
              o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c("u" = 0.01, "z" = 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              tmp <- oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("m > u > d")
                      )

tmp


set.seed(103)
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9),
                        noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = 0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d") )
    ossI <- oncoSimulSample(4, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 2,
                        onlyCancer = TRUE,
                        initSize = 500,
                        initMutant = c("z > d"),
                        seed = NULL,
                        thresholdWhole = 1 ## check presence of initMutant
                        )

ossI

set.seed(127)



    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9,
                            "M > A"     = 0.25,
                            "A > H"     = 0.2,
                            "A > G"     = 0.3),
                        noIntGenes = c("u" = 0.1, 
                                       "v" = 0.2,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "A" = "a",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d",
                              "H" = "h",
                              "G" = "g") )
    ossI <- oncoSimulSample(4, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 3,
                        onlyCancer = TRUE,
                        initSize = 500,
                        initMutant = c("z > a"),
                        seed = NULL,
                        thresholdWhole = 1 ## check presence of initMutant
                        )
ossI


set.seed(987)
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.2,
                            "M > D"     = 0.9),
                        noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d") )
    ospI <- oncoSimulPop(4, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 3,
                        onlyCancer = TRUE,
                        keepPhylog = TRUE,
                        initSize = 500,
                        seed = NULL,
                        initMutant = c("d > m > y"),
                        mc.cores = 2
                        )
ospI

set.seed(56)

    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9),
                            noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d") )
    ospI <- oncoSimulPop(4, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 70,
                        detectionDrivers = 4, ## yes, reach end
                        onlyCancer = FALSE,
                        keepPhylog = TRUE,
                        initSize = 100,
                        initMutant = c("m > v > d"),
                        mc.cores = 2,
                        seed = NULL
                        )
ospI
## Up to, and including, test.init-mutant.R

set.seed(8)
oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.3, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2") )
set.seed(9)
oncoSimulIndiv(oi, 
               detectionSize = 1e4,
               onlyCancer = FALSE)

set.seed(9)
oncoSimulPop(4, oi, 
             detectionSize = 1e4,
             seed = NULL,
               onlyCancer = FALSE)


set.seed(10)
pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                                        "TP53", "TP53", "MLL3"),
                                                    child = c("KRAS","SMAD4", "CDNK2A", 
                                                        "TP53", "MLL3",
                                                        rep("PXDN", 3), rep("TGFBR2", 2)),
                                                    s = 0.05,
                                                    sh = -0.3,
                                      typeDep = "MN"))
oncoSimulSample(2, pancr, seed = NULL)
set.seed(11)
oncoSimulSample(2,
                pancr,
                seed = NULL,
                typeSample = "single")

set.seed(25)
    np <- 20
    s <- 0.015
    spp <- 0.01
    nd <- 5
    mcf1 <- allFitnessEffects(noIntGenes = rep(spp, np),
                              drvNames = integer(0))
    mcf2 <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                              drvNames = character(0))
set.seed(26)
oncoSimulSample(2, mcf1, seed = NULL)
set.seed(27)
oncoSimulSample(5, mcf2, seed = NULL)

## up to test.oncoSimulIndiv-miscell







data(examplesFitnessEffects)

## sometimes cancer is not reached. No problem.

## Very rarely, popSize > 1e15, and we get an exception. Decrease
## sampleEvery. And e2 only has two genes.

for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    if (names(examplesFitnessEffects)[16] == "e2") {
        detectionDrv <- 2
        sE <- 0.05
    } else {
        detectionDrv <- 4
        sE <- 2
    }
    set.seed(16)
    tmp <-  oncoSimulIndiv(examplesFitnessEffects[[i]],
                           model = "Bozic", 
                           mu = 1e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = detectionDrv,
                           sampleEvery = sE,
                           max.num.tries = 100,
                           initSize = 2000,
                           onlyCancer = FALSE)
    print(tmp)
}

for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    if (names(examplesFitnessEffects)[16] == "e2") {
        detectionDrv <- 2
        sE <- 0.05
    } else {
        detectionDrv <- 4
        sE <- 2
    }
    set.seed(17)
    tmp <-  oncoSimulIndiv(examplesFitnessEffects[[i]],
                           model = "Exp", 
                           mu = 1e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = detectionDrv,
                           sampleEvery = sE,
                           max.num.tries = 100,
                           initSize = 2000,
                           onlyCancer = FALSE)
    print(tmp)
}


for(i in 1:length(examplesFitnessEffects)) {
        set.seed(18)
    tmp <-  oncoSimulIndiv(examplesFitnessEffects[[i]],
                           model = "McFL", 
                           mu = 5e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = 2,
                           sampleEvery = 0.025,
                           max.num.tries = 10,
                           initSize = 2000,
                           finalTime = 15000,
                           onlyCancer = FALSE)
    print(tmp)
}



for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
        cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    if (names(examplesFitnessEffects)[16] == "e2") {
        sE <- 0.05
    } else {
        sE <- 1
    }
        set.seed(26)
    tmp <-  oncoSimulSample(4, examplesFitnessEffects[[i]],
                            onlyCancer = FALSE, seed = NULL,
                            sampleEvery = sE)
    print(tmp)
}



for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    if (names(examplesFitnessEffects)[16] == "e2") {
        detectionDrv <- 2
        sE <- 0.05
    } else {
        detectionDrv <- 4
        sE <- 2
    }
        set.seed(416)
    tmp <-  oncoSimulPop(4, examplesFitnessEffects[[i]],
                         onlyCancer = FALSE,
                         detectionDrivers = detectionDrv,
                         sampleEvery = sE, seed = NULL, 
                         mc.cores = 2)
    print(tmp)
    set.seed(419)
    tmp2 <- samplePop(tmp)
    print(tmp2)
}








