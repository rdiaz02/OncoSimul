## Do:
##  - finish testing formally
##       - equivalence of tree transformations
##       - equivalence of keepEvery = -9
##  - launch simuls with only single sample
##  - work on keepTheseMany








library(OncoSimulR)

data(examplePosets)
p701 <- examplePosets[["p701"]]

set.seed(1); b1 <- oncoSimulIndiv(p701)
set.seed(1); b2 <- oncoSimulIndiv(p701, keepEvery = -9, verbosity = 1)

summary(b1)
summary(b2)


set.seed(23)
m1 <- oncoSimulIndiv(p701,
                     model = "McFL",
                     mu = 5e-7,
                     initSize = 4000,
                     sampleEvery = 0.025,
                     finalTime = 15000,
                     keepEvery = 5)
set.seed(23)
m2 <- oncoSimulIndiv(p701,
                     model = "McFL",
                     mu = 5e-7,
                     initSize = 4000,
                     sampleEvery = 0.025,
                     finalTime = 15000,
                     keepEvery = -9)
summary(m1)
summary(m2)




set.seed(2)
m1 <- oncoSimulIndiv(p701,
                     model = "McFL",
                     mu = 5e-7,
                     initSize = 4000,
                     sampleEvery = 0.025,
                     finalTime = 15000,
                     keepEvery = 5)
set.seed(2)
m2 <- oncoSimulIndiv(p701,
                     model = "McFL",
                     mu = 5e-7,
                     initSize = 4000,
                     sampleEvery = 0.025,
                     finalTime = 15000,
                     keepEvery = -9)
summary(m1)
summary(m2)




cucu <- OncoSimulR:::oncoSimulSample(5, p701, mc.cores = 1)

coco <- OncoSimulR:::oncoSimulSample(5, p701, mc.cores = 1, verbosity = 0)





m2 <- oncoSimulIndiv(p701,
                     model = "McFL",
                     mu = 5e-7,
                     initSize = 4000,
                     sampleEvery = 0.025,
                     endTimeEvery = 5 * 0.025,
                     finalTime = 15000,
                     keepEvery = -9)
summary(m2)
m2$pops.by.time





set.seed(1); b1 <- oncoSimulIndiv(p701)

set.seed(2); b2 <- oncoSimulIndiv(p701)
set.seed(3); b3 <- oncoSimulIndiv(p701)
set.seed(4); b4 <- oncoSimulIndiv(p701)
set.seed(5); b5 <- oncoSimulIndiv(p701)
set.seed(6); b6 <- oncoSimulIndiv(p701)



ff <- function(seed) {
    cat("\n Bozic \n")

    set.seed(seed); b <- oncoSimulIndiv(p701)
    set.seed(seed); bb <- oncoSimulIndiv(p701, keepEvery = -9)

    print(summary(b))
    print(summary(bb))
    print(b$pops.by.time[nrow(b$pops.by.time), ])
    print(bb$pops.by.time)

    cat("\n mcfl \n")
    set.seed(seed); m <- oncoSimulIndiv(p701,
                                        model = "McFL",
                                        mu = 5e-7,
                                        initSize = 4000,
                                        sampleEvery = 0.025,
                                        finalTime = 15000,
                                        keepEvery = 5)
    
    set.seed(seed); mm <- oncoSimulIndiv(p701,
                                         model = "McFL",
                                         mu = 5e-7,
                                         initSize = 4000,
                                         sampleEvery = 0.025,
                                         finalTime = 15000,
                                         keepEvery = -9)

    print(summary(m))
    print(summary(mm))
    print(m$pops.by.time[nrow(m$pops.by.time), ])
    print(mm$pops.by.time)
    return(list(b, bb, m, mm))
   
}






