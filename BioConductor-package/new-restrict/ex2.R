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
noint <- rexp(25, 10)
names(noint) <- paste0("n", 1:25)
fea2 <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                         noIntGenes = noint, geneToModule = modules)
set.seed(1)
fs <- oncoSimulIndiv(fea2, model = "Exp")


unix.time({
set.seed(1)
fea2_s1 <- oncoSimulPop(50, fea2, model = "Exp",
                        keepEvery = 10,
                        mc.cores = 1)
set.seed(2)
fea2_s12 <- oncoSimulPop(10, fea2, model = "McFL",
                        keepEvery = 10,
                        mc.cores = 1)
})
## 60.03. con auto const: 62

unix.time({
set.seed(1)
fea2_s1x <- oncoSimulPop(50, fea2, model = "Exp",
                        keepEvery = 10,
                        mc.cores = 2)
set.seed(2)
fea2_s12x <- oncoSimulPop(10, fea2, model = "McFL",
                        keepEvery = 10,
                        mc.cores = 2)
})
## 89 and 55. con auto const: 84 and 45

unix.time({
set.seed(3)
fea2_s13 <- oncoSimulPop(50, fea2, model = "Exp",
                        keepEvery = 10,
                        mc.cores = 4)
set.seed(4)
fea2_s13x <- oncoSimulPop(10, fea2, model = "McFL",
                        keepEvery = 10,
                        mc.cores = 4)
})
## 118 and 48. con auto const: 108 and 38 or 230 and 33 or 81 and 31

unix.time({
set.seed(99)
fea2_s15 <- oncoSimulPop(500, fea2, model = "Exp",
                        keepEvery = 10,
                        mc.cores = 1)
set.seed(100)
fea2_s15x <- oncoSimulPop(100, fea2, model = "McFL",
                        keepEvery = 10,
                        mc.cores = 1)
})
## auto const: 690 and 691

unix.time({
set.seed(999)
fea2_s15_2 <- oncoSimulPop(200, fea2, model = "Exp",
                        keepEvery = 10,
                        mc.cores = 2)
set.seed(1000)
fea2_s15_2x <- oncoSimulPop(50, fea2, model = "McFL",
                        keepEvery = 10,
                        mc.cores = 2)
})
## auto const 282 and 216



load("examplesFitnessEffects.RData")

for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    tmp <-  oncoSimulIndiv(fE = examplesFitnessEffects[[i]],
                           model = "Bozic", 
                           mu = 1e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = 4,
                           sampleEvery = 2,
                           max.num.tries = 100,
                           initSize = 2000,
                           onlyCancer = TRUE)
}

for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    tmp <-  oncoSimulIndiv(fE = examplesFitnessEffects[[i]],
                           model = "Exp", 
                           mu = 1e-6,
                           detectionSize = 1e8, 
                           detectionDrivers = 4,
                           sampleEvery = 2,
                           max.num.tries = 100,
                           initSize = 2000,
                           onlyCancer = TRUE)
}


for(i in 1:length(examplesFitnessEffects)) {
    cat(paste("\n Doing i = ", i , " name = ",
              names(examplesFitnessEffects)[i], "\n"))
    tmp <-  oncoSimulIndiv(fE = examplesFitnessEffects[[i]],
                           model = "McFL", 
                           mu = 5e-7,
                           detectionSize = 1e8, 
                           detectionDrivers = 2,
                           sampleEvery = 0.025,
                           max.num.tries = 10,
                           initSize = 2000,
                           finalTime = 15000,
                           onlyCancer = TRUE)
}




