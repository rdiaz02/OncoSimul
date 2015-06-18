p30 <- cbind(0, 30)
plotPoset(p30, addroot = TRUE)

set.seed(1)
m1 <- oncoSimulIndiv(p30, model = "McFL", numPassengers = 0, mu = 1e-6,
                     detectionDrivers = 6, keepEvery = 5,
                     sampleEvery = 0.025)
plot(m1, addtot = TRUE)


## the following clearly shows clonal interference. There are several
## clones with 3 drivers, and with 4 drivers and several with 5.
set.seed(2)
m2 <- oncoSimulIndiv(p30, model = "McFL", numPassengers = 0, mu = 1e-5,
                     detectionDrivers = 8, keepEvery = 2,
                     sampleEvery = 0.025, initSize = 500)
plot(m2, addtot = FALSE, plotClones = TRUE, plotDrivers = TRUE,
     lwdClone = 0.9)



## even more in here
set.seed(3)
m3 <- oncoSimulIndiv(p30, model = "McFL", numPassengers = 0, mu = 1e-5,
                     detectionDrivers = 8, keepEvery = 2,
                     sampleEvery = 0.025, initSize = 1000)
plot(m3, addtot = FALSE, plotClones = TRUE, plotDrivers = TRUE,
     lwdClone = 0.9)


## a little in here
set.seed(4)
m4 <- oncoSimulIndiv(p30, model = "McFL", numPassengers = 0, mu = 1e-6,
                     detectionDrivers = 8, keepEvery = 2,
                     sampleEvery = 0.025, initSize = 1000)
plot(m4, addtot = FALSE, plotClones = TRUE, plotDrivers = TRUE,
     lwdClone = 0.9)


set.seed(4)
m5 <- oncoSimulIndiv(p30, model = "McFL", numPassengers = 0, mu = 1e-6,
                     detectionDrivers = 8, keepEvery = 2,
                     sampleEvery = 0.025, initSize = 4000)
plot(m5, addtot = FALSE, plotClones = TRUE, plotDrivers = TRUE,
     lwdClone = 0.9)

plot(m5, addtot = FALSE, plotClones = TRUE, plotDrivers = TRUE,
     lwdClone = 0.9, log = "")



## more here
set.seed(4)
m5 <- oncoSimulIndiv(p30, model = "McFL", numPassengers = 0, mu = 5e-6,
                     detectionDrivers = 8, keepEvery = 2,
                     sampleEvery = 0.025, initSize = 1000)

plot(m5, addtot = FALSE, plotClones = TRUE, plotDrivers = TRUE,
     lwdClone = 0.9)


domcf <- function(mu, s, initSize, keepEvery = 2) {
    m <- oncoSimulIndiv(p30, model = "McFL", numPassengers = 0, mu = mu,
                        s = s,
                     detectionDrivers = 8, keepEvery = keepEvery,
                        sampleEvery = 0.025, initSize = initSize)
    plot(m, addtot = TRUE, plotClones = TRUE, plotDrivers = TRUE,
     lwdClone = 0.9, log = "")
}


domcf(1e-6, 0.1, 500, 3)

domcf(1e-6, 0.1, 1000, 3)

domcf(1e-6, 0.07, 1000, 3)

domcf(1e-5, 0.1, 1000, 1)

domcf(5e-6, 0.1, 1000, 1) ## yes, clonal interference clear

domcf(5e-6, 0.1, 2000, 1) ## yes, clonal interference clear, but too short stasis


domcf(5e-6, 0.2, 1000, 1)

domcf(5e-6, 0.08, 1000, 1) ## yes, clonal interference clear


domcf(5e-6, 0.05, 2000, 1) ## yes, clonal interference clear


domcf(5e-6, 0.05, 5000, 1) ## yes, clonal interference clear

domcf(5e-6, 0.2, 5000, 1) ## yes, clonal interference clear

domcf(5e-6, 0.1, 5000, 1) ## yes, clonal interference clear




set.seed(1)
domcf(5e-6, 0.1, 2000, 1) 


set.seed(1)
domcf(5e-6, 0.1, 4000, 1) 


set.seed(1)
domcf(5e-6, 0.1, 6000, 1) 


set.seed(1)
domcf(5e-6, 0.05, 6000, 1) 

## these are neat
set.seed(1)
domcf(1e-6, 0.05, 6000, 1) 
set.seed(2)
domcf(1e-6, 0.05, 6000, 1) 


set.seed(3)
domcf(1e-6, 0.05, 4000, 1) 

set.seed(6)
domcf(1e-6, 0.05, 4000, 1) 

set.seed(6)
domcf(1e-6, 0.045, 2000, 1) 


## beautiful too
set.seed(4)
domcf(1e-6, 0.075, 3000, 1) 
