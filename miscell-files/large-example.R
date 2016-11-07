## A few examples to illustrate the problems of large population sizes
## and their effects on ti.

## They are commented out for the sake of speed in routing building of the
## vignette.


## Note that in all these cases we are storing the complete relationships
## of clones, with interval keepEvery = 1.


nd <- 10
se <- .2
maxsize <- 1e11 ## Around 1e11, we start getting Recoverable exception ti
## set to DBL_MIN.  as time to next mutation is less than smallest
## positive non-zero float.  You can make that more unlikely if you sample
## more often but eventually you cannot avoid it. Sampling more often also
## means slower execution.  You can also decrease it by making mu smaller
u1 <- allFitnessEffects(noIntGenes = rep(0.05, nd))
unix.time(E1 <-  oncoSimulIndiv(u1, model = "Exp",
                                mu = 1e-7,
                                detectionSize = maxsize, 
                                detectionDrivers = NA,
                                sampleEvery = se, 
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 1000,
                                finalTime = 2000,
                                onlyCancer = TRUE))
summary(E1)

## runs of the example above often throw several recoverable exceptions
## but eventually complete in about 30 seconds on my laptop.
## The number of clones is somewhere around 70 and 110 in most runs. 


### What if we increase the number of drivers to 100?
## Since we know 1e11 leads to the exception, we will decrease detection
## size.

## Note that the number of clones goes up quickly. This contributes a lot to
## increasing running times.
nd <- 15
se <- .2
maxsize <- 5e10 
u2 <- allFitnessEffects(noIntGenes = rep(0.05, nd))
unix.time(E2 <-  oncoSimulIndiv(u2, model = "Exp",
                                mu = 1e-7,
                                detectionSize = maxsize, 
                                detectionDrivers = NA,
                                sampleEvery = se, 
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 1000,
                                finalTime = 2000,
                                onlyCancer = TRUE))
summary(E2)

## With 15 genes it often takes about 5 to 20 seconds in my laptop, and we
## often get between 200 and 400 clones.
print(object.size(E2), units = "MB") ## about 4.5 to 5 MB




nd <- 25
se <- .2
maxsize <- 5e10 
u3 <- allFitnessEffects(noIntGenes = rep(0.05, nd))
unix.time(E3 <-  oncoSimulIndiv(u3, model = "Exp",
                                mu = 1e-7,
                                detectionSize = maxsize, 
                                detectionDrivers = NA,
                                sampleEvery = se, 
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 1000,
                                finalTime = 2000,
                                onlyCancer = TRUE))
summary(E3)

## Now we are up to about a minute or two. And now we get around 1000 or
## more clones. The object is of a size of about
print(object.size(E3), units = "MB")
## 12 MB



nd <- 50
se <- .5
maxsize <- 5e10 
u4 <- allFitnessEffects(noIntGenes = rep(0.05, nd))
unix.time(E4 <-  oncoSimulIndiv(u4, model = "Exp",
                                mu = 1e-7,
                                detectionSize = maxsize, 
                                detectionDrivers = NA,
                                sampleEvery = se, 
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 1000,
                                finalTime = 2000,
                                onlyCancer = TRUE))
summary(E4)
print(object.size(E4), units = "MB")



## Let's try many more genes, most passengers with mildly deleterious effects.
## Set maximum size to 1e10 for speed. We sample every 1 unit, not 0.2
## for speed too. We need to increase max.wall.time, as this will take a while.
nd <- 50
np <- 100
maxsize <- 1e10 
u5 <- allFitnessEffects(noIntGenes = c(rep(0.1, nd), runif(np, -.5, 0)))
unix.time(E5 <-  oncoSimulIndiv(u5, model = "Exp",
                                mu = 1e-7,
                                detectionSize = maxsize, 
                                detectionDrivers = NA,
                                sampleEvery = 1,
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 1000,
                                finalTime = 2000,
                                max.wall.time = 1200,
                                onlyCancer = TRUE))
summary(E5)
print(object.size(E5), units = "MB")

## Simulations take between 300 and over 500 seconds. More than 6500
## clones are produced.


## Use McFarland's model. As the McFarland model has strong carrying
## capacity limits, we set the stopping criterion to be number of
## drivers. Say 0.9*number of drivers. We can only reach large population sizes starting
## from large initial sizes; otherwise, we exit much sooner.
nd <- 50
np <- 100
maxsize <- 1e10 
u6 <- allFitnessEffects(noIntGenes = c(rep(0.1, nd), runif(np, -.2, 0)),
                        drvNames = seq.int(nd))

unix.time(M1 <-  oncoSimulIndiv(u6, model = "McFL",
                                mu = 1e-7,
                                detectionSize = maxsize, 
                                detectionDrivers = 10,
                                sampleEvery = .02,
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 2000,
                                finalTime = 4000,
                                max.wall.time = 500,
                                onlyCancer = TRUE))
summary(M1)
print(object.size(M1), units = "MB")

## That takes less than a second and gives an object of size 4.5 MB. We
## see over 150 to 170 clones.


## Let's start from a larger initSize: 2e6 instead of 2000.
nd <- 50
np <- 100
u6 <- allFitnessEffects(noIntGenes = c(rep(0.1, nd), runif(np, -.2, 0)),
                        drvNames = seq.int(nd))

unix.time(M2 <-  oncoSimulIndiv(u6, model = "McFL",
                                mu = 1e-7,
                                detectionSize = maxsize, 
                                detectionDrivers = 10,
                                sampleEvery = .02,
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 2e6,
                                finalTime = 4000,
                                max.wall.time = 500,
                                onlyCancer = TRUE))
summary(M2)
print(object.size(M2), units = "MB")

## It runs in about 7 seconds, produces an object of 20 to 30 MB, and most
## runs around 4500 to 9000 clones.

## What if we did not keep the clone history? In these cases, the
## variability between runs swamps the effects of keeping or not the
## phylogeny in terms of time and object size.




########################################################

## Let us now use models with more genes. We will start with the McFL
## model. 20000 genes with mildly deleterious genes and 100 genes with
## positive fitness effects.


nd <- 100
np <- 20000
u7 <- allFitnessEffects(noIntGenes = c(rep(0.1, nd), runif(np, -.2, 0)),
                        drvNames = seq.int(nd))
unix.time(M3 <-  oncoSimulIndiv(u7, model = "McFL",
                                mu = 1e-7,
                                detectionSize = maxsize, 
                                detectionDrivers = 10,
                                sampleEvery = .02,
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 2000,
                                finalTime = 4000,
                                max.wall.time = 500,
                                onlyCancer = TRUE))
summary(M3)
print(object.size(M3), units = "GB")

## Takes about 5 seconds, produces and object size of around 1.4 to 2 GB
## and gives about 15000 to 20000 clones.




## Same as before, multiply by 10 the initial population size.
nd <- 100
np <- 20000
u7 <- allFitnessEffects(noIntGenes = c(rep(0.1, nd), runif(np, -.2, 0)),
                        drvNames = seq.int(nd))
unix.time(M5 <-  oncoSimulIndiv(u7, model = "McFL",
                                mu = 1e-7,
                                detectionSize = NA,
                                detectionDrivers = 10,
                                sampleEvery = .02,
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 2e4,
                                finalTime = 4000,
                                max.wall.time = 500,
                                onlyCancer = TRUE))
summary(M5)
print(object.size(M5), units = "GB")
## This takes longer, about 20 seconds, gives an object of size over 4 GB
## and we are up to about 50000 clones. Can we plot it? Sure, but it is a
## very busy plot; let us make it slightly less busy:
plotClonePhylog(M5, N = 100)
## And note that by default we only plot those clones with the N satisfied
## at the last time period (so we do not even try to plot the 50000
## clones). But you could do 
plotClonePhylog(M5, N = 10, t = c(20, 900))
## and you'd get a much, much busier plot.

## Same as before, but multiply by another factor of 10 the initial
## population size.
nd <- 100
np <- 20000
u7 <- allFitnessEffects(noIntGenes = c(rep(0.1, nd), runif(np, -.2, 0)),
                        drvNames = seq.int(nd))
unix.time(M6 <-  oncoSimulIndiv(u7, model = "McFL",
                                mu = 1e-7,
                                detectionSize = NA, 
                                detectionDrivers = 10,
                                sampleEvery = .02,
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 2e5,
                                finalTime = 4000,
                                max.wall.time = 2000,
                                onlyCancer = TRUE))
summary(M6)
print(object.size(M6), units = "GB")
## This takes much, much longer, about 15 minutes (not seconds), gives an
## object of size over 20 GB and we get over 250000 clones in many runs.
## This might not work with a laptop unless you have plenty of RAM.

## What if we keep information every 10 time periods, instead of 1?

nd <- 100
np <- 20000
u7 <- allFitnessEffects(noIntGenes = c(rep(0.1, nd), runif(np, -.2, 0)),
                        drvNames = seq.int(nd))
unix.time(M7 <-  oncoSimulIndiv(u7, model = "McFL",
                                mu = 1e-7,
                                detectionSize = NA, 
                                detectionDrivers = 10,
                                sampleEvery = .02,
                                keepEvery = 10,
                                keepPhylog = TRUE,
                                initSize = 2e5,
                                finalTime = 4000,
                                max.wall.time = 2000,
                                onlyCancer = TRUE))
summary(M7)
print(object.size(M7), units = "GB")

## It takes about 15 minutes, the object is much smaller (6 GB) and we
## only see about 80000 to 90000 clones since the resolution is much
## coarser and many clones are likely to have arisen and gone extinct in
## between those periods of 10 unit times.





## Let us return to more reasonable population sizes in the McFarland
## model but increase the number of passengers, and let's increase the
## number of drivers for detection:
nd <- 100
np <- 50000
u8 <- allFitnessEffects(noIntGenes = c(rep(0.1, nd), runif(np, -.2, 0)),
                        drvNames = seq.int(nd))
unix.time(M8 <-  oncoSimulIndiv(u8, model = "McFL",
                                mu = 1e-7,
                                detectionSize = NA,
                                detectionDrivers = 15,
                                sampleEvery = .02,
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 2e3,
                                finalTime = 4000,
                                max.wall.time = 1200,
                                onlyCancer = TRUE))
summary(M8)
print(object.size(M8), units = "GB")

## That took about 50 seconds, created an object of about 18 GB with
## generally over 80000 to 90000 different clones (of which, yes, we are
## keeping the history also).



nd <- 100
np <- 100000
u9 <- allFitnessEffects(noIntGenes = c(rep(0.1, nd), runif(np, -.2, 0)),
                        drvNames = seq.int(nd))
unix.time(M9 <-  oncoSimulIndiv(u9, model = "McFL",
                                mu = 1e-7,
                                detectionSize = NA,
                                detectionDrivers = 15,
                                sampleEvery = .02,
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 2e3,
                                finalTime = 4000,
                                max.wall.time = 1200,
                                onlyCancer = TRUE))
summary(M9)
print(object.size(M9), units = "GB")









## about here






##########################################################


pancr <- allFitnessEffects(
    data.frame(parent = c("Root", rep("KRAS", 4), 
                   "SMAD4", "CDNK2A", 
                   "TP53", "TP53", "MLL3"),
               child = c("KRAS","SMAD4", "CDNK2A", 
                   "TP53", "MLL3",
                   rep("PXDN", 3), rep("TGFBR2", 2)),
               s = 0.1,
               sh = -0.9,
               typeDep = "MN"),
    drvNames = c("KRAS", "SMAD4", "CDNK2A", "TP53", 
                 "MLL3", "TGFBR2", "PXDN"))

se <- .5
maxsize <- 1e11 ## Around 1e11, we start getting Recoverable exception ti
## set to DBL_MIN.  as time to next mutation is less than smallest
## positive non-zero float.  You can make that more unlikely if you sample
## more often but eventually you cannot avoid it. Sampling more often also
## means slower execution.  
unix.time(E6 <-  oncoSimulIndiv(pancr,
                                model = "Exp",
                                mu = 1e-7,
                                detectionSize = maxsize, 
                                detectionDrivers = NA,
                                sampleEvery = se, 
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 1000,
                                finalTime = 2000,
                                onlyCancer = TRUE))
summary(E6)


## Let us run the McFarland model. Because of the strong carrying capacity
## limits with the McFarland model, if we want to reach really large
## population sizes we must start from a large population size. Otherwise,
## we reach the largest possible genotype well before a large population size.


initsize <- 1e10 ## Beyond this, we often get the ti set to DBL_MIN issue.
unix.time(M2 <-  oncoSimulIndiv(pancr,
                                model = "McFL",
                                mu = 1e-7,
                                detectionSize = NA, 
                                detectionDrivers = 7,
                                sampleEvery = 0.03, 
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = initsize,
                                finalTime = 2000,
                                onlyCancer = TRUE))
summary(M2)




#######################
nd <- 5
se <- .1
maxsize <- 1e12 ## Around 1e11, we start getting Recoverable exception ti
## set to DBL_MIN.  as time to next mutation is less than smallest
## positive non-zero float.  You can make that more unlikely if you sample
## more often but eventually you cannot avoid it. Sampling more often also
## means slower execution.  You can also decrease it by making mu smaller
u1 <- allFitnessEffects(noIntGenes = rep(0.05, nd))
unix.time(E3 <-  oncoSimulIndiv(u1, model = "Exp",
                                mu = 1e-7,
                                detectionSize = maxsize, 
                                detectionDrivers = NA,
                                sampleEvery = se, 
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 1000,
                                finalTime = 2000,
                                onlyCancer = TRUE))
summary(E3)












nd <- 10
np <- 0
u1 <- allFitnessEffects(noIntGenes = c(rep(0.05, nd), runif(np, -.2, 0)))

unix.time(E1 <-  oncoSimulIndiv(u1, model = "Exp",
                                mu = 1e-8,
                                detectionSize = 1e12, 
                                detectionDrivers = NA,
                                sampleEvery = .001, ## avoid the ti DBL_MIN
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 1000,
                                finalTime = 2000,
                                ## mutationPropGrowth = FALSE,
                                onlyCancer = TRUE))
summary(E1)



set.seed(123)

nd <- 100
np <- 400
u1 <- allFitnessEffects(noIntGenes = c(rep(0.1, nd), runif(np, -.5, 0)))

unix.time(E1 <-  oncoSimulIndiv(u1, model = "Exp",
                                mu = 1e-7,
                                detectionSize = 5e9, 
                                detectionDrivers = NA,
                                sampleEvery = 1,
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 1000,
                                finalTime = 2000,
                                onlyCancer = TRUE))
summary(E1)








set.seed(123)
nd <- 70  
np <- 50000 
s <- 0.1  
sp <- 1e-4 ## as we have many more passengers
spp <- -sp/(1 + sp)
mcfL <- allFitnessEffects(noIntGenes = c(rep(s, nd), rep(spp, np)),
                          drvNames = seq.int(nd))







## The McFarland model cannot reach 1e14 sizes
unix.time(mcL1 <-  oncoSimulIndiv(mcfL,
                         model = "McFL", 
                         mu = 1e-7,
                         detectionSize = 5e14, 
                         detectionDrivers = NA,
                         sampleEvery = 0.02,
                         keepEvery = 1,
                         keepPhylog = TRUE,
                         initSize = 1000,
                         finalTime = 2000,
                         onlyCancer = FALSE))
## 7 seconds

## Exponential

unix.time(E1 <-  oncoSimulIndiv(mcfL,
                         model = "Exp", 
                         mu = 1e-7,
                         detectionSize = 5e6, 
                         detectionDrivers = NA,
                         sampleEvery = 0.02,
                         keepEvery = 1,
                         keepPhylog = TRUE,
                         initSize = 1000,
                         finalTime = 2000,
                         onlyCancer = FALSE))


## 4.5 seconds

unix.time(E2 <-  oncoSimulIndiv(mcfL, model = "Exp",
                                mu = 1e-7,
                                detectionSize = 5e6, 
                                detectionDrivers = NA,
                                sampleEvery = 1,
                                keepEvery = 1,
                                keepPhylog = TRUE,
                                initSize = 10000,
                                finalTime = 2000,
                                onlyCancer = FALSE))
