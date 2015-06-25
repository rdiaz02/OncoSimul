s1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = 0.1,
                 sh = c(rep(0, 4), c(-.9, -.9), rep(-.95, 3)),
                 typeDep = "SM")
smn1 <- allFitnessEffects(s1)
plot(smn1)
set.seed(123)
tmp <- oncoSimulIndiv(smn1, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("b, a")
                      )
plotClonePhylog(tmp, N = 0)


plot(tmp)



tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
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
plotClonePhylog(tmp, N = 0)


s1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = 0.1,
                 sh = c(rep(0, 4), c(-.9, -.9), rep(-.95, 3)),
                 typeDep = "SM")
smn1 <- allFitnessEffects(s1)
plot(smn1)


set.seed(123)


par(ask = TRUE)
par(mfrow = c(1, 2))

while(TRUE) {
    tmp <- oncoSimulIndiv(smn1, model = "McFL",
                          mu = 5e-5, finalTime = 500,
                          detectionDrivers = 3,
                          onlyCancer = FALSE,
                          initSize = 1000, keepPhylog = TRUE)
    plotClonePhylog(tmp, N = 0)
    plot(tmp)
}



## animation like
for(i in seq(from = 800, to = 1000, by = 50))
    plotClonePhylog(mcf1s, N = 1, t = c(i, i + 5))



## a silly thing to catch clones that appear in phylog but never
## in popps.by.time
for(i in 1:100){
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
                           finalTime = 20000,
                           onlyCancer = FALSE,
                           extraTime = 1500,
                           keepPhylog = TRUE);
    a <- sort(tmp$GenotypesLabels)
    b <- sort(unique(as.character(unlist(tmp$other$PhylogDF[, c(1, 2)]))))
    if(!(all(a == b))) stop("here")
}



library(OncoSimulR)

data(examplesFitnessEffects)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
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
summary(tmp)



library(OncoSimulR)
data(examplesFitnessEffects)
tmp <-  oncoSimulIndiv(examplesFitnessEffects[["o3"]],
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



## many epistasis; is this silly?



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
see <- oncoSimulIndiv(ee, model = "Exp", detectionDrivers = 1000,
                      sampleEvery = 10,
                      keepEvery = -9)
summary(see)



##### This is silly in a vignette
## %% \section{Running all}

## %% <<>>=
## %% for(i in 1:length(examplesFitnessEffects)) {
## %%     cat(paste("\n Doing i = ", i , " name = ",
## %%               names(examplesFitnessEffects)[i], "\n"))
## %%     tmp <-  oncoSimulIndiv(examplesFitnessEffects[[i]],
## %%                            model = "Bozic", 
## %%                            mu = 1e-6,
## %%                            detectionSize = 1e8, 
## %%                            detectionDrivers = 4,
## %%                            sampleEvery = 2,
## %%                            max.num.tries = 100,
## %%                            initSize = 2000,
## %%                            onlyCancer = TRUE)
## %% }

## %% for(i in 1:length(examplesFitnessEffects)) {
## %%     cat(paste("\n Doing i = ", i , " name = ",
## %%               names(examplesFitnessEffects)[i], "\n"))
## %%     tmp <-  oncoSimulIndiv(examplesFitnessEffects[[i]],
## %%                            model = "Exp", 
## %%                            mu = 1e-6,
## %%                            detectionSize = 1e8, 
## %%                            detectionDrivers = 4,
## %%                            sampleEvery = 2,
## %%                            max.num.tries = 100,
## %%                            initSize = 2000,
## %%                            onlyCancer = TRUE)
## %% }


## %% for(i in 1:length(examplesFitnessEffects)) {
## %%     cat(paste("\n Doing i = ", i , " name = ",
## %%               names(examplesFitnessEffects)[i], "\n"))
## %%     tmp <-  oncoSimulIndiv(examplesFitnessEffects[[i]],
## %%                            model = "McFL", 
## %%                            mu = 5e-7,
## %%                            detectionSize = 1e8, 
## %%                            detectionDrivers = 2,
## %%                            sampleEvery = 0.025,
## %%                            max.num.tries = 10,
## %%                            initSize = 2000,
## %%                            finalTime = 15000,
## %%                            onlyCancer = TRUE)
## %% }
## %% @ 
