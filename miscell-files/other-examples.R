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
