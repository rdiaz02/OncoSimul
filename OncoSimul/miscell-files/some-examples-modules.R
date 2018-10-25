## Do modules work? Of course they do. And it can be made equivalent
## to stopping by size.

library(OncoSimulR)
library(graph)
library(igraph)

# Genes per module
genesA <-paste0("a",1:3, collapse =",")
genesB <-paste0("b",1:2, collapse =",")
genesC <-paste0("c",1:4, collapse =",")
genesD <-paste0("d",1:1, collapse =",")
genesE <-paste0("e",1:2, collapse =",")
genesF <-paste0("f",1:3, collapse =",")
genesG <-paste0("g",1:2, collapse =",")


namesA <-paste0("a",1:3)
namesB <-paste0("b",1:2)
namesC <-paste0("c",1:4)
namesD <-paste0("d",1:1)
namesE <-paste0("e",1:2)
namesF <-paste0("f",1:3)
namesG <-paste0("g",1:2)


names <- c(namesA, namesB, namesC, namesD, namesE, namesF, namesG)

#Genero los modulos
modules <- c("Root" = "Root", "A" = genesA,
             "B" = genesB, "C" = genesC,
             "D" = genesD, "E" = genesE,
             "F" = genesF, "G" = genesG)



fdA <- 0.1 
fdB <- 0.15
fdC <- 0.12
fdD <- 0.08
fdE <- 0.075
fdF <- 0.125
fdG <- 0.17


fe1A <- allFitnessEffects(epistasis = c("A" = fdA,
                                        "B" = fdB, "C" = fdC,
                                        "D" = fdD, "E" = fdE,
                                        "F" = fdF, "G" = fdG),
                          geneToModule = modules,
                          drvNames = names)


pi <- 4000 

mcflEv <- function(p, s, initSize) {
  ## expects vectors for p and s
  K <- initSize/(exp(1) - 1)
  
  ## Expected number at equilibrium
  return( K * (exp(prod((1 + s)^p)) - 1))
}



tt <-mcflEv(rep(1, 7), c(fdA, fdB, fdC, fdD, fdE, fdF, fdG), pi)

ttf <-tt-(tt*0.01)



sim1A <- oncoSimulPop(Nindiv = 100, ## tr
                     fp = fe1A,
                     model = "McFL",
                     initSize = pi,
                     detectionSize = ttf,
                     detectionDrivers = NA,
                     detectionProb = NA,
                     sampleEvery = 0.025,
                     keepPhylog = FALSE, 
                     onlyCancer = TRUE,
                     keepEvery = 5,
                     mc.cores = 8) 
summary(sim1A)



## There are smarter ways to do this
## samplePop output -> samplePop output collapsed by module
## but note we can get values > 1. Done on purpose
collapse_by_module <- function(x, listmodules) {
    if(is.null(names(listmodules))) stop("listmodules must have names")
    mm <- do.call(cbind,
                  lapply(listmodules,
                         function(m) {
                             the_cols <- which(colnames(x) %in% m)
                             rowSums(x[, the_cols, drop = FALSE])}
                         )
                  )
    colnames(mm) <- names(listmodules)
    return(mm)
}

single_largest_last_pop <- function(x) {
    last <- nrow(x$pops.by.time)
    largest <- which.max(x$pops.by.time[last, , drop = FALSE]) - 1
    genot <- x[["GenotypesLabels"]][largest]
    strsplit(genot, split = ", ")[[1]]
}


## like samplePop(x, timeSample = "last") but return the single most abundant genotype
largest_last_pop <- function(x) {
    y <- lapply(x, single_largest_last_pop)
    allg <- sort(unique(unlist(y)))
    m <- t(vapply(y, function(z) {as.integer(allg %in%  z) },
                  integer(length(allg)) ))
    colnames(m) <- allg
    return(m)
}



## threshold is important, and the larges clone might not have a frequency
## of 100%
p1 <- samplePop(sim1A, timeSample = "last",
                typeSample = "whole",
                thresholdWhole = 0.9)[, 1:17]
## seems to make sense
colSums(p1)

## Yes, makes sense, though some modules not always present
## in the sampled genotype (probably because its frequency < 0.9)
p1Modules <- collapse_by_module(p1,
                                list(A = namesA, B = namesB,
                                     C = namesC, D = namesD,
                                     E = namesE, F = namesF,
                                     G = namesG))

p1Modules


p2 <- largest_last_pop(sim1A)[, 1:17]
p2Modules <- collapse_by_module(p2,
                                list(A = namesA, B = namesB,
                                     C = namesC, D = namesD,
                                     E = namesE, F = namesF,
                                     G = namesG))
## Hummm.. not always all modules in the most common genotype (depends on run)
p2Modules
which(p2Modules < 1, arr.ind = TRUE) ## nothing here



## Redo simulations allowing for a lag until sampling. Give time
## for the fittest population to really do a clonal sweep.
## This is not necessarily important for the code itself, just to check
## things. Looking at output from simulations above,
## extra time of 100 seems plenty

sim2A <- oncoSimulPop(Nindiv = 200, ## tr
                     fp = fe1A,
                     model = "McFL",
                     initSize = pi,
                     detectionSize = ttf,
                     detectionDrivers = NA,
                     detectionProb = NA,
                     sampleEvery = 0.025,
                     keepPhylog = FALSE, 
                     onlyCancer = TRUE,
                     extraTime = 100,
                     keepEvery = 5,
                     mc.cores = 4) 
summary(sim2A)


p1_2 <- samplePop(sim2A, timeSample = "last",
                typeSample = "whole",
                thresholdWhole = 0.9)[, 1:17]
p1Modules_2 <- collapse_by_module(p1_2,
                                list(A = namesA, B = namesB,
                                     C = namesC, D = namesD,
                                     E = namesE, F = namesF,
                                     G = namesG))
## still a few modules sometimes missing
p1Modules_2


p2_2 <- largest_last_pop(sim2A)[, 1:17]
p2Modules_2 <- collapse_by_module(p2_2,
                                list(A = namesA, B = namesB,
                                     C = namesC, D = namesD,
                                     E = namesE, F = namesF,
                                     G = namesG))
## Always all modules in the most common genotype.
p2Modules_2
which(p2Modules_2 < 1, arr.ind = TRUE) ## nothing here
