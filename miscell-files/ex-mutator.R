## This looks innocent, but can easily blow up:
pops <- 5
ft <- 500
lni <- 100 
no <- 1e3
ni <- c(0.1, 0.1, 0.1, rep(0, lni))
## scramble around names
names(ni) <- c("a",
               "b",
               "c",
               replicate(lni,
                         paste(sample(letters, 12), collapse = "")))
fe <- allFitnessEffects(noIntGenes = ni)
pg1 <- runif(lni + 3, min = 1e-8, max = 1e-6) 
names(pg1) <- names(ni)
m1 <- allMutatorEffects(noIntGenes = c("b" = 25))
m1.pg1.b <- oncoSimulPop(pops,
                         fe,
                         mu = pg1,
                         muEF = m1,
                         finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no,
                         onlyCancer = FALSE)
m1.pg1.b


## as fewer genes, many fewer possible clones to keep track of
## though, of course, occasionally it will grow huge.
pops <- 5
ft <- 500
lni <- 10 
no <- 1e3
ni <- c(0.1, 0.1, 0.1, rep(0, lni))
## scramble around names
names(ni) <- c("a",
               "b",
               "c",
               replicate(lni,
                         paste(sample(letters, 12), collapse = "")))
fe <- allFitnessEffects(noIntGenes = ni)
pg1 <- runif(lni + 3, min = 1e-8, max = 1e-6) 
names(pg1) <- names(ni)
m1 <- allMutatorEffects(noIntGenes = c("b" = 25))
m1.pg1.b <- oncoSimulPop(pops,
                         fe,
                         mu = pg1,
                         muEF = m1,
                         finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no,
                         onlyCancer = FALSE)
m1.pg1.b
