data(examplePosets)


bozic <- function(poset) oncoSimulIndiv(poset)
bozic9 <- function(poset) oncoSimulIndiv(poset, keepEvery = -9)
    
Exp <- function(poset) oncoSimulIndiv(poset, model = "Exp")
Exp9 <- function(poset) oncoSimulIndiv(poset, model = "Exp", keepEvery = -9)

mc <- function(poset) oncoSimulIndiv(poset, model = "McFL",
                                     mu = 5e-7,
                                     initSize = 4000,
                                     sampleEvery = 0.025,
                                     finalTime = 15000,
                                     keepEvery = 5)
mc9 <- function(poset) oncoSimulIndiv(poset, model = "McFL",
                                     mu = 5e-7,
                                     initSize = 4000,
                                     sampleEvery = 0.025,
                                     finalTime = 15000,
                                      keepEvery = -9)


popsNoZero <- function(x) {
    keep <- which(x$pops.by.time[nrow(x$pops.by.time), ] > 0)
    pops <- x$pops.by.time[nrow(x$pops.by.time), keep]
    genots <- x$Genotypes[, keep[-1] - 1] ## keep contains time too in first col.
    return(list(pops = pops, genots = genots))
}


the.bozics <- function(seed) {
    set.seed(seed); b1 <- bozic(Poset)
    set.seed(seed); b2 <- bozic9(Poset)
    return(list(all = b1, last = b2))
}

the.exps <- function(seed) {
    set.seed(seed); b1 <- Exp(Poset)
    set.seed(seed); b2 <- Exp9(Poset)
    return(list(all = b1, last = b2))
}

the.mcfls <- function(seed) {
    set.seed(seed); b1 <- mc(Poset)
    set.seed(seed); b2 <- mc9(Poset)
    return(list(all = b1, last = b2))
}


## A lot of repetition below. Could make it simpler. Not worth it now?

for(i in 1:length(examplePosets)) {
    s1 <- round(runif(1) * 10000)
    Poset <- examplePosets[[i]]
    
    test_that("Sampling only last same for Bozic", {
        bb <- the.bozics(s1)
        b1 <- bb$all
        b2 <- bb$last
        popsGenots <- popsNoZero(b1)
        expect_equal(b1$TotalPopSize, b2$TotalPopSize)
        expect_equal(b1$FinalTime, b2$FinalTime)
        expect_equal(b1$NumIter, b2$NumIter)
        expect_equal(popsGenots$pops, b2$pops.by.time[1, ])
        expect_equal(popsGenots$genots, b2$Genotypes)
        rm(bb, b1, b2)
    })

    test_that("Sampling only last same for Exp", {
        bb <- the.exps(s1)
        b1 <- bb$all
        b2 <- bb$last
        popsGenots <- popsNoZero(b1)
        expect_equal(b1$TotalPopSize, b2$TotalPopSize)
        expect_equal(b1$FinalTime, b2$FinalTime)
        expect_equal(b1$NumIter, b2$NumIter)
        expect_equal(popsGenots$pops, b2$pops.by.time[1, ])
        expect_equal(popsGenots$genots, b2$Genotypes)
        rm(bb, b1, b2)
    })
    
    test_that("Sampling only last same for Bozic", {
        bb <- the.mcfls(s1)
        b1 <- bb$all
        b2 <- bb$last
        
        popsGenots <- popsNoZero(b1)
        expect_equal(b1$TotalPopSize, b2$TotalPopSize)
        expect_equal(b1$FinalTime, b2$FinalTime)
        expect_equal(b1$NumIter, b2$NumIter)
        expect_equal(popsGenots$pops, b2$pops.by.time[1, ])
        expect_equal(popsGenots$genots, b2$Genotypes)
        rm(bb, b1, b2)
    })
}
