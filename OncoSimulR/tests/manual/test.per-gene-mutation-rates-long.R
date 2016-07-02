cat("\n Starting per-gene-mutation-rates long at", date(), "\n") ## whole file takes about 30 seconds
## When submitting, probably move half of the tests (mcfl?) to the "long"
## file.

## FIXME wrap some of the p-value based tests on a loop to catch
## occasional mistakes. See, e.g., test.mutator-oncoSimulSample.R.

mutsPerClone <- function(x, per.pop.mean = TRUE) {
    perCl <- function(z)
        unlist(lapply(z$GenotypesWDistinctOrderEff, length))
    perCl2 <- function(z)
        mean(unlist(lapply(z$GenotypesWDistinctOrderEff, length)))

    if(per.pop.mean)    
        unlist(lapply(x, function(u) perCl2(u)))
    else
        lapply(x, function(u) perCl(u))
}



mutsPerCloneLast <- function(x, per.pop.mean = TRUE) {
    ## Only clones alive in the last period
    aliveLast <- function(u) {
        pbt <- u$pops.by.time
        which(pbt[nrow(pbt), -1] >= 1)
    }
    perCl <- function(z) {
        this <- aliveLast(z)
        unlist(lapply(z$GenotypesWDistinctOrderEff[this],  length))
    }
    perCl2 <- function(z) {
        this <- aliveLast(z)
        mean(unlist(lapply(z$GenotypesWDistinctOrderEff[this], length)))
    }
    if(per.pop.mean)    
        unlist(lapply(x, function(u) perCl2(u)))
    else
        lapply(x, function(u) perCl(u))
}

## we could have used this below . Oh well
## totalind <- function(out) {
##     ## total num indivs
##   sum(unlist(lapply(out, function(x) x$TotalPopSize)))  
## }


## RNGkind("L'Ecuyer-CMRG") ## for the mclapplies
## RNGkind("Mersenne-Twister")

p.value.threshold <- 0.01

date()
test_that("Same freqs, chisq, when s=0", {
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
   
    
    cat("\n s08: a runif is", runif(1), "\n")
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = 0.001, max.wall.time = 900,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ## colSums(OncoSimulR:::geneCounts(bb))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)))
        ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
         
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()

date()
test_that("Same freqs, chisq, when s", {
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n s09: a runif is", runif(1), "\n")
    muvar2 <- c("U" = 1e-5, "z" = 1e-5, "e" = 1e-5, "m" = 1e-5, "D" = 1e-5)
    ni1 <- rep(0.02, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = 0.001,max.wall.time = 900,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ## colSums(OncoSimulR:::geneCounts(bb))
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
     })
date()

date()
test_that("Different freqs as they should be ordered and chisq, when s=0",{
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n s010: a runif is", runif(1), "\n")
    ## muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ## too small mu: observed often 0 and expected below 1. Make larger
    muvar2 <- c("U" = 1e-3, "z" = 5e-5, "e" = 5e-4, "m" = 5e-3, "D" = 1e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 5e7
    reps <- 800
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = .00001, max.wall.time = 900,
                       seed =NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
    ## the ordering criterion could fail even if things are OK
    TTT <- c(TTT, identical(
                      order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC)))    
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
     })
date()

date()
test_that("Different freqs as they should be ordered, when s and t > 1", {
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n s011: a runif is", runif(1), "\n")
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 5e5
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = 2, max.wall.time = 900,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    ## A chisq will not work as we increase finalTime. But ordering of
    ## freqs. should.
    TTT <- c(TTT, identical(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC)))
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     TTT <- c(TTT, all(TTT))
})
date()

date()
test_that("Different freqs as they should be ordered, when s and t> 1, again", {
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n s012: a runif is", runif(1), "\n")
    muvar2 <- c("U" = 1e-3, "z" = 5e-5, "e" = 5e-4, "m" = 5e-3, "D" = 1e-4)
    ## we increase s and finalTime
    ni2 <- rep(0.1, 5)
    names(ni2) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 1e5
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = 4, max.wall.time = 900,
                       mutationPropGrowth = FALSE, ## cleaner, tough no real effect
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    TTT <- c(TTT, identical(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC)))
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()

date()
test_that("Complex fitness specification, s diffs, tiny finalTime, systematic mu", {
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n s1: a runif is", runif(1), "\n")
    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                     child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                     s = c(0.1, 0.2, 0.3, 0.4, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                     sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                     typeDep = c(rep("--", 4), 
                                 "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    oe <- c("C > F" = -0.1, "H > I" = 0.12)
    sm <- c("I:J"  = -.1)
    sv <- c("-K:M" = -.5, "K:-M" = -.5)
    epist <- c(sm, sv)
    modules <- c("Root" = "Root", "A" = "a1",
                 "B" = "b1, b2", "C" = "c1",
                 "D" = "d1, d2", "E" = "e1",
                 "F" = "f1, f2", "G" = "g1",
                 "H" = "h1, h2", "I" = "i1",
                 "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")
    noint <- runif(5, min = 0.051, max = 0.1)
    names(noint) <- paste0("n", 1:5)
    drvN <- paste0(letters[c(1:11, 13, 2, 4, 6, 8, 10, 11)],
               c(rep(1, 12), rep(2, 6)))
    fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                             noIntGenes = noint, geneToModule = modules,
                             drvNames = drvN)
    nfea <- OncoSimulR:::allNamedGenes(fea)$Gene
    ## systematic spacing
    muvar <- sample(seq(from = 5e-6, to = 1e-3, length.out = length(nfea)))
    names(muvar) <- nfea
    no <- 5e7
    reps <- 300
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = 0.0001, max.wall.time = 900,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar)
    colSums(OncoSimulR:::geneCounts(bb))
    ## expectedC - colSums(OncoSimulR:::geneCounts(bb))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##                        p = expectedC/sum(expectedC))
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
    ## Order will not necessarily be OK, because of random fluctuations,
    ## when muvar diffs are tiny
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))),
    ##     order(expectedC))
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()



date()
test_that("Complex fitness specification, tiny s diffs", {
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n s2: a runif is", runif(1), "\n")
    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                     child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                     s = c(0.00001, 0.00002, 0.00003, 0.00004, 0.00001, 0.00001, 0.00002, 0.00002, 0.00003, 0.00003),
                     sh = c(rep(0, 4), c(-.0000009, -.0000009), c(-.00000095, -.00000095), c(-.00000099, -.00000099)),
                     typeDep = c(rep("--", 4), 
                                 "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    oe <- c("C > F" = -0.00001, "H > I" = 0.000012)
    sm <- c("I:J"  = -.00001)
    sv <- c("-K:M" = -.000005, "K:-M" = -.000005)
    epist <- c(sm, sv)
    modules <- c("Root" = "Root", "A" = "a1",
                 "B" = "b1, b2", "C" = "c1",
                 "D" = "d1, d2", "E" = "e1",
                 "F" = "f1, f2", "G" = "g1",
                 "H" = "h1, h2", "I" = "i1",
                 "J" = "j1, j2", "K" = "k1, k2", "M" = "m1")
    noint <- runif(5, min = 0.0000051, max = 0.00001)
    names(noint) <- paste0("n", 1:5)
    drvN <- paste0(letters[c(1:11, 13, 2, 4, 6, 8, 10, 11)],
                   c(rep(1, 12), rep(2, 6)))
    fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                             noIntGenes = noint, geneToModule = modules,
                             drvNames = drvN)
    nfea <- OncoSimulR:::allNamedGenes(fea)$Gene
    ## Now, random muvar
    ## muvar <- sample(seq(from = 5e-6, to = 1e-3, length.out = length(nfea)))
    muvar <- runif(length(nfea), min = 5e-6, max = 1e-3)
    names(muvar) <- nfea
    no <- 5e7
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = .0001,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar)
    colSums(OncoSimulR:::geneCounts(bb))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##                        p = expectedC/sum(expectedC))
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value > p.fail)
    ## Even with systematic spacing, you need huge reps to even out the
    ## sampling effects on order. And ordering tested above several
    ## times. This is an overkill.
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))),
    ##     order(expectedC))
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()




date()
test_that("Init mutant  with tiny mutation always present", {
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n s11: a runif is", runif(1), "\n")
    ## We check the two init mutants have same frequency, and it is the
    ## largest. So we also verify no additional mutation to
    ## either one (since their freq is same) and we verify present in all.
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
    noint <- runif(5, min = 0.01, max = 0.1)
    names(noint) <- paste0("n", 1:5)
    drvN <- paste0(letters[c(1:11, 13, 2, 4, 6, 8, 10, 11)],
                   c(rep(1, 12), rep(2, 6)))
    fea <- allFitnessEffects(rT = p4, epistasis = epist, orderEffects = oe,
                             noIntGenes = noint, geneToModule = modules,
                             drvNames = drvN)
    nfea <- OncoSimulR:::allNamedGenes(fea)$Gene
    ## muvar <- runif(length(nfea), min = 1e-7, max = 1e-3) ## too tiny
    ## diffs sometimes for order comp
    muvar <- sample(seq(from = 1e-7, to = 1e-4, length.out = length(nfea)))
    names(muvar) <- nfea
    muvar["h2"] <- 1e-13
    muvar["i1"] <- 2e-13
    no <- 1e5
    reps <- 20
    bb <- oncoSimulPop(reps,
                       fea, mu = muvar,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       sampleEvery = 0.03,
                       keepEvery = 2,
                       finalTime = 10000, ## yes, huge; we only get close to 4 or 5.
                       detectionDrivers = 4,
                       mutationPropGrowth = FALSE, ## yes, exclude this possible effect
                       initMutant = "h2 > i1",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar)
    ccs <- colSums(OncoSimulR:::geneCounts(bb))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["h2"] == ccs["i1"])
    expect_true(ccs["h2"] == totalindivs)
    TTT <- c(TTT, all(ccs["h2"] > ccs[!(names(ccs) %in% c("h2", "i1"))]))
    ## of course, no agreement with chi-square or ordering
    p.fail <- 1e-6
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
                           p = expectedC/sum(expectedC))$p.value < p.fail)
    TTT <- c(TTT, !(
        identical(
            order(colSums(OncoSimulR:::geneCounts(bb))),
            order(expectedC))))
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()



date()
test_that("Different freqs as they should be ordered and chisq, when s  and a tiny mu", {
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n s13: a runif is", runif(1), "\n")
    muvar2 <- c("U" = 1e-13, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    ni2["U"] <- 0.5
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 1e6
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = 5,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    TTT <- c(TTT, colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    TTT <- c(TTT, identical(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC)))
    ## A chisq will not work as we increase finalTime.
    
    
    cat("\n s13b: a runif is", runif(1), "\n")
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = .001,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    TTT <- c(TTT, colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    ## This will fail sometimes
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-1],
                           p = expectedC[-1]/sum(expectedC))$p.value > p.fail)
    ## could fail because very small freqs
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))),
         ##     order(expectedC))
          if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
}
)
date()



date()
test_that("Different freqs as they should be ordered and chisq, when s=0, and initMutant",{
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n s17: a runif is", runif(1), "\n")
    muvar2 <- c("U" = 1e-3, "z" = 5e-5, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = .001,
                       seed =NULL,
                       initMutant = "e",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ccs <- colSums(OncoSimulR:::geneCounts(bb))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
    TTT <- c(TTT, identical(
        order(colSums(OncoSimulR:::geneCounts(bb))[-3]),
        order(expectedC[-3])))
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()




date()
test_that("Different freqs as they are expected with chisq, when s=0 and initMutant, many genotypes",{
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n s19: a runif is", runif(1), "\n")
    ft <- 0.001 ## yes, small
    lni <- 100  ## 16
    muvar2 <- runif(lni, min = 1e-4, max = 1e-3)
    names(muvar2) <- c(replicate(lni,
                                 paste(sample(letters, 12), collapse = "")))
    names(muvar2)[3] <- "e"
    muvar2[3] <- 1e-9
    ni1 <- rep(0, lni)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 100
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = ft,
                       seed =NULL,
                       initMutant = "e",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    (ccs <- colSums(OncoSimulR:::geneCounts(bb)))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()




date()
test_that("McFL: Different freqs as they should be ordered and chisq, when s=0, and initMutant",{
    ## More on the above, with less variation. But yet another set of tests.
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n mcs20: a runif is", runif(1), "\n")
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ## moderately small mu
    muvar2[] <- 1e-5
    muvar2["e"] <- 1e-3
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = 0.002,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       initMutant = "e",
                       model = "McFL",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ccs <- colSums(OncoSimulR:::geneCounts(bb))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ## relatively large mu
    muvar2[] <- 1e-3
    muvar2["e"] <- 1e-6
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    
    
    cat("\n mcs20b: a runif is", runif(1), "\n")
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = 0.002,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       initMutant = "e",
                       model = "McFL",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    ccs <- colSums(OncoSimulR:::geneCounts(bb))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
    ## nope, as many are equal
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))[-3]),
    ##     order(expectedC[-3]))
    ## chisq.test(colSums(OncoSimulR:::geneCounts(bb)),
    ##            p = expectedC/sum(expectedC))
    ## some moderate, one very large
    muvar2 <- c("U" = 1e-3, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    muvar2[] <- 1e-4
    muvar2["e"] <- 1e-6
    muvar2[4] <- 1e-2
    ni1 <- rep(0, 5)
    names(ni1) <- names(muvar2)
    fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e7
    reps <- 200
    
    
    cat("\n mcs20c: a runif is", runif(1), "\n")
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = 0.002,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       initMutant = "e",
                       model = "McFL",
                       mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    (ccs <- colSums(OncoSimulR:::geneCounts(bb)))
    totalindivs <- sum(unlist(lapply(bb, function(x) x$TotalPopSize)))
    expect_true(ccs["e"] == totalindivs)
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-3],
                           p = expectedC[-3]/sum(expectedC[-3]))$p.value > p.fail)
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()



date()
test_that("Num clones, muts per clone for different per-gene-mut",{
    ## Like previous, but larger finalTime, so no longer chi-square test
    ## here.
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n df2: a runif is", runif(1), "\n")
    ng <- 40
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-6, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 2
    no <- 1e5
    reps <- 20
    
    
    cat("\n df2a: a runif is", runif(1), "\n")
    b1 <- oncoSimulPop(reps,
                       fe1,
                       mu = m1,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       mc.cores = 2, max.wall.time = 900
                       )
    
    
    cat("\n df2b: a runif is", runif(1), "\n")
    b2 <- oncoSimulPop(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       mc.cores = 2, max.wall.time = 900
                       )
    TTT <- c(TTT,  wilcox.test(summary(b2)$NumClones, 
                 summary(b1)$NumClones, alternative = "greater")$p.value < p.value.threshold)
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
    TTT <- c(TTT,  t.test(mutsPerClone(b2),
                 mutsPerClone(b1), alternative = "greater")$p.value < p.value.threshold)
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()





date()
test_that("McFL: Expect freqs, num clones, muts per clone for different per-gene-mut",{
## More on the above, with less variation. But yet another set of tests.
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n mcdf1: a runif is", runif(1), "\n")
    ng <- 10
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-6, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 0.01
    no <- 1e7
    reps <- 80
    
    
    cat("\n mcdf1a: a runif is", runif(1), "\n")
    b1 <- oncoSimulPop(reps,
                       fe1,
                       mu = m1,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       model = "McFL",
                       mc.cores = 2, max.wall.time = 900
                       )
    
    
    cat("\n mcdf1b: a runif is", runif(1), "\n")
    b2 <- oncoSimulPop(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       model = "McFL",
                       mc.cores = 2, max.wall.time = 900
                       )
    (expected1 <- no*reps*m1)
    (expected2 <- no*reps*m2)
    (cc1 <- colSums(OncoSimulR:::geneCounts(b1)))
    (cc2 <- colSums(OncoSimulR:::geneCounts(b2)))    
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(b1)),
                           p = expected1/sum(expected1))$p.value > p.fail)
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(b2)),
                           p = expected2/sum(expected2))$p.value > p.fail)
    TTT <- c(TTT,  wilcox.test(summary(b2)$NumClones, 
                 summary(b1)$NumClones, alternative = "greater")$p.value < p.value.threshold)
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
    TTT <- c(TTT,  t.test(mutsPerClone(b2),
                 mutsPerClone(b1), alternative = "greater")$p.value < p.value.threshold)
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()



date()
test_that(" And mutPropGrowth, 3",{
    ## This is extremely variable, so the pattern is hard to catch as few genes.
    ## And one can even see apparently counterintuitive patterns if
    ## if you make the number of genes tiny and s very large.
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
         cat("\n sz033: a runif is", runif(1), "\n")
         
         muvar2 <- c("U" = 1e-4, "Z" = 5e-5, "E" = 5e-4, "M" = 5e-3, "D" = 1e-3)
         muvar2 <- c(muvar2, rep(1e-7, 10))
         names(muvar2)[6:15] <- letters[1:10]
         ## muvar2 <- c("U" = 5e-5, "z" = 5e-5, "e" = 5e-5, "m" = 5e-5, "D" = 5e-5)
         ni1 <- rep(.5, 5)  ## 1.9
         ni1 <- c(ni1, rep(0, 10))
         names(ni1) <- names(muvar2)
         fe1 <- allFitnessEffects(noIntGenes = ni1)
    no <- 1e6
    reps <- 25 
    ft <- 36 ## irrelevant, we stop on size
    cat("\n sz033a: a runif is", runif(1), "\n")
    b1 <- oncoSimulPop(reps,
                       fe1, mu = muvar2,
                       mutationPropGrowth = FALSE,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = ft,
                       detectionSize = 1e8,
                       sampleEvery = 0.05,
                       keepEvery = 1,
                       seed =NULL,
                       mc.cores = 2, max.wall.time = 900
                       )
    cat("\n sz033b: a runif is", runif(1), "\n")
    b2 <- oncoSimulPop(reps,
                       fe1, mu = muvar2,
                       mutationPropGrowth = TRUE,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = ft,
                       detectionSize = 1e8,
                       sampleEvery = 0.05,
                       keepEvery = 1,
                       seed =NULL,
                       mc.cores = 2, max.wall.time = 900
                       )
         summary(b1)[, c(1:3, 8:9)]
         summary(b2)[, c(1:3, 8:9)]
         mean(mutsPerClone(b1));mean(mutsPerClone(b2))
    median(summary(b1)$NumClones)
    median(summary(b2)$NumClones)

         ## More mutations in mutationPropGrowth
    TTT <- c(TTT,  t.test(mutsPerClone(b2),
                 mutsPerClone(b1), alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(summary(b2)$NumClones, 
                 summary(b1)$NumClones, alternative = "greater")$p.value < p.value.threshold)
    ## But frequency of mutations about the same? Nope: since very fast
    ## growth and thus non-indep, huge variation in geneCounts in each
    ## run, etc. so hard to compare geneCounts
    ## Just for reference, here
    ## First, look at run to run variation
    ## OncoSimulR:::geneCounts(b1)
    ## OncoSimulR:::geneCounts(b2)
    ## ## The next makes sense
    ## fb1 <- colSums(OncoSimulR:::geneCounts(b1))
    ## fb2 <- colSums(OncoSimulR:::geneCounts(b2))
    ## fb1
    ## fb2
    ## fb2/fb1
    ## fb1/sum(fb1)
    ## fb2/sum(fb2)
    ## (fb2/sum(fb2))/(fb1/sum(fb1))
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()




date()
test_that("McFL: More mutpropgrowth, in modules of s", {
    ## From a similar test in mutPropGrowth, but we have a vector mu

    ## With McFL, since total size is bounded from above, the effects of
    ## mutPropGrowth are rarely those of popSize. But we stop on popSize
    ## anyway here.
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n mcmpgs3: a runif is", runif(1), "\n")
    pops <- 20
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e3
    ft <-  15 ## 5
    s3 <- 3.0
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    names(ni) <- paste0("ni", 1:lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f3 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    mu <- runif(fni + lni, min = 1e-7, max = 1e-4)
    names(mu) <- c(paste0("a", 1:fni), names(ni))
    
    
    cat("\n mcmpgs3a: a runif is", runif(1), "\n")
    s3.ng <- oncoSimulPop(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          sampleEvery = 0.01,
                          detectionSize = 2.5e4,
                          detectionDrivers = 9999,
                          initSize = no,
                          onlyCancer = FALSE, detectionProb = NA,
                          model = "McFL", max.wall.time = 900,
                          seed = NULL, mc.cores = 2)
    
    
    cat("\n mcmpgs3b: a runif is", runif(1), "\n")
    s3.g <- oncoSimulPop(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         sampleEvery = 0.01,
                         detectionSize = 2.5e4,
                         detectionDrivers = 9999,
                         initSize = no,
                         onlyCancer = FALSE, detectionProb = NA,
                         model = "McFL",  max.wall.time = 900,
                         seed = NULL, mc.cores = 2)
    summary(s3.g)[, c(1, 2, 3, 8, 9)]
    summary(s3.ng)[, c(1, 2, 3, 8, 9)]
    summary(summary(s3.ng)[, 2])
    summary(summary(s3.g)[, 2])
         TTT <- c(TTT,  t.test(mutsPerClone(s3.g),
                               mutsPerClone(s3.ng), alternative = "greater")$p.value < p.value.threshold)
         TTT <- c(TTT,  wilcox.test(summary(s3.g)$NumClones, 
                                    summary(s3.ng)$NumClones, alternative = "greater")$p.value < p.value.threshold)
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()




date()
test_that("oncoSimulSample: expected vs. observed for different per-gene-mut",{
    ## Here, we test that freqs as they should, but so that the test is
    ## not eternal, we use different settings of reps and no

    ## We want to get about 2 mutations in each population, but not more
    ## than 2, as that would mean we are not picking diffrences between
    ## muts with different mut rate, because we have
    ## wholePopulationSample.
    
    ## We probably want about a mean or median number of clones of about 2
    ## or so. Though if fewer, better but then to have power in the
    ## chi-square we need much larger reps (as usual, if ft increase, etc,
    ## we increase the reproduction/death events, which then screws up
    ## simple expectations for chi-square)
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n oss11: a runif is", runif(1), "\n")
    ng <- 10
    ni <- rep(0, ng)
    m1 <- seq(from = 1e-7, to = 1e-4, length.out = ng)
    m2 <- runif(ng, min = 1e-6, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 5e-3
    no <- 5e5 # delicate as if this is huge, we get the cc1 or cc2 below
              # to be equal to reps in many genes, because they are
              # present in at least one cell in all populations
    reps <- 5000 ## large because few events with small mut
                 ## freqs. o.w. chi has many small cells.
    
    
    x <- 1e-20
    cat("\n oss1a: a runif is", runif(1), "\n")
    b1 <- oncoSimulSample(reps,
                          fe1,
                          mu = m1,
                          onlyCancer = FALSE, detectionProb = NA,
                          initSize = no,
                          finalTime = ft,
                          mutationPropGrowth = FALSE, ## cleaner, though no real effect
                          seed =NULL,
                          sampleEvery = 0.01, thresholdWhole = 1e-20,
                          detectionSize = 1e9,
                          detectionDrivers = 9999,
                          )
    b1$popSummary[, c(1:3, 8:9)]
    summary(b1$popSummary[, "NumClones"])
    (expected1 <- no*reps*m1)
    (cc1 <- colSums(b1$popSample))
    if( (any(cc1 == reps)) )
        warning("The test is likely to fail because reps == cc1 or cc2")
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(cc1,
                           p = expected1/sum(expected1))$p.value > p.fail)

    reps <- 500
    no <- 1e4
    ft <- 0.03
    
    
    cat("\n oss1b: a runif is", runif(1), "\n")
    b2 <- oncoSimulSample(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       sampleEvery = 0.01, thresholdWhole = 1e-20,
                       detectionSize = 1e9,
                       detectionDrivers = 9999,
                       )
    summary(b2$popSummary[, "NumClones"])
    ## we detect anything that is present in at least one case.
    ## Not exactly the same as what we did in oncoSimulPop
    (expected2 <- no*reps*m2)
    (cc2 <- colSums(b2$popSample))
    if( (any(cc2 == reps)))
        warning("The test is likely to fail because reps == cc1 or cc2")
    ## It will fail with prob ~ p.fail
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(cc2,
                           p = expected2/sum(expected2))$p.value > p.fail)


         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()





date()
test_that("oncoSimulSample comparing different per-gene-mut",{
    ## No attempt to compare against expected (other tests do that). We
    ## just verify that larger mutations rates lead to more total
    ## mutations and clones.
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n oss2: a runif is", runif(1), "\n")
    ng <- 10
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-6, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- .05 ## if you make it too large, then all pops will have at
             ## least one cell with one of the genes mutated. You can see
             ## this when cc1 or cc2 have most/all entries equal to reps.
    no <- 1e5 # delicate as if this is huge, we get the cc1 or cc2 below
              # to be equal to reps in many genes, because they are
              # present in at least one cell in all populations
    reps <- 500
    x <- 1e-20
    
    
    cat("\n oss2a: a runif is", runif(1), "\n")
    b1 <- oncoSimulSample(reps,
                          fe1,
                          mu = m1,
                          onlyCancer = FALSE, detectionProb = NA,
                          initSize = no,
                          finalTime = ft,
                          mutationPropGrowth = FALSE, ## cleaner, though no real effect
                          seed =NULL,
                          thresholdWhole = x
                          )
    
    
    cat("\n oss2b: a runif is", runif(1), "\n")
    b2 <- oncoSimulSample(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = FALSE, ## cleaner, though no real effect
                       seed =NULL,
                       thresholdWhole = x
                       )
    ## we detect anything that is present in at least one case.
    ## Not exactly the same as what we did in oncoSimulPop
    (cc1 <- colSums(b1$popSample))
    (cc2 <- colSums(b2$popSample))
    TTT <- c(TTT, sum(cc2) > sum(cc1))
    ## This is very similar to above, like assimilating a pop to a clone
    mutsPerClone1 <- rowSums(b1$popSample)
    mutsPerClone2 <- rowSums(b2$popSample)
    TTT <- c(TTT,  t.test(mutsPerClone2, 
                          mutsPerClone1, alternative = "greater")$p.value < p.value.threshold)
         TTT <- c(TTT,  wilcox.test(b2$popSummary[, "NumClones"], 
                 b1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.value.threshold)
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()




## From similar tests in mutPropGrwoth-long, but here mu is a vector

cat("\n", date(), "\n")
test_that("oncoSimulSample Without initmutant and modules, fixed size", {
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n osSFPS: a runif is", runif(1), "\n")
    pops <- 60
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e4 
    ft <- 9  #4 
    s3 <- 2.5 
    mu <- 1e-5 
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    names(ni) <- paste0("ni", 1:lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste0("a", 1:fni)
    mu <- runif(lni + fni, min = 1e-7, max = 1e-4)
    names(mu) <- c(gn, names(ni))
    gn <- paste(gn, collapse = ", ")
    f3 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
   
    x <- 1e-9 ## so basically anything that appears once
    
    
    cat("\n osSFPSa: a runif is", runif(1), "\n")
    b1 <- oncoSimulSample(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE, detectionProb = NA,
                          sampleEvery = 0.01,
                          detectionSize = 6e4,
                          detectionDrivers = 99, max.wall.time = 900,
                          seed =NULL,
                          thresholdWhole = x)
    
    
    cat("\n osSFPSb: a runif is", runif(1), "\n")
    b2 <- oncoSimulSample(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE, detectionProb = NA,
                         sampleEvery = 0.01,
                          detectionSize = 6e4,
                          detectionDrivers = 99, max.wall.time = 900,
                          seed =NULL,
                         thresholdWhole = x)
    b1$popSummary[1:5, c(1:3, 8:9)]
    summary(b1$popSummary[, "NumClones"])
    summary(b1$popSummary[, "TotalPopSize"])
    b2$popSummary[1:5, c(1:3, 8:9)]
    summary(b2$popSummary[, "NumClones"])
    summary(b2$popSummary[, "TotalPopSize"])
    ## cc1 and cc2 should all be smaller than pops, or you are maxing
    ## things and not seeing patterns
    (cc1 <- colSums(b1$popSample))
    (cc2 <- colSums(b2$popSample))
    ## Of course, these are NOT really mutationsPerClone: we collapse over
    ## whole population.
    (mutsPerClone1 <- rowSums(b1$popSample))
    (mutsPerClone2 <- rowSums(b2$popSample))
    summary(mutsPerClone1)
    summary(mutsPerClone2)
    TTT <- c(TTT,  t.test(mutsPerClone2,
                          mutsPerClone1, alternative = "greater")$p.value < p.value.threshold)
    TTT <- c(TTT,  wilcox.test(b2$popSummary[, "NumClones"], 
                 b1$popSummary[, "NumClones"], alternative = "greater")$p.value < p.value.threshold)
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
cat("\n", date(), "\n")




date()
test_that("Different freqs as they should be ordered and chisq, when s  and a tiny mu", {
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n mpg s13: a runif is", runif(1), "\n")
    muvar2 <- c("U" = 1e-13, "z" = 1e-7, "e" = 1e-6, "m" = 1e-5, "D" = 1e-4)
    ni2 <- rep(0.01, 5)
    names(ni2) <- names(muvar2)
    ni2["U"] <- 0.5
    fe1 <- allFitnessEffects(noIntGenes = ni2)
    no <- 1e6
    reps <- 600
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = 5,
                       mutationPropGrowth = TRUE,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    TTT <- c(TTT, colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    TTT <- c(TTT, identical(
        order(colSums(OncoSimulR:::geneCounts(bb))),
        order(expectedC)))
    ## A chisq will not work as we increase finalTime.
    
    
    cat("\n mpg s13b: a runif is", runif(1), "\n")
    bb <- oncoSimulPop(reps,
                       fe1, mu = muvar2, onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = .001,
                       mutationPropGrowth = TRUE,
                       seed = NULL, mc.cores = 2
                       )
    (expectedC <- no*reps*muvar2)
    colSums(OncoSimulR:::geneCounts(bb))
    TTT <- c(TTT, colSums(OncoSimulR:::geneCounts(bb))[1] == 0)
    ## This will fail sometimes
    p.fail <- 1e-3
    TTT <- c(TTT, chisq.test(colSums(OncoSimulR:::geneCounts(bb))[-1],
                           p = expectedC[-1]/sum(expectedC))$p.value > p.fail)
    ## expect_equal(
    ##     order(colSums(OncoSimulR:::geneCounts(bb))),
         ##     order(expectedC))
          if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
}
)
date()




date()
test_that("Num clones, muts per clone for different per-gene-mut",{
    ## Like previous, but larger finalTime, so no longer chi-square test
    ## here.
     max.tries <- 4
     for(tries in 1:max.tries) {
         TTT <- NULL
    
    
    cat("\n mpg df2: a runif is", runif(1), "\n")
    ng <- 40
    ni <- rep(0, ng)
    m1 <- runif(ng, min = 1e-6, max = 1e-5)
    m2 <- runif(ng, min = 1e-4, max = 1e-3)
    names(ni) <- names(m1) <- names(m2) <- c(replicate(ng,
                                 paste(sample(letters, 12), collapse = "")))
    fe1 <- allFitnessEffects(noIntGenes = ni)
    ft <- 2
    no <- 1e5
    reps <- 20
    
    
    cat("\n mpg df2a: a runif is", runif(1), "\n")
    b1 <- oncoSimulPop(reps,
                       fe1,
                       mu = m1,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = TRUE, 
                       seed =NULL, max.wall.time = 900,
                       mc.cores = 2
                       )
    
    
    cat("\n mpg df2b: a runif is", runif(1), "\n")
    b2 <- oncoSimulPop(reps,
                       fe1,
                       mu = m2,
                       onlyCancer = FALSE, detectionProb = NA,
                       initSize = no,
                       finalTime = ft,
                       mutationPropGrowth = TRUE, 
                       seed =NULL, max.wall.time = 900,
                       mc.cores = 2
                       )
    TTT <- c(TTT,  wilcox.test(summary(b2)$NumClones, 
                 summary(b1)$NumClones, alternative = "greater")$p.value < p.value.threshold)
    ## Note the short time, so this is not always very different as few
    ## have double or triple mutants
    TTT <- c(TTT,  t.test(mutsPerClone(b2), 
                          mutsPerClone(b1), alternative = "greater")$p.value < p.value.threshold)
         if( all(TTT) ) break;
     }
     cat(paste("\n done tries", tries, "\n"))
     expect_true(all(TTT))
})
date()


cat("\n Ending per-gene-mutation-rates long at", date(), "\n")
