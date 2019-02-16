cat(paste("\n Starting sample-only-last tests", date(), "\n"))

data(examplePosets)
RNGkind("Mersenne-Twister") 

bozic <- function(poset) oncoSimulIndiv(poset, sampleEvery = 0.03,
                                        keepEvery = 1)
bozic9 <- function(poset) oncoSimulIndiv(poset, sampleEvery = 0.03,
                                         keepEvery = -9)
    
Exp <- function(poset) oncoSimulIndiv(poset, sampleEvery = 0.03, keepEvery = 1,
                                      model = "Exp")
Exp9 <- function(poset) oncoSimulIndiv(poset, model = "Exp",
                                       sampleEvery = 0.03,
                                       keepEvery = -9)

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


## A list, where each element is a list
## with two members, both expected to give same result
FexpectedSame <- list(
    Bozics = list(bozic = bozic, bozic9 = bozic9),
    Exps = list(Exp = Exp, Exp9 = Exp9),
    McFLs = list(mc = mc, mc9 = mc9)
)


popsNoZero <- function(x) {
    keep <- which(x$pops.by.time[nrow(x$pops.by.time), ] > 0)
    pops <- x$pops.by.time[nrow(x$pops.by.time), keep]
    genots <- x$Genotypes[, keep[-1] - 1, drop = FALSE] ## keep contains time too in first col.
    return(list(pops = pops, genots = genots))
}

runBothFuncts <- function(seed, Poset, functionPair) {
    set.seed(seed); b1 <- functionPair[[1]](Poset)
    set.seed(seed); b2 <- functionPair[[2]](Poset)
    return(list(all = b1, last = b2))
}

## We only use one from each of 11, 9, 7.
## Whole collection tested in long tests
examplePosets <- examplePosets[c(1, 5, 9)]
for(i in 1:length(examplePosets)) {
    s1 <- round(runif(1) * 10000) ## do better as: as.integer(runif(1, 1, 1e9))
    Poset <- examplePosets[[i]]
    attributes(Poset)$namePoset <- names(examplePosets)[[i]]
    for(ffs in FexpectedSame) {
        test_that(paste("Sampling only last same for ",
                        paste(names(ffs), collapse = " ")), {
                            ## comment next cat later
                            cat(paste("\n ",
                                      " Seed = ", s1, " ",
                                      paste(names(ffs), collapse = " "),
                                      ". Poset = ",
                                      attributes(Poset)$namePoset,
                                      "\n"))
                            
                            bb <- runBothFuncts(s1, Poset, ffs)
                            b1 <- bb$all
                            b2 <- bb$last
                            popsGenots <- popsNoZero(b1)
                            expect_equal(b1$TotalPopSize, b2$TotalPopSize)
                            expect_equal(b1$FinalTime, b2$FinalTime)
                            expect_equal(b1$NumIter, b2$NumIter)

                            expect_equal(b1$NumDriversLargestPop,
                                         b2$NumDriversLargestPop)
                            expect_equal(b1$MaxDriversLast, b2$MaxDriversLast)
                            expect_equal(b1$PropLargestPopLast,
                                         b2$PropLargestPopLast)                            
                            expect_equal(b1$LargestClone, b2$LargestClone)
                            ## these need not be the same as those
                            ## accumulate over all samples

                            ## expect_equal(b1$MaxNumDrivers,
                            ## b2$MaxNumDrivers)
                            ## expect_equal(b1$TotalPresentDrivers,
                            ## b2$TotalPresentDrivers)

                            expect_equal(popsGenots$pops, b2$pops.by.time[1, ])
                            expect_equal(popsGenots$genots, b2$Genotypes)
                            expect_false(
                                all(dim(b2$pops.by.time) ==
                                        dim(b1$pops.by.time)))
                            rm(bb, b1, b2)
                        })   
    }
}

set.seed(NULL)
cat(paste("\n Ending sample-only-last tests", date(), "\n"))
