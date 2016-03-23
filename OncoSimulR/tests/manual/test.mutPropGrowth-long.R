
## Why this does not really reflect what we want, and why number of clones
## is better that capture the idea of "more mutations". NumClones reflects
## the creation of a new clone, something that happens whenever there is a
## mutation (that does not land you on a pre-existing clone).

######################################################################
######################################################################

## ## The functions below measure number of mutated positions by summing
## ## number of alleles and dividing by pop.size. But only at final time.
## ## And large pops of clones with few muts remain and swamp.

## popS <- function(out) unlist(lapply(out, function(x) x$TotalPopSize))

## muts <- function(out) {
##     popSize <- popS(out)
##     gc <- rowSums(OncoSimulR:::geneCounts(out))
##     Muts.per.indiv <- na.omit(gc/popSize)
##     return(list(PopSize = popSize,
##                 Muts = gc,
##                 Muts.per.indiv = Muts.per.indiv,
##                 Muts.per.indiv.no.0 = Muts.per.indiv[gc > 0]))
## }


## This is better, but still you need time to allow accumulation of clones
## with many mutations, and those with few remain

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


######################################################################
######################################################################



RNGkind("L'Ecuyer-CMRG") ## for the mclapplies
## If crashes I want to see where: thus output seed.

## The tests below can occasionally fail (but that probability decreases
## as we increase number of pops), as they should.
## Same of these tests take some time. For now, show times.
cat(paste("\n Starting at", date(), "\n"))

## Beware of exiting because max number of subjects reached, and that
## could happen sooner for the faster growing, so less time to accumulate
## mutations. Likewise, differences between nca and nca2 depend on large
## enough differences in mutation, and thus a relatively large multiplier
## factor, so a large s.

cat("\n", date(), "\n")
test_that("Ordering of number of clones and mutsPerClone with mutpropgrowth, 1", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpc1: the seed is", pseed, "\n")
    ft <- 2.5
    pops <- 200
    lni <- 400
    no <- 10
    ni <- c(5, 3, rep(0, lni))
    names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
    fe <- allFitnessEffects(noIntGenes = ni)
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ncb <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no,
                         initMutant = "a",
                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no,
                         initMutant = "b",
                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(ncb)$NumClones) > 1e-4)
    expect_true(var(summary(nca2)$NumClones) > 1e-4)
    expect_true(var(summary(ncb2)$NumClones) > 1e-4)
    ## The real comparison
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(ncb)$NumClones))
    expect_true( median(summary(ncb)$NumClones) >
                 median(summary(ncb2)$NumClones))
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(ncb)))
    expect_true( mean(mutsPerClone(ncb)) >
                 mean(mutsPerClone(ncb2)))
    ## These can fail in this case, since small diffs. as small mutlipliers
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(nca2)))
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(nca2)$NumClones))
    ## In this cases, we would expect differences in total population size
    ## between a and a2, but minor or non detectable between b and b2. In
    ## a we expect them because, since those affect a lot the mutation
    ## rate, we expect them to get b faster, and thus grow faster noticeably.
})
cat("\n", date(), "\n")


cat("\n", date(), "\n")
test_that("Ordering of number of clones and mutsPerClone with mutpropgrowth, 2", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpc2: the seed is", pseed, "\n")
    ## The s coefficient is small, and so small differences between nca and
    ## nca2.
    ft <- 15 ## going beyond 16 or so, gets it to bail because of reaching max
    ## pop
    pops <- 200
    lni <- 300
    no <- 10
    ni <- c(1, 0.5, rep(0, lni))
    names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
    fe <- allFitnessEffects(noIntGenes = ni)
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ncb <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        initSize = no,
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no,
                         initMutant = "a",
                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no,
                         initMutant = "b",
                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(ncb)$NumClones) > 1e-4)
    expect_true(var(summary(nca2)$NumClones) > 1e-4)
    expect_true(var(summary(ncb2)$NumClones) > 1e-4)
    ## The real comparison
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(ncb)$NumClones))
    expect_true( median(summary(ncb)$NumClones) >
                 median(summary(ncb2)$NumClones))
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(ncb)))
    expect_true( mean(mutsPerClone(ncb)) >
                 mean(mutsPerClone(ncb2)))
    ## These can fail in this case, since small diffs. as small mutlipliers
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(nca2)))
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(nca2)$NumClones))
})



cat("\n", date(), "\n")
test_that("Ordering of number of clones and mutsPerClone with mutpropgrowth, 3", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpc3: the seed is", pseed, "\n")
    ## The s coefficient is small, and so small differences. Here, much large
    ## mu
    ft <- 12 ## going beyond 13 or so, gets it to bail because of reaching max
    ## pop in nca
    pops <- 200
    lni <- 50
    no <- 10
    ni <- c(1.0, 0.8, rep(0, lni))
    mu <- 1e-5
    names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
    fe <- allFitnessEffects(noIntGenes = ni)
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        mu = mu,
                        initSize = no,
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ncb <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        mu = mu,
                        initSize = no,
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mutationPropGrowth = FALSE,
                         mu = mu,
                         initSize = no,
                         initMutant = "a",
                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mutationPropGrowth = FALSE,
                         mu = mu,                     
                         initSize = no,
                         initMutant = "b",
                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(ncb)$NumClones) > 1e-4)
    expect_true(var(summary(nca2)$NumClones) > 1e-4)
    expect_true(var(summary(ncb2)$NumClones) > 1e-4)
    ## The real comparison
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(ncb)$NumClones))
    expect_true( median(summary(ncb)$NumClones) >
                 median(summary(ncb2)$NumClones))
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(ncb)))
    expect_true( mean(mutsPerClone(ncb)) >
                 mean(mutsPerClone(ncb2)))
    ## These can fail in this case, since small diffs. as small mutlipliers
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(nca2)))
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(nca2)$NumClones))
})


cat("\n", date(), "\n")
test_that("McFL: Ordering of number of clones and mutsPerClone with mutpropgrowth, 1", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpcmcf1: the seed is", pseed, "\n")
    ft <- 20 ## unless large you rarely get triple, etc, mutatns
    pops <- 100
    lni <- 50 
    no <- 1e3
    ni <- c(3, 1.5, rep(0, lni))
    names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
    fe <- allFitnessEffects(noIntGenes = ni)
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        initSize = no, model = "McFL",
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ncb <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        initSize = no, model = "McFL",
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no, model = "McFL",
                         initMutant = "a",
                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no, model = "McFL",
                         initMutant = "b",
                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(ncb)$NumClones) > 1e-4)
    expect_true(var(summary(nca2)$NumClones) > 1e-4)
    expect_true(var(summary(ncb2)$NumClones) > 1e-4)
    ## The real comparison
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(ncb)$NumClones))
    expect_true( median(summary(ncb)$NumClones) >
                 median(summary(ncb2)$NumClones))
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(ncb)))
    expect_true( mean(mutsPerClone(ncb)) >
                 mean(mutsPerClone(ncb2)))
    ## These can fail in this case, since small diffs. as small mutlipliers
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(nca2)))
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(nca2)$NumClones))
})
cat("\n", date(), "\n")

## A variation of the former
cat("\n", date(), "\n")
test_that("McFL: Ordering of number of clones and mutsPerClone with mutpropgrowth, 2", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpcmcf2: the seed is", pseed, "\n")
    ## Increase ft
    ft <- 50 
    pops <- 200
    lni <- 30
    no <- 1e3
    ni <- c(2, 0.5, rep(0, lni))
    names(ni) <- c("a", "b", paste0("n", seq.int(lni)))
    fe <- allFitnessEffects(noIntGenes = ni)
    nca <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        initSize = no, model = "McFL",
                        initMutant = "a",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ncb <- oncoSimulPop(pops, fe, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        initSize = no, model = "McFL",
                        initMutant = "b",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    nca2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no, model = "McFL",
                         initMutant = "a",
                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ncb2 <- oncoSimulPop(pops, fe, finalTime = ft,
                         mutationPropGrowth = FALSE,
                         initSize = no, model = "McFL",
                         initMutant = "b",
                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(ncb)$NumClones) > 1e-4)
    expect_true(var(summary(nca2)$NumClones) > 1e-4)
    expect_true(var(summary(ncb2)$NumClones) > 1e-4)
    ## The real comparison
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(ncb)$NumClones))
    expect_true( median(summary(ncb)$NumClones) >
                 median(summary(ncb2)$NumClones))
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(nca2)$NumClones))
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(ncb)))
    expect_true( mean(mutsPerClone(ncb)) >
                 mean(mutsPerClone(ncb2)))
    ## These can fail in this case, since small diffs. as small mutlipliers
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(nca2)))
})



## When diffs are very tiny in s, as we increase ft, it is easy to get a
## mutant in the second gene with non-zero s and thus the growth rates
## will be the same and there will be no differences. This is exacerbated
## with McFL as the less fit clone disappears

## We take this further below. When we init from a1 below, we start from a
## clone with smaller growth rate. But there are many b to jump
## to. However, if you start from b, you can only increase birth rate
## mutating exactly one a. Thus, over moderate finalTimes, starting a a
## will lead to more clones, etc

cat("\n", date(), "\n")
test_that("McFL: Ordering of number of clones and mutsPerClone with mutpropgrowth, and different mmodules",{
    
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n mpcmcf3: the seed is", pseed, "\n")
    ft <- 10 
    pops <- 40
    mu <- 1e-5
    lni <- 10
    fni <- 50
    no <- 1e4
    s1 <- 0.9
    s2 <- 1
    gn <- paste(paste0("b", 1:fni), collapse = ", ")
    ni <- rep(0, lni)
    names(ni) <- paste0("n", seq.int(lni))
    f1 <- allFitnessEffects(epistasis = c("A" = s1,
                                          "B" = s2),
                            geneToModule = c("A" = "a1",
                                             "B" = gn),
                            noIntGenes = ni)
    nca <- oncoSimulPop(pops, f1, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        mu = mu,
                        initSize = no, model = "McFL",
                        initMutant = "a1",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ncb <- oncoSimulPop(pops, f1, finalTime = ft,
                        mutationPropGrowth = TRUE,
                        mu = mu,
                        initSize = no, model = "McFL",
                        initMutant = "b1",
                        onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ## OK, but it also happens below just because init a1 eventually grows
    ## faster, so larger pop, so more mutants, etc
    ## We could try just counting the number of mutation events, but we would
    ## need to standardize by total time of life over all cells.
    ## As growth rates are faster (those that start with a, here, and then
    ## get a b), even when mutationPropGrowth is false, the population
    ## gets larger. And thus, it is easier to get a clone with three
    ## mutations, etc.
    nca2 <- oncoSimulPop(pops, f1, finalTime = ft,
                         mutationPropGrowth = FALSE,
                         mu = mu,
                         initSize = no, model = "McFL",
                         initMutant = "a1",
                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    ncb2 <- oncoSimulPop(pops, f1, finalTime = ft,
                         mutationPropGrowth = FALSE,
                         mu = mu,                     
                         initSize = no, model = "McFL",
                         initMutant = "b1",
                         onlyCancer = FALSE, seed = NULL, mc.cores = 2)
    expect_true(var(summary(nca)$NumClones) > 1e-4)
    expect_true(var(summary(ncb)$NumClones) > 1e-4)
    expect_true(var(summary(nca2)$NumClones) > 1e-4)
    expect_true(var(summary(ncb2)$NumClones) > 1e-4)
    ## The real comparison
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(ncb)$NumClones))
    expect_true( median(summary(ncb)$NumClones) >
                 median(summary(ncb2)$NumClones))
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(ncb)))
    expect_true( mean(mutsPerClone(ncb)) >
                 mean(mutsPerClone(ncb2)))
    ## These can fail in this case, since small diffs. as small mutlipliers
    expect_true( mean(mutsPerClone(nca)) >
                 mean(mutsPerClone(nca2)))
    expect_true( median(summary(nca)$NumClones) >
                 median(summary(nca2)$NumClones))
})



cat("\n", date(), "\n")
test_that("Without initmutant", {
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s3: the seed is", pseed, "\n")
    pops <- 40
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e3
    ft <- 5
    s3 <- 3.0
    mu <- 5e-5 ## easier to see
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f3 <- allFitnessEffects(epistasis = c("A" = s3),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    s3.ng <- oncoSimulPop(pops,
                          f3,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE,
                          seed = NULL, mc.cores = 2)
    s3.g <- oncoSimulPop(pops,
                         f3,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE,
                         seed = NULL, mc.cores = 2)
    summary(s3.g)[, c(1, 2, 3, 8, 9)]
    summary(s3.ng)[, c(1, 2, 3, 8, 9)]
    expect_true( mean(mutsPerClone(s3.g)) >
                 mean(mutsPerClone(s3.ng)))
    expect_true( median(summary(s3.g)$NumClones) >
                 median(summary(s3.ng)$NumClones))
})
cat("\n", date(), "\n")

cat("\n", date(), "\n")
test_that("Without initmutant, 2", {
    ## More of the above. Use smaller s2 and smaller mutation, but then to
    ## see it reliably you need large ft and we also increase
    ## init. pop. size.
    pseed <- sample(9999999, 1)
    set.seed(pseed)
    cat("\n s2: the seed is", pseed, "\n")
    s2 <- 1.0
    ft <- 15
    pops <- 40
    lni <- 1 ## no fitness effects genes
    fni <- 50 ## fitness effects genes
    no <- 1e4
    mu <- 5e-6 ## easier to see
    ## noInt have no fitness effects, but can accumulate mutations
    ni <- rep(0, lni)
    ## Those with fitness effects in one module, so
    ## neither fitness nor mut. rate blow up
    gn <- paste(paste0("a", 1:fni), collapse = ", ")
    f2 <- allFitnessEffects(epistasis = c("A" = s2),
                            geneToModule = c("A" = gn),
                            noIntGenes = ni)
    s2.ng <- oncoSimulPop(pops,
                          f2,
                          mu = mu,
                          mutationPropGrowth = FALSE,
                          finalTime =ft,
                          initSize = no,
                          onlyCancer = FALSE,
                          seed = NULL, mc.cores = 2)
    s2.g <- oncoSimulPop(pops,
                         f2,
                         mu = mu,
                         mutationPropGrowth = TRUE,
                         finalTime =ft,
                         initSize = no,
                         onlyCancer = FALSE,
                         seed = NULL, mc.cores = 2)
    summary(s2.g)[, c(1, 2, 3, 8, 9)]
    summary(s2.ng)[, c(1, 2, 3, 8, 9)]
    expect_true( mean(mutsPerClone(s2.g)) >
                 mean(mutsPerClone(s2.ng)))
    expect_true( median(summary(s2.g)$NumClones) >
                 median(summary(s2.ng)$NumClones))
})
cat("\n", date(), "\n")



## ##     A way to check is to see the output from the C++ code with the
## ##     verbosity option.

## set.seed("Mersenne-Twister")
## ni <- rep(0.4, 20)
## names(ni) <- c("a", "b", "c", "d", paste0("n", 1:16))
## fe <- allFitnessEffects(noIntGenes = ni)
## set.seed(5) 
## oncoSimulIndiv(fe, finalTime =30,
##                mutationPropGrowth = TRUE,
##                initSize = 1e4,
##                mu = 1e-06,
##                verbosity = 6,
##                onlyCancer = FALSE)

## ###### Iteration 30.
## ## mutation
## ## child
## 1.4 * 1e-06 * 19

## ni <- rep(0.4, 20)
## names(ni) <- c("a", "b", "c", "d", paste0("n", 1:16))
## fe <- allFitnessEffects(noIntGenes = ni)
## set.seed(25) 
## oncoSimulIndiv(fe, finalTime =40,
##                mutationPropGrowth = TRUE,
##                initSize = 1e4,
##                mu = 1e-06,
##                verbosity = 6,
##                onlyCancer = FALSE)
## ## Iteration 48.
## ## Birth of child:
## 1.4 * 1.4
## ## Mutation of child
## 1.96 * 1e-06 * 18

cat(paste("\n Ending at", date(), "\n"))
