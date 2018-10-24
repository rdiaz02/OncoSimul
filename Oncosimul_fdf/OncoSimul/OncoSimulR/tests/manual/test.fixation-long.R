## This is not to be run by default in BioC servers, only
## on our, or user's, machines
## But using all cores might be bad, and also precludes my logic
## of launching many at the same time.
cat("\n Starting long fixation  at", date(), "\n")

## Since we run the long tests after simply loading the package
## list_g_matches_fixed <- OncoSimulR:::list_g_matches_fixed

test_that("Check output is correct", {
    initS <- 200
    u <- 0.2; i <- -0.02; vi <- 0.6; ui <- uv <- -Inf
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:-i" = -Inf, "v:i" = vi))
    ## drvNames = c("u", "v"))
    odm <- allMutatorEffects(noIntGenes = c("i" = 50))
    evalAllGenotypesFitAndMut(od, odm, addwt = TRUE)
    initS <- 2000
    cat("\n\n ************** fixation 1:  ***********\n")
    op <- oncoSimulPop(2000, od, muEF = odm, model = "McFL",
                        mu = 1e-4, 
                        onlyCancer = TRUE, finalTime = 5000, detectionSize = NA, detectionProb = NA,
                        initSize = initS, 
                        keepEvery = NA,
                        fixation = c("u", "v"),
                        mc.cores = 2
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("u", "v")))
    expect_false(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    cat("\n\n ************** fixation 2:  ***********\n")
    op <- oncoSimulPop(2000, od, muEF = odm, model = "McFL",
                        mu = 1e-4, 
                        onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                        initSize = initS, 
                        keepEvery = NA,
                        fixation = c("u", "i"),
                        mc.cores = 2
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("u", "i")))
    expect_false(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("ui")))
    ## this takes much longer, of course, so only 100
    ## and increase mu and its fitness
    u <- 0.8; i <- -0.02; vi <- 0.6; ui <- uv <- -Inf
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:-i" = -Inf, "v:i" = vi))
    ## drvNames = c("u", "v"))
    odm <- allMutatorEffects(noIntGenes = c("i" = 50))
    evalAllGenotypesFitAndMut(od, odm, addwt = TRUE)
    initS <- 2000
    cat("\n\n ************** fixation 3:  ***********\n")    
    op <- oncoSimulPop(2000, od, muEF = odm, model = "McFL",
                        mu = 1e-3, 
                        onlyCancer = TRUE, finalTime = 5000, detectionSize = NA, detectionProb = NA,
                        initSize = initS,
                        keepEvery = NA,
                        fixation = c("u"),
                        mc.cores = 2
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("u")))
    expect_false(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    ## this takes much longer, of course, 
    ## and increase mu and its fitness
    u <- 0.8; i <- -0.02; vi <- 2.6; ui <- uv <- -Inf
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:-i" = -Inf, "v:i" = vi))
    ## drvNames = c("u", "v"))
    odm <- allMutatorEffects(noIntGenes = c("i" = 50))
    evalAllGenotypesFitAndMut(od, odm, addwt = TRUE)
    initS <- 50
    cat("\n\n ************** fixation 4:  ***********\n")    
    op <- oncoSimulPop(200, od, muEF = odm, model = "McFL",
                       mu = 1e-3, 
                       onlyCancer = TRUE, finalTime = 1000, detectionSize = NA, detectionProb = NA,
                       initSize = initS,
                       keepEvery = NA,
                       max.num.tries = 5000,
                       fixation = c("i"),
                       mc.cores = 2
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("i")))
    expect_false(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("iv")))
    ## very slow
    cat("\n\n ************** fixation 5:  ***********\n")    
    op <- oncoSimulPop(200, od, model = "McFL",
                       mu = 1e-3, 
                       onlyCancer = TRUE, finalTime = 1000, detectionSize = NA, detectionProb = NA,
                       initSize = 50,
                       keepEvery = NA,
                       fixation = c("v"),
                       max.num.tries = 10000,
                       mc.cores = 2
                       )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("v")))
    expect_false(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("iv")))
    initS <- 100
    u <- 0.2; i <- -0.02; vi <- 0.6; ui <- uv <- 1.2
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:i" = vi))
    odm <- allMutatorEffects(noIntGenes = c("i" = 50))
    evalAllGenotypesFitAndMut(od, odm, addwt = TRUE)
    initS <- 1000
    cat("\n\n ************** fixation 6:  ***********\n")    
    op <- oncoSimulPop(2000, od, muEF = odm, model = "McFL",
                        mu = 1e-3, 
                        onlyCancer = TRUE, finalTime = 5000, detectionSize = NA, detectionProb = NA,
                        initSize = initS,
                        keepEvery = NA,
                        fixation = c("u, v"),
                        mc.cores = 2
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("u,v")))
    expect_false(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    cat("\n\n ************** fixation 7:  ***********\n")    
    op <- oncoSimulPop(1000, od, model = "McFL",
                        mu = 1e-3, 
                        onlyCancer = TRUE, finalTime = 1000, detectionSize = NA, detectionProb = NA,
                        initSize = initS,
                        keepEvery = NA,
                        fixation = c("u,v"),
                        mc.cores = 2
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("u,v")))
    expect_false(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    cat("\n\n ************** fixation 8:  ***********\n")    
    op <- oncoSimulPop(200, od, model = "McFL",
                        mu = 1e-3, 
                        onlyCancer = TRUE, finalTime = 1000, detectionSize = NA, detectionProb = NA,
                        initSize = 30,
                        keepEvery = NA,
                        fixation = c("i,v"),
                        mc.cores = 2
                       )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("i,v")))
    expect_false(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("uv")))

    
    cat("\n\n ************** fixation 9:  ***********\n")    
    op <- oncoSimulPop(2000, od, model = "McFL",
                       mu = 1e-3, 
                       onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                       initSize = initS,
                       keepEvery = NA,
                       fixation = c("u, i, v"),
                       mc.cores = 2
                       )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("v, i, u")))
    expect_false(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("uv")))


    
    u <- 0.2; i <- -0.02; vi <- 1.6; ui <- uv <- -Inf
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:i" = vi))
    cat("\n\n ************** fixation 10:  ***********\n")    
    op <- oncoSimulPop(2000, od, model = "McFL",
                       mu = 1e-3, 
                       onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                       initSize = initS,
                       keepEvery = NA,
                       fixation = c("u", "i", "v"),
                       mc.cores = 2
                       )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("v", "i", "u")))
    expect_false(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("uv")))


    u <- 0.2; i <- -0.02; vi <- 1.6; 
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = -.9, "v" = 1.2,
                      "u:v" = 2, "i" = i,
                      "v:i" = vi, "u:i:v" = 10))
    evalAllGenotypes(od, addwt = TRUE)
    cat("\n\n ************** fixation 11:  ***********\n")    
    op <- oncoSimulPop(2000, od, model = "McFL",
                       mu = 1e-3, 
                       onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                       initSize = 100,
                       keepEvery = NA,
                       fixation = c("u,i,v"),
                       mc.cores = 2
                       )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("v, i, u")))
    expect_false(OncoSimulR:::list_g_matches_fixed(sg[, "Genotype"], c("uv")))



    
})




## the really relevant is r2. The other two if they succeed in few
## will succeed in more.
system.time(
    test_that("long Local max: not stopping, stopping, and tolerance", {
         cat("\n\n ************** fixation: long local max  ***********\n")  
        initS <- 4000
    r1 <- matrix(0, ncol = 6, nrow = 9)
    colnames(r1) <- c(LETTERS[1:5], "Fitness")
    r1[1, 6] <- 1
    r1[cbind((2:4), c(1:3))] <- 1
    r1[2, 6] <- 1.4
    r1[3, 6] <- 1.32
    r1[4, 6] <- 1.32
    r1[5, ] <- c(0, 1, 0, 0, 1, 1.5)
    r1[6, ] <- c(0, 0, 1, 1, 0, 1.54)
    r1[7, ] <- c(1, 0, 1, 1, 0, 1.65)
    r1[8, ] <- c(1, 1, 1, 1, 0, 1.75)
    r1[9, ] <- c(1, 1, 1, 1, 1, 1.85)
    class(r1) <- c("matrix", "genotype_fitness_matrix")
    ## plot(r1) ## to see the fitness landscape
    local_max_g <- c("A", "B, E", "A, B, C, D, E")
    local_max <- paste0("_,", local_max_g)
    fr1 <- allFitnessEffects(genotFitness = r1, drvNames = LETTERS[1:5])
    ## we pass sets of genes, so not stopping at genotypes
    r1 <- oncoSimulPop(200,
                       fp = fr1,
                       model = "McFL",
                       initSize = initS,
                       mu = 1e-4,
                       detectionSize = NA,
                       sampleEvery = .03,
                       keepEvery = 1, 
                       finalTime = 50000,
                       fixation = local_max_g, 
                       detectionDrivers = NA,
                       detectionProb = NA,
                       onlyCancer = TRUE,
                       max.num.tries = 500,
                       max.wall.time = 20, 
                       errorHitMaxTries = TRUE,
                       keepPhylog = FALSE,
                       mc.cores = 2)
    sr1 <- summary(r1)
    expect_true(!all(sr1$TotalPopSize == sr1$LargestClone))
    sp1 <- samplePop(r1, "last", "singleCell")
    sgsp1 <- sampledGenotypes(sp1)
    expect_true(!all(sgsp1$Genotype %in% local_max_g))
    ## genotypes, exactly
    mm <- rep(1e-4, 5)
    mm[3] <- 2e-4
    mm[1] <- 1e-5
    names(mm) <- LETTERS[1:5]
    r2 <- oncoSimulPop(2000,
                       fp = fr1,
                       model = "McFL",
                       initSize = initS,
                       mu = 1e-4,
                       detectionSize = NA,
                       sampleEvery = .03,
                       keepEvery = 1, 
                       finalTime = 50000,
                       fixation = local_max, 
                       detectionDrivers = NA,
                       detectionProb = NA,
                       onlyCancer = TRUE,
                       max.num.tries = 500,
                       max.wall.time = 20, 
                       errorHitMaxTries = TRUE,
                       keepPhylog = FALSE,
                       mc.cores = 2)
    sr2 <- summary(r2)
    expect_true(all(sr2$TotalPopSize == sr2$LargestClone))
    sp2 <- samplePop(r2, "last", "singleCell")
    sgsp2 <- sampledGenotypes(sp2)
    expect_true(all(sgsp2$Genotype %in% local_max_g))
    ## tolerance
    r3 <- oncoSimulPop(200,
                       fp = fr1,
                       model = "McFL",
                       initSize = initS,
                       mu = 1e-4,
                       detectionSize = NA,
                       sampleEvery = .03,
                       keepEvery = 1, 
                       finalTime = 50000,
                       fixation = c(local_max, fixation_tolerance = 0.1),
                       detectionDrivers = NA,
                       detectionProb = NA,
                       onlyCancer = TRUE,
                       max.num.tries = 500,
                       max.wall.time = 20, 
                       errorHitMaxTries = TRUE,
                       keepPhylog = FALSE,
                       mc.cores = 2)
    sr3 <- summary(r3)
    expect_true(!all(sr3$TotalPopSize == sr3$LargestClone))
    sp3 <- samplePop(r3, "last", "singleCell")
    sgsp3 <- sampledGenotypes(sp3)
    expect_true(!all(sgsp3$Genotype %in% local_max_g))
})
)


cat("\n Ending long fixation  at", date(), "\n") 

