inittime <- Sys.time()
cat("\n Starting fixation  at", date(), "\n") ## about 4 seconds

test_that("Minimal run", {
    initS <- 2000
    u <- 0.2; i <- -0.02; vi <- 0.6; ui <- uv <- -Inf
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:-i" = -Inf, "v:i" = vi))
    ## drvNames = c("u", "v"))
    odm <- allMutatorEffects(noIntGenes = c("i" = 50))
    ## evalAllGenotypesFitAndMut(od, odm, addwt = TRUE)
    initS <- 1000
    ## these should all run
    oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE, finalTime = 15000, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = c("u", "v")
                   )
    oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE, finalTime = 15000, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = list("u", "v")
                   )
    oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE, finalTime = 15000, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = list("u", "v", fixation_tolerance = 0.01)
                   )
    oncoSimulPop(2, od, muEF = odm, model = "McFL",
                 mu = 1e-4, 
                 onlyCancer = TRUE, finalTime = 15000, detectionSize = NA, detectionProb = NA,
                 initSize = initS, 
                 keepEvery = NA,
                 fixation = c("u", "v"),
                 mc.cores = 2
                 )
    oncoSimulPop(2, od, muEF = odm, model = "McFL",
                 mu = 1e-4, 
                 onlyCancer = TRUE, finalTime = 15000, detectionSize = NA, detectionProb = NA,
                 initSize = initS, 
                 keepEvery = NA,
                 fixation = list("u", "v"),
                 mc.cores = 2
                 )
    oncoSimulSample(2, od, muEF = odm, model = "McFL",
                    mu = 1e-4, 
                    onlyCancer = TRUE, finalTime = 15000, detectionSize = NA, detectionProb = NA,
                    initSize = 100, 
                    fixation = c("u", "v")
                    )
    ## A gene without fitness effects
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:-i" = -Inf, "v:i" = vi),
        noIntGenes = c("m" = 0))
    odm <- allMutatorEffects(noIntGenes = c("i" = 50, "m" = 5))
    oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE, finalTime = 15000, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = c("u", "v", "m")
                   )
})

test_that("Catch errors", {
    initS <- 2000
    u <- 0.2; i <- -0.02; vi <- 0.6; ui <- uv <- -Inf
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:-i" = -Inf, "v:i" = vi))
    ## drvNames = c("u", "v"))
    odm <- allMutatorEffects(noIntGenes = c("i" = 50))
    ## evalAllGenotypesFitAndMut(od, odm, addwt = TRUE)
    initS <- 1000
    ## these should all fail
    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = c("u", "u > v")
                   ),
                 "Order effects not allowed", fixed = TRUE)
    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = c("u", "n")
                   ),
                   "The 'fixation' list contains genes that are not present",
                   fixed = TRUE)
    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = matrix(1:5, ncol = 5),
                   ),
                 "'fixation' must be a list or a vector", fixed = TRUE)
    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = 98
                   ),
                 "At least one element of 'fixation' must be a single element character vector",
                 fixed = TRUE)
    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                                mu = 1e-4, 
                                onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = list("u", c("u", "v"))
                   ),
                 "Each element of 'fixation' must be a single element vector",
                 fixed = TRUE)
    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                                mu = 1e-4, 
                                onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = list("u", list("u", "v"))
                   ),
                 "Each element of 'fixation' must be a single element vector",
                 fixed = TRUE)
    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                                mu = 1e-4, 
                                onlyCancer = TRUE, finalTime = 15000, detectionSize = NA, 
                                initSize = initS,
                                keepEvery = NA,
                                fixation = c("u", "v"),
                                detectionDrivers = 1,
                                detectionProb = "default",
                                AND_DrvProbExit = TRUE
                   ),
                 "It makes no sense to pass AND_DrvProbExit and a fixation list",
                 fixed = TRUE)
    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    expect_error(oncoSimulIndiv(p701, model = "McFL",
                                mu = 1e-4, 
                                onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                                initSize = initS, 
                                keepEvery = NA,
                                fixation = c("u", "v")
                                ),
                 "'fixation' cannot be specified with the old poset format",
                 fixed = TRUE)
    expect_error(oncoSimulIndiv(od, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = c("u",fixation_tolerance = -0.3)
                   ),
                 "Impossible range for fixation tolerance", fixed = TRUE)
    expect_error(oncoSimulIndiv(od,  model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = c("u",fixation_tolerance = 1.1)
                   ),
                 "Impossible range for fixation tolerance", fixed = TRUE)
})

test_that("Internal checking function works", {
    x <- c("i, u, g", "g, u", "i")
    y <- c("i", "g") ## TRUE
    y2 <- c("g") ## FALSE
    y3 <- c("g, i", "i") ## FALSE
    y4 <- c("i, u", "i, g", "u, g", "u") ## FALSE
    y5 <- c("i", "u, g") ## TRUE
    y6 <- c("i, u", "i, g", "u") ## FALSE
    y7 <- c("u", "g", "i") ## TRUE
    y8 <- c("g", "i") ## FALSE
    y9 <- c("i, u", "i, g", "u, g", "i") ## TRUE
    y10 <- c("i, u", "i, g", "u") ## FALSE
    y11 <- c("g", "u") ## FALSE
    y12 <- c("i, u", "i, g", "u", "i") ## TRUE
    expect_true(list_g_matches_fixed(x, y))
    expect_false(list_g_matches_fixed(x, y2))
    expect_false(list_g_matches_fixed(x, y3)) 
    expect_false(list_g_matches_fixed(x, y4))
    expect_true(list_g_matches_fixed(x, y5))
    expect_false(list_g_matches_fixed(x, y6))
    expect_true(list_g_matches_fixed(x, y7))
    expect_true(list_g_matches_fixed(x, y8))  
    expect_true(list_g_matches_fixed(x, y9))
    expect_false(list_g_matches_fixed(x, y10))
    expect_false(list_g_matches_fixed(x, y11))
    expect_true(list_g_matches_fixed(x, y12))    
    x2 <- c("i, u, g", "g, u", "i, u")
    y <- c("i", "g") ## TRUE
    y2 <- c("g") ## FALSE
    y3 <- c("g, i", "i") ## FALSE
    y4 <- c("i, u", "i, g", "u, g", "u") ## TRUE
    y5 <- c("i", "u, g") ## TRUE
    y6 <- c("i, u", "i, g", "u") ## TRUE
    y7 <- c("u", "g", "i") ## TRUE
    y8 <- c("g", "i") ## TRUE
    y9 <- c("i, u", "i, g", "u, g", "i") ## TRUE
    y10 <- c("i, u", "i, g", "u") ## TRUE
    y11 <- c("g", "u") ## TRUE
    y12 <- c("i, u", "i, g", "u", "i") ## TRUE
    expect_true(list_g_matches_fixed(x2, y))
    expect_false(list_g_matches_fixed(x2, y2))
    expect_false(list_g_matches_fixed(x2, y3)) 
    expect_true(list_g_matches_fixed(x2, y4))
    expect_true(list_g_matches_fixed(x2, y5))
    expect_true(list_g_matches_fixed(x2, y6))
    expect_true(list_g_matches_fixed(x2, y7))
    expect_true(list_g_matches_fixed(x2, y8))  
    expect_true(list_g_matches_fixed(x2, y9))
    expect_true(list_g_matches_fixed(x2, y10))
    expect_true(list_g_matches_fixed(x2, y11))
    expect_true(list_g_matches_fixed(x2, y12))    
    x3 <- c("i, u", "g")
    y <- c("i", "g") ## TRUE
    y2 <- c("g") ## FALSE
    y3 <- c("g, i", "i") ## FALSE
    y4 <- c("i, u", "i, g", "u, g", "u") ## FALSE
    y5 <- c("i", "u, g") ## FALSE
    y6 <- c("i, u", "i, g", "u") ## FALSE
    y7 <- c("u", "g", "i") ## TRUE
    y8 <- c("g", "i") ## TRUE
    y9 <- c("i, u", "i, g", "u, g", "i") ## FALSE
    y10 <- c("i, u", "i, g", "u") ## FALSE
    y11 <- c("g", "u") ## TRUE
    y12 <- c("i, u", "i, g", "u", "i") ## FALSE
    expect_true(list_g_matches_fixed(x3, y))
    expect_false(list_g_matches_fixed(x3, y2))
    expect_false(list_g_matches_fixed(x3, y3)) 
    expect_false(list_g_matches_fixed(x3, y4))
    expect_false(list_g_matches_fixed(x3, y5))
    expect_false(list_g_matches_fixed(x3, y6))
    expect_true(list_g_matches_fixed(x3, y7))
    expect_true(list_g_matches_fixed(x3, y8))  
    expect_false(list_g_matches_fixed(x3, y9))
    expect_false(list_g_matches_fixed(x3, y10))
    expect_true(list_g_matches_fixed(x3, y11))
    expect_false(list_g_matches_fixed(x3, y12))    
})

date()
test_that("Check output is correct", {
    initS <- 2000
    u <- 0.2; i <- -0.02; vi <- 0.6; ui <- uv <- -Inf
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:-i" = -Inf, "v:i" = vi))
    ## drvNames = c("u", "v"))
    odm <- allMutatorEffects(noIntGenes = c("i" = 50))
    evalAllGenotypesFitAndMut(od, odm, addwt = TRUE)
    initS <- 2000
    op <- oncoSimulPop(100, od, muEF = odm, model = "McFL",
                        mu = 1e-4, 
                       onlyCancer = TRUE, finalTime = 15000,
                       detectionSize = NA, detectionProb = NA,
                        initSize = initS, 
                        keepEvery = NA,
                        fixation = c("u", "v"),
                        mc.cores = 2
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("u", "v")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    op <- oncoSimulPop(100, od, muEF = odm, model = "McFL",
                        mu = 1e-4, 
                       onlyCancer = TRUE, finalTime = NA,
                       detectionSize = NA, detectionProb = NA,
                        initSize = initS, 
                        keepEvery = NA,
                        fixation = c("u", "i"),
                        mc.cores = 2
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("u", "i")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("ui")))
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
    op <- oncoSimulPop(100, od, muEF = odm, model = "McFL",
                        mu = 1e-3, 
                        onlyCancer = TRUE, finalTime = 15000, detectionSize = NA, detectionProb = NA,
                        initSize = initS,
                        keepEvery = NA,
                        fixation = c("u"),
                        mc.cores = 2
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("u")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    ## ## very slow; in long tests; just a minimal thing
    ## ## this takes much longer, of course, 
    ## ## and increase mu and its fitness
    u <- 0.2; i <- 0.1; vi <- 2.6; ui <- uv <- -Inf
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:-i" = -Inf, "v:i" = vi))
    ## drvNames = c("u", "v"))
    odm <- allMutatorEffects(noIntGenes = c("i" = 50))
    evalAllGenotypesFitAndMut(od, odm, addwt = TRUE)
    initS <- 20
    op <- oncoSimulPop(10, ## 50
                       od, model = "McFL",
                       mu = 1e-3, 
                       onlyCancer = TRUE,
                       finalTime = 1000, detectionSize = NA, detectionProb = NA,
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
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("i")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("iv")))
    initS <- 2000
    u <- 0.2; i <- 0; vi <- 0.6; ui <- uv <- 1.2
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:-i" = -Inf, "v:i" = vi))
    odm <- allMutatorEffects(noIntGenes = c("i" = 50))
    evalAllGenotypesFitAndMut(od, odm, addwt = TRUE)
    initS <- 2000
    op <- oncoSimulPop(100, od, muEF = odm, model = "McFL",
                        mu = 1e-3, 
                        onlyCancer = TRUE, finalTime = 15000, detectionSize = NA, detectionProb = NA,
                        initSize = initS,
                        keepEvery = NA,
                        fixation = c("u, v"),
                        mc.cores = 2
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("u,v")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    op <- oncoSimulPop(100, od, model = "McFL",
                        mu = 1e-3, 
                        onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                        initSize = initS,
                        keepEvery = NA,
                        fixation = c("u,v"),
                        mc.cores = 2
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("u,v")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    op <- oncoSimulPop(100, od, model = "McFL",
                        mu = 1e-3, 
                        onlyCancer = TRUE, finalTime = 1500, detectionSize = NA, detectionProb = NA,
                        initSize = initS,
                        keepEvery = NA,
                        fixation = c("i,v"),
                        mc.cores = 2
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("i,v")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    op <- oncoSimulPop(100, od, model = "McFL",
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
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("v, i, u")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("uv")))
})
date()


## ## This need not ever finish successfully, as there is
## ## no competition in the Exp model and we might need to
## ## run maaaaany times to filter for that.
## op5 <- oncoSimulPop(100, od, model = "Exp",
##                     mu = 1e-3, 
##                     onlyCancer = TRUE,
##                     initSize = 50,
##                     keepEvery = NA,
##                     fixation = c("v"),
##                     mc.cores = 2
##                     )
## sp5 <- samplePop(op5)
## rsop5 <- rowSums(sp5)
## stopifnot(all(rsop5 >= 1))
## sg5 <- sampledGenotypes(sp5)
## expect_true(list_g_matches_fixed(sg5[, "Genotype"], c("v")))




test_that("Local max: not stopping, stopping, and tolerance", {
        initS <- 2000
    rl1 <- matrix(0, ncol = 6, nrow = 9)
    colnames(rl1) <- c(LETTERS[1:5], "Fitness")
    rl1[1, 6] <- 1
    rl1[cbind((2:4), c(1:3))] <- 1
    rl1[2, 6] <- 1.4
    rl1[3, 6] <- 1.32
    rl1[4, 6] <- 1.32
    rl1[5, ] <- c(0, 1, 0, 0, 1, 1.5)
    rl1[6, ] <- c(0, 0, 1, 1, 0, 1.54)
    rl1[7, ] <- c(1, 0, 1, 1, 0, 1.65)
    rl1[8, ] <- c(1, 1, 1, 1, 0, 1.75)
    rl1[9, ] <- c(1, 1, 1, 1, 1, 1.85)
    class(rl1) <- c("matrix", "genotype_fitness_matrix")
    ## plot(rl1) ## to see the fitness landscape
    local_max_g <- c("A", "B, E", "A, B, C, D, E")
    local_max <- paste0("_,", local_max_g)
    fr1 <- allFitnessEffects(genotFitness = rl1, drvNames = LETTERS[1:5])
    ## we pass sets of genes, so not stopping at genotypes
    r1 <- oncoSimulPop(100,
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
    r2 <- oncoSimulPop(100,
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
    ## yes, this can occasionally fail, because all are in
    ## the list of local_max_g
        
    r3 <- oncoSimulPop(200,
                       fp = fr1,
                       model = "McFL",
                       initSize = initS,
                       mu = 1e-4,
                       detectionSize = NA,
                       sampleEvery = .03,
                       keepEvery = 1, 
                       finalTime = 50000,
                       fixation = c(local_max, fixation_tolerance = 0.5),
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
    ## sum(sgsp3$Genotype %in% local_max_g)
    ## sum(!(sgsp3$Genotype %in% local_max_g))
})


cat("\n Ending fixation  at", date(), "\n")
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
