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
                 "Each element of 'fixation' must be a single element character vector",
                 fixed = TRUE)
    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                                mu = 1e-4, 
                                onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = list("u", c("u", "v"))
                   ),
                 "Each element of 'fixation' must be a single element character vector",
                 fixed = TRUE)
    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                                mu = 1e-4, 
                                onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = list("u", list("u", "v"))
                   ),
                 "Each element of 'fixation' must be a single element character vector",
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
                        onlyCancer = TRUE, finalTime = 15000, detectionSize = NA, detectionProb = NA,
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
    op <- oncoSimulPop(50, od, model = "McFL",
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
    ## ## very slow; in long tests
    ## op <- oncoSimulPop(20, od, model = "McFL",
    ##                    mu = 1e-3, 
    ##                    onlyCancer = TRUE, finalTime = NA, detectionSize = NA, detectionProb = NA,
    ##                    initSize = initS,
    ##                    keepEvery = NA,
    ##                    fixation = c("v"),
    ##                    max.num.tries = 5000,
    ##                    mc.cores = 2
    ##                    )
    ## sp <- samplePop(op)
    ## rsop <- rowSums(sp)
    ## stopifnot(all(rsop >= 1))
    ## sg <- sampledGenotypes(sp)
    ## expect_true(list_g_matches_fixed(sg[, "Genotype"], c("v")))
    ## expect_false(list_g_matches_fixed(sg[, "Genotype"], c("iv")))
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

cat("\n Ending fixation  at", date(), "\n")
