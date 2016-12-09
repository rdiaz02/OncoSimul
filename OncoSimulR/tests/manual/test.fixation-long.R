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
cat("\n Ending long fixation  at", date(), "\n") 

