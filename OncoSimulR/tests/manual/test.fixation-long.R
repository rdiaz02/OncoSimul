## I set cores to the number of cores.
## This is not to be run by default in BioC servers, only
## on our, or user's, machines
cat("\n Starting long fixation  at", date(), "\n") 
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
    op <- oncoSimulPop(2000, od, muEF = odm, model = "McFL",
                        mu = 1e-4, 
                        onlyCancer = TRUE,
                        initSize = initS, 
                        keepEvery = NA,
                        fixation = c("u", "v"),
                        mc.cores = parallel::detectCores()
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("u", "v")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    op <- oncoSimulPop(2000, od, muEF = odm, model = "McFL",
                        mu = 1e-4, 
                        onlyCancer = TRUE,
                        initSize = initS, 
                        keepEvery = NA,
                        fixation = c("u", "i"),
                        mc.cores = parallel::detectCores()
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
    op <- oncoSimulPop(2000, od, muEF = odm, model = "McFL",
                        mu = 1e-3, 
                        onlyCancer = TRUE,
                        initSize = initS,
                        keepEvery = NA,
                        fixation = c("u"),
                        mc.cores = parallel::detectCores()
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("u")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("uv")))
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
    op <- oncoSimulPop(200, od, muEF = odm, model = "McFL",
                       mu = 1e-3, 
                       onlyCancer = TRUE,
                       initSize = initS,
                       keepEvery = NA,
                       max.num.tries = 5000,
                       fixation = c("i"),
                       mc.cores = parallel::detectCores()
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("i")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("iv")))
    ## very slow
    op <- oncoSimulPop(200, od, model = "McFL",
                       mu = 1e-3, 
                       onlyCancer = TRUE,
                       initSize = initS,
                       keepEvery = NA,
                       fixation = c("v"),
                       max.num.tries = 5000,
                       mc.cores = parallel::detectCores()
                       )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("v")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("iv")))
    initS <- 100
    u <- 0.2; i <- -0.02; vi <- 0.6; ui <- uv <- 1.2
    od <- allFitnessEffects(
        epistasis = c("u" = u,  "u:i" = ui,
                      "u:v" = uv, "i" = i,
                      "v:-i" = -Inf, "v:i" = vi))
    odm <- allMutatorEffects(noIntGenes = c("i" = 50))
    evalAllGenotypesFitAndMut(od, odm, addwt = TRUE)
    initS <- 1000
    op <- oncoSimulPop(2000, od, muEF = odm, model = "McFL",
                        mu = 1e-3, 
                        onlyCancer = TRUE,
                        initSize = initS,
                        keepEvery = NA,
                        fixation = c("u, v"),
                        mc.cores = parallel::detectCores()
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("u,v")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    op <- oncoSimulPop(2000, od, model = "McFL",
                        mu = 1e-3, 
                        onlyCancer = TRUE,
                        initSize = initS,
                        keepEvery = NA,
                        fixation = c("u,v"),
                        mc.cores = parallel::detectCores()
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("u,v")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    op <- oncoSimulPop(2000, od, model = "McFL",
                        mu = 1e-3, 
                        onlyCancer = TRUE,
                        initSize = initS,
                        keepEvery = NA,
                        fixation = c("i,v"),
                        mc.cores = parallel::detectCores()
                        )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("i,v")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("uv")))
    op <- oncoSimulPop(2000, od, model = "McFL",
                       mu = 1e-3, 
                       onlyCancer = TRUE,
                       initSize = initS,
                       keepEvery = NA,
                       fixation = c("u, i, v"),
                       mc.cores = parallel::detectCores()
                       )
    sp <- samplePop(op)
    rsop <- rowSums(sp)
    stopifnot(all(rsop >= 1))
    sg <- sampledGenotypes(sp)
    expect_true(list_g_matches_fixed(sg[, "Genotype"], c("v, i, u")))
    expect_false(list_g_matches_fixed(sg[, "Genotype"], c("uv")))
})
cat("\n Ending long fixation  at", date(), "\n") 
