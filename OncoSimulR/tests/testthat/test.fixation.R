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
                   onlyCancer = TRUE,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = c("u", "v")
                   )
    oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = list("u", "v")
                   )
    oncoSimulPop(2, od, muEF = odm, model = "McFL",
                 mu = 1e-4, 
                 onlyCancer = TRUE,
                 initSize = initS, 
                 keepEvery = NA,
                 fixation = c("u", "v"),
                 mc.cores = 2
                 )
    oncoSimulPop(2, od, muEF = odm, model = "Exp",
                 mu = 1e-4, 
                 onlyCancer = TRUE,
                 initSize = initS, 
                 keepEvery = NA,
                 fixation = list("u", "v"),
                 mc.cores = 2
                 )
    oncoSimulSample(2, od, muEF = odm, model = "Exp",
                    mu = 1e-4, 
                    onlyCancer = TRUE,
                    initSize = 100, 
                    fixation = c("u", "v")
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
    ## these should all run

    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = c("u", "u > v")
                   ),
                 "Order effects not allowed", fixed = TRUE)

    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = c("u", "n")
                   ),
                   "The 'fixation' list contains genes that are not present",
                   fixed = TRUE)

    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = matrix(1:5, ncol = 5),
                   ),
                 "'fixation' must be a list or a vector", fixed = TRUE)

    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                   mu = 1e-4, 
                   onlyCancer = TRUE,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = 98
                   ),
                 "Each element of 'fixation' must be a single element character vector",
                 fixed = TRUE)

    
    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                                mu = 1e-4, 
                                onlyCancer = TRUE,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = list("u", c("u", "v"))
                   ),
                 "Each element of 'fixation' must be a single element character vector",
                 fixed = TRUE)

    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                                mu = 1e-4, 
                                onlyCancer = TRUE,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = list("u", list("u", "v"))
                   ),
                 "Each element of 'fixation' must be a single element character vector",
                 fixed = TRUE)

    expect_error(oncoSimulIndiv(od, muEF = odm, model = "McFL",
                                mu = 1e-4, 
                                onlyCancer = TRUE,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = c("u", "v"),
                   AND_DrvProbExit = TRUE
                   ),
                 "It makes no sense to pass AND_DrvProbExit and a fixation list",
                 fixed = TRUE)

    data(examplePosets)
    p701 <- examplePosets[["p701"]]
    expect_error(oncoSimulIndiv(p701, model = "McFL",
                                mu = 1e-4, 
                                onlyCancer = TRUE,
                   initSize = initS, 
                   keepEvery = NA,
                   fixation = c("u", "v"),
                   AND_DrvProbExit = TRUE
                   ),
                 "'fixation' cannot be specified with the old poset format",
                 fixed = TRUE)


    })



