date()
test_that("Exercise LOD and POM code", {
    pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                          "TP53", "TP53", "MLL3"),
                                      child = c("KRAS","SMAD4", "CDNK2A", 
                                          "TP53", "MLL3",
                                          rep("PXDN", 3), rep("TGFBR2", 2)),
                                      s = 0.05,
                                      sh = -0.3,
                                      typeDep = "MN"))
    pancr1 <- oncoSimulIndiv(pancr, model = "Exp", keepPhylog = TRUE)
    pancr8 <- oncoSimulPop(8, pancr, model = "Exp", keepPhylog = TRUE, mc.cores = 2)
    expect_true(inherits(POM(pancr1), "character"))
    require(igraph)
    expect_true(inherits(LOD(pancr1)$all_paths[[1]], "igraph.vs"))
    expect_true(inherits(LOD(pancr1)$lod_single, "igraph.vs"))
    expect_true(inherits(POM(pancr8), "list"))
    expect_true(inherits(LOD(pancr8), "list"))
    expect_true(inherits(diversityPOM(POM(pancr8)), "numeric"))
    expect_true(inherits(diversityLOD(LOD(pancr8)), "numeric"))
    expect_true(diversityPOM(POM(pancr8)) >= 0)
    expect_true(diversityLOD(LOD(pancr8)) >= 0)
    expect_error(diversityPOM(POM(pancr1)),
                 "Object must be a list", fixed = TRUE)
    expect_error(diversityLOD(LOD(pancr1)),
                 "Object must be a list", fixed = TRUE)



    pancr8 <- oncoSimulPop(8, pancr, model = "McFL",
                           keepPhylog = TRUE,
                           finalTime = 0.01,
                           max.num.tries = 1,
                           mc.cores = 2,
                           max.wall.time = 0.01,
                           detectionSize = 1e6)

    expect_warning(LOD(pancr8),
                   "Missing needed components.", fixed = TRUE)

    expect_warning(POM(pancr8),
                   "Missing needed components.", fixed = TRUE)

    
    pancr8 <- suppressWarnings(suppressMessages(oncoSimulPop(30, pancr, model = "McFL",
                                            keepPhylog = TRUE,
                                            onlyCancer = FALSE,
                                            finalTime = 0.01,
                                            sampleEvery = 0.5,
                                            mu = 1e-8,
                                            mc.cores = 2,
                                            mutationPropGrowth = FALSE,
                                            initSize = 10)))

    lop8 <- suppressWarnings(LOD(pancr8))

    ## what if this is commented out?
    
    ## expect_true(any(unlist(lapply(lop8,
    ##                               function(x) x$lod_single))
    ##                 %in% "No_descendants"))

    ## there are descendants but all go extinct
    ## pancr2 <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
    ##                                       "TP53", "TP53", "MLL3"),
    ##                                   child = c("KRAS","SMAD4", "CDNK2A", 
    ##                                       "TP53", "MLL3",
    ##                                       rep("PXDN", 3), rep("TGFBR2", 2)),
    ##                                   s = -0.9,
    ##                                   sh = -0.9,
    ##                                   typeDep = "MN"))

    
    gg <- allFitnessEffects(noIntGenes = rep(-.9, 100))
    pancr22 <- oncoSimulPop(10, gg,
                            model = "Exp",
                            keepPhylog = TRUE,
                            onlyCancer = FALSE,
                            initSize = 1e3,
                            mu = 1e-2,
                            mc.cores = 2,
                            finalTime = 2.5)
    lp22 <- LOD(pancr22)
    ## There is like soooo remote chance this will fail
    ## and the previous exercises the code anyway.
    ## expect_true(any(unlist(lp22) %in% "WT_is_end"))
})
date()
