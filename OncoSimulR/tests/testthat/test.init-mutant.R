test_that("initMutant crashes",
          {
              o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c(0.01, 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              expect_error(oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("b, a") ## also "a" and "b" separate
                      ),
                      "genes not in the fitness table")
              expect_error(oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("b", "a") ## also "a" and "b" separate
                      ),
                      "genes not in the fitness table")
              expect_error(oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("1", "2") ## also "a" and "b" separate
                      ),
                      "genes not in the fitness table")
          })


test_that("initMutant lexicog order",
          {
              o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c(0.01, 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              tmp <- oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("d > m")
                      )
              cn <- colnames(igraph::get.adjacency(plotClonePhylog(tmp, N = 0,
                                                                   returnGraph = TRUE),
                                                   sparse = FALSE))
              expect_true( "d > m _ " %in% cn )
              expect_true( "d > m > f _ " %in% cn )
              expect_false( "m > d _" %in% cn )
              expect_false( "m > d > f _" %in% cn )
          })


test_that("initMutant lexicog order with noint",
          {
              o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c("u" = 0.01, "z" = 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              tmp <- oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("d > m > z")
                      )
              cn <- colnames(igraph::get.adjacency(plotClonePhylog(tmp, N = 0,
                                                                   returnGraph = TRUE),
                                                   sparse = FALSE))
              expect_true( "d > m _ z" %in% cn )
              expect_true( "d > m > f _ z" %in% cn )
              expect_false( "m > d _" %in% cn )
              expect_false( "m > d > f _" %in% cn )
              ## expect_true( "d, m_z" %in% cn )
              ## expect_true( "d, m, f_z" %in% cn )
              ## expect_false( "m, d_" %in% cn )
              ## expect_false( "m, d, f_" %in% cn )
          })


test_that("initMutant non lexicog order",
          {
              o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c("u" = 0.01, "z" = 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              tmp <- oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("m > d")
                      )
              cn <- colnames(igraph::get.adjacency(plotClonePhylog(tmp, N = 0,
                                                                   returnGraph = TRUE),
                                                   sparse = FALSE))
              expect_false( "d > m _" %in% cn )
              expect_false( "d > m > f _" %in% cn )
              expect_true( "m > d _ " %in% cn )
              expect_true( "m > d > f _ " %in% cn )
          })


test_that("initMutant non lexicog order",
          {
              o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c("u" = 0.01, "z" = 0.01),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              tmp <- oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("m > u > d")
                      )
              cn <- colnames(igraph::get.adjacency(plotClonePhylog(tmp, N = 0,
                                                                   returnGraph = TRUE),
                                                   sparse = FALSE))
              expect_false( "d > m _" %in% cn )
              expect_false( "d > m > f _" %in% cn )
              expect_true( "m > d _ u" %in% cn )
              expect_true( "m > d > f _ u" %in% cn )
          })

## FIXME: we could use stronger test: we will never see M > D
test_that("initMutant with oncoSimulSample", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9),
                        noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = 0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d") )
    ossI <- oncoSimulSample(4, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 2,
                        onlyCancer = TRUE,
                        initSize = 500,
                        initMutant = c("z > d"),
                        thresholdWhole = 1 ## check presence of initMutant
                        )
    ssp <- ossI$popSample
    expect_equal(ssp[, c("d", "z")],
                 matrix(1, nrow = 4, ncol = 2,
                        dimnames = list(NULL, c("d", "z"))))
    expect_false(sum(ssp) == prod(dim(ssp))) ## we don't just have all of
                                             ## them, which would turn the
                                             ## previous into irrelevant
    expect_equal(grep("d",
                      as.character(ossI$popSummary$OccurringDrivers)), 1:4)
})


test_that("initMutant with oncoSimulSample, 2", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9,
                            "M > A"     = 0.25,
                            "A > H"     = 0.2,
                            "A > G"     = 0.3),
                        noIntGenes = c("u" = 0.1, 
                                       "v" = 0.2,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "A" = "a",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d",
                              "H" = "h",
                              "G" = "g") )
    ossI <- oncoSimulSample(4, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 3,
                        onlyCancer = TRUE,
                        initSize = 500,
                        initMutant = c("z > a"),
                        thresholdWhole = 1 ## check presence of initMutant
                        )
    ssp <- ossI$popSample
    expect_equal(ssp[, c("a", "z")],
                 matrix(1, nrow = 4, ncol = 2,
                        dimnames = list(NULL, c("a", "z"))))
    expect_false(sum(ssp) == prod(dim(ssp))) ## we don't just have all of
                                             ## them, which would turn the
                                             ## previous into irrelevant
    expect_equal(grep("a",
                      as.character(ossI$popSummary$OccurringDrivers)), 1:4)
})


test_that("initMutant with oncoSimulPop", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.2,
                            "M > D"     = 0.9),
                        noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d") )
    ospI <- oncoSimulPop(4, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 5000,
                        detectionDrivers = 3,
                        onlyCancer = TRUE,
                        keepPhylog = TRUE,
                        initSize = 500,
                        initMutant = c("d > m > y"),
                        mc.cores = 2
                        )
    genepos <- match(c("d", "m", "y"), ospI[[1]]$geneNames)
    expect_true( all(
        unlist(lapply(ospI,
                      function(x) x$Genotypes[genepos, , drop = FALSE])) == 1))
    ## make sure not all 1, by error, which would render previous useless
    expect_false( all( 
        unlist(lapply(ospI,
                      function(x) x$Genotypes)) == 1))
    expect_true(all(unlist(lapply(ospI,
                                  function(x) (
                                      lapply(x$GenotypesWDistinctOrderEff,
                                             function(z) genepos %in% z))))))
    ## Like before, check no silly things
    allgenepos <- seq_along(ospI[[1]]$geneNames)
    expect_false(all(unlist(lapply(ospI,
                                  function(x) (
                                      lapply(x$GenotypesWDistinctOrderEff,
                                             function(z) allgenepos %in% z))))))
    expect_true(all(
        lapply(ospI,
               function(x)
                   as.character(x$other$PhylogDF[1, 1])) == "d > m _ y"))
})



test_that("initMutant with oncoSimulPop, 2", {
    o3init <- allFitnessEffects(orderEffects = c(
                            "M > D > F" = 0.99,
                            "D > M > F" = 0.2,
                            "D > M"     = 0.1,
                            "M > D"     = 0.9),
                            noIntGenes = c("u" = 0.01, 
                                       "v" = 0.01,
                                       "w" = 0.001,
                                       "x" = 0.0001,
                                       "y" = -0.0001,
                                       "z" = -0.001),
                        geneToModule =
                            c("Root" = "Root",
                              "M" = "m",
                              "F" = "f",
                              "D" = "d") )
    ospI <- oncoSimulPop(4, 
                        o3init, model = "Exp",
                        mu = 5e-5, finalTime = 70,
                        detectionDrivers = 4, ## yes, reach end
                        onlyCancer = FALSE,
                        keepPhylog = TRUE,
                        initSize = 100,
                        initMutant = c("m > v > d"),
                        mc.cores = 2
                        )
    genepos <- match(c("d", "m", "v"), ospI[[1]]$geneNames)
    expect_true( all(
        unlist(lapply(ospI,
                      function(x) x$Genotypes[genepos, , drop = FALSE])) == 1))
    ## make sure not all 1, by error, which would render previous useless
    expect_false( all( 
        unlist(lapply(ospI,
                      function(x) x$Genotypes)) == 1))
    expect_true(all(unlist(lapply(ospI,
                                  function(x) (
                                      lapply(x$GenotypesWDistinctOrderEff,
                                             function(z) genepos %in% z))))))
    ## Like before, check no silly things
    allgenepos <- seq_along(ospI[[1]]$geneNames)
    expect_false(all(unlist(lapply(ospI,
                                  function(x) (
                                      lapply(x$GenotypesWDistinctOrderEff,
                                             function(z) allgenepos %in% z))))))
    expect_true(all(
        lapply(ospI,
               function(x)
                   as.character(x$other$PhylogDF[1, 1])) == "m > d _ v"))
##                   as.character(x$other$PhylogDF[1, 1])) == "m, d_v"))
})
