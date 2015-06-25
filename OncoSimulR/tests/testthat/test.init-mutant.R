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
                                          ))
              expect_error(oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("b", "a") ## also "a" and "b" separate
                                          ))
              expect_error(oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("1", "2") ## also "a" and "b" separate
                                          ))
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
              expect_true( "d, m_" %in% cn )
              expect_true( "d, m, f_" %in% cn )
              expect_false( "m, d_" %in% cn )
              expect_false( "m, d, f_" %in% cn )
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
              expect_true( "d, m_z" %in% cn )
              expect_true( "d, m, f_z" %in% cn )
              expect_false( "m, d_" %in% cn )
              expect_false( "m, d, f_" %in% cn )
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
              expect_false( "d, m_" %in% cn )
              expect_false( "d, m, f_" %in% cn )
              expect_true( "m, d_" %in% cn )
              expect_true( "m, d, f_" %in% cn )
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
              expect_false( "d, m_" %in% cn )
              expect_false( "d, m, f_" %in% cn )
              expect_true( "m, d_u" %in% cn )
              expect_true( "m, d, f_u" %in% cn )
          })
