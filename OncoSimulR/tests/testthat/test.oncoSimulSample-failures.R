test_that("oncoSimulSample out of time triggered", {

              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              
              expect_message(out <- oncoSimulSample(100, oi,
                                                    max.wall.time.total = 1),
                          "Run out of time")

              expect_true(out$HittedWallTime)
          })

test_that("oncoSimulSample out of time triggered, 2", {
              pancr <- allFitnessEffects(
                  data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4",
                                 "CDNK2A", 
                                 "TP53", "TP53", "MLL3"),
                             child = c("KRAS","SMAD4", "CDNK2A", 
                                 "TP53", "MLL3",
                                 rep("PXDN", 3), rep("TGFBR2", 2)),
                             s = 0.05,
                             sh = -0.3,
                             typeDep = "MN"))
              expect_message(out <- oncoSimulSample(100, pancr,
                                                    max.wall.time.total = 1),
                          "Run out of time")
              expect_true(out$HittedWallTime)
          })



test_that("oncoSimulSample out of attempts triggered", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              expect_message(out <- oncoSimulSample(10, oi,
                                                 max.num.tries.total = 11),
                          "Run out of attempts (in C++)", fixed = TRUE)
              expect_true(out$HittedMaxTries)
          })

test_that("oncoSimulSample out of attempts triggered, 2", {
              pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                          "TP53", "TP53", "MLL3"),
                                      child = c("KRAS","SMAD4", "CDNK2A", 
                                          "TP53", "MLL3",
                                          rep("PXDN", 3), rep("TGFBR2", 2)),
                                      s = 0.05,
                                      sh = -0.3,
                                      typeDep = "MN"))
              expect_message(out <- oncoSimulSample(10, pancr,
                                                 max.num.tries.total = 12),
                          "Run out of attempts (in C++)", fixed = TRUE)
              expect_true(out$HittedMaxTries)
          })




test_that("oncoSimulSample out of time in C++ triggered", {
              pancr <- allFitnessEffects(data.frame(parent = c("Root",
                                                        rep("KRAS", 4),
                                                        "SMAD4", "CDNK2A", 
                                                        "TP53", "TP53", "MLL3"),
                                                    child = c("KRAS","SMAD4",
                                                        "CDNK2A", 
                                                        "TP53", "MLL3",
                                                        rep("PXDN", 3),
                                                        rep("TGFBR2", 2)),
                                                    s = 0.05,
                                                    sh = -0.3,
                                                    typeDep = "MN"))
              
              expect_message(out <- oncoSimulSample(1, pancr, "McFL",
                                                 detectionSize = 1e9,
                                                 detectionDrivers = 4,
                                                 max.wall.time.total = 0.5),
                             "Run out of time (in C++)", fixed = TRUE)
              
              expect_true(out$HittedWallTime)
          })


test_that("oncoSimulSample C++ exception triggered", {
              pancr <- allFitnessEffects(data.frame(parent = c("Root",
                                                        rep("KRAS", 4),
                                                        "SMAD4", "CDNK2A", 
                                                        "TP53", "TP53", "MLL3"),
                                                    child = c("KRAS","SMAD4",
                                                        "CDNK2A", 
                                                        "TP53", "MLL3",
                                                        rep("PXDN", 3),
                                                        rep("TGFBR2", 2)),
                                                    s = 0.05,
                                                    sh = -0.3,
                                                    typeDep = "MN"))
              expect_message(out <- oncoSimulSample(1, pancr, "McFL", mu = 0),
                             "Unrecoverable exception (in C++)", fixed = TRUE)
              expect_true(out$UnrecoverExcept)
          })

