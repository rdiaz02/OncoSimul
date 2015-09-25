test_that("samplePop with oncoSimulIndiv object", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              oiI1 <- oncoSimulIndiv(oi, model = "Exp")
              expect_message(out <- samplePop(oiI1),
                             "Subjects by Genes matrix of 1 subjects and 10 genes.",
                             fixed = TRUE)
              expect_true(ncol(out) == 10)
              expect_true(nrow(out) == 1)
          })


test_that("McFl warning with small sampleEvery", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              expect_warning(oiI1 <- oncoSimulIndiv(oi, model = "McFL",
                                                    sampleEvery = 1,
                                                    onlyCancer = FALSE),
                             "With the McFarland model you often want smaller sampleEvery")
          })

test_that("mu < 0 error", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              expect_error(oiI1 <- oncoSimulIndiv(oi, mu = -1))                                                     
          })


test_that("keepEevery and sampleEvery consitency", {
              oi <- allFitnessEffects(orderEffects =
               c("F > D" = -0.3, "D > F" = 0.4),
               noIntGenes = rexp(5, 10),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2, f3",
                                "D" = "d1, d2") )
              expect_warning(oiI1 <- oncoSimulIndiv(oi, keepEvery = 2,
                                                    sampleEvery = 5),
                             "setting keepEvery <- sampleEvery",
                             fixed = TRUE)                                                     
          })


test_that("using old poset format: error if need param null", {
              data(examplePosets)
              p701 <- examplePosets[["p701"]]
              expect_error(oncoSimulIndiv(p701, s = NULL),
                           "You must specify all of")                                                     
          })

test_that("using old poset format: error if need param null", {
              data(examplePosets)
              p701 <- examplePosets[["p701"]]
              expect_error(oncoSimulIndiv(p701, sh = NULL),
                           "You must specify all of")                                                     
          })

test_that("using old poset format: error if need param null", {
              data(examplePosets)
              p701 <- examplePosets[["p701"]]
              expect_error(oncoSimulIndiv(p701, numPassengers = NULL),
                           "You must specify all of")                                                     
          })

test_that("verbosity options", {
              oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.3, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2") )
              expect_output(oncoSimulIndiv(oi, verbosity = 2,
                                           detectionSize = 1e4,
                                           onlyCancer = FALSE),
                            "Total Pop Size = ")
          })


test_that("printing oncosimul object", {
              oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.3, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2") )
              out <- oncoSimulIndiv(oi, 
                                    detectionSize = 1e4,
                                    onlyCancer = FALSE)
              expect_output(out,
                            "Individual OncoSimul trajectory with call")
          })


test_that("printing oncosimul pop object", {
              oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.3, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2") )
              out <- oncoSimulPop(4,
                                  oi, 
                                  detectionSize = 1e4,
                                  onlyCancer = FALSE)
              expect_output(out,
                            "Population of OncoSimul trajectories of 4 individuals")
          })



test_that("exercising oncoSimulSample, old format", {
              data(examplePosets)
              p701 <- examplePosets[["p701"]]
              expect_message(ofw <- oncoSimulSample(2, p701),
                             "Successfully sampled 2 individuals")
              expect_message(ofs <- oncoSimulSample(2, p701,
                                                    typeSample = "single"),
                             "Successfully sampled 2 individuals")
              expect_equal(dim(ofw$popSample), c(2, 7))
              expect_equal(dim(ofs$popSample), c(2, 7))
          })


test_that("exercising oncoSimulSample, new format", {
              pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                                        "TP53", "TP53", "MLL3"),
                                                    child = c("KRAS","SMAD4", "CDNK2A", 
                                                        "TP53", "MLL3",
                                                        rep("PXDN", 3), rep("TGFBR2", 2)),
                                                    s = 0.05,
                                                    sh = -0.3,
                                                    typeDep = "MN"))
              expect_message(pS <- oncoSimulSample(2, pancr),
                             "Successfully sampled 2 individuals")
              expect_message(
                  pSs <- oncoSimulSample(2,
                                                pancr,
                                                typeSample = "single"),
                  "Successfully sampled 2 individuals")
              expect_equal(dim(pS$popSample), c(2, 7))
              expect_equal(dim(pSs$popSample), c(2, 7))
          })


