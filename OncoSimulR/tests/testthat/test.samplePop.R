test_that("exercising the sampling code, v1 objects", {
              data(examplePosets)
              p705 <- examplePosets[["p705"]]
              r1 <- oncoSimulIndiv(p705)
              p1 <- oncoSimulPop(4, p705, mc.cores = 2)
              expect_message(samplePop(p1),
                            "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(p1, typeSample = "single"),
                             "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(p1, typeSample = "whole", timeSample = "unif"),
                             "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(p1, typeSample = "single", timeSample = "unif"),
                            "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "whole", timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "single", timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "whole", timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "single", timeSample = "uniform"),
                            "Subjects by Genes matrix of 1 subjects and 37 genes")
          })

test_that("exercising the sampling code, v1 objects, possibly NAs", {
              data(examplePosets)
              p705 <- examplePosets[["p705"]]
              r1 <- oncoSimulIndiv(p705, onlyCancer = FALSE)
              p1 <- oncoSimulPop(4, p705, mc.cores = 2, onlyCancer = FALSE)
              expect_message(samplePop(p1),
                            "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(p1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(p1, typeSample = "whole",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(p1, typeSample = "single",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(p1, typeSample = "single",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(p1, typeSample = "single"),
                             "Subjects by Genes matrix of 4 subjects and 37 genes")
              expect_message(samplePop(r1),
                             "Subjects by Genes matrix of 1 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "whole",
                                       timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 37 genes")
              expect_message(samplePop(r1, typeSample = "single",
                                       timeSample = "uniform"),
                            "Subjects by Genes matrix of 1 subjects and 37 genes")
          })

test_that("exercising the sampling code, v1 objects, possibly NAs, more", {
              data(examplePosets)
              p1101 <- examplePosets[["p1101"]]
              r1 <- oncoSimulIndiv(p1101, onlyCancer = FALSE)
              p1 <- oncoSimulPop(4, p1101, mc.cores = 2, onlyCancer = FALSE)
              expect_message(samplePop(p1),
                            "Subjects by Genes matrix of 4 subjects and 41 genes")
              expect_message(samplePop(p1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 4 subjects and 41 genes")
              expect_message(samplePop(p1, typeSample = "whole",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 4 subjects and 41 genes")
              expect_message(samplePop(p1, typeSample = "single",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 4 subjects and 41 genes")
              expect_message(samplePop(p1, typeSample = "single",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 4 subjects and 41 genes")
              expect_message(samplePop(p1, typeSample = "single"),
                             "Subjects by Genes matrix of 4 subjects and 41 genes")
              expect_message(samplePop(r1),
                             "Subjects by Genes matrix of 1 subjects and 41 genes")
              expect_message(samplePop(r1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 41 genes")
              expect_message(samplePop(r1, typeSample = "whole",
                                       timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 41 genes")
              expect_message(samplePop(r1, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 41 genes")
              expect_message(samplePop(r1, typeSample = "single",
                                       timeSample = "uniform"),
                            "Subjects by Genes matrix of 1 subjects and 41 genes")
          })

test_that("exercising the sampling code, v2 objects, possibly NAs in output", {
              oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.3, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2") )
              o1 <- oncoSimulIndiv(oi, detectionSize = 1e4,
                                   onlyCancer = FALSE)
              o4 <- oncoSimulPop(2,
                                 oi, 
                                 detectionSize = 1e4,
                                 onlyCancer = FALSE)
              ## many of them are generating warnings, because sampling
              ## with pop size of 0. That is OK.
              expect_message(samplePop(o1),
                             "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o4),
                                           "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "uniform"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, timeSample = "unif",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, timeSample = "last",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
          })

test_that("exercising the sampling code, v2 objects, possibly NAs in output, more", {
              cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                                child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                                s = 0.1,
                                sh = -0.9,
                                typeDep = "MN")
              cbn1 <- allFitnessEffects(cs)
              o1 <- oncoSimulIndiv(cbn1, detectionSize = 1e4,
                                   onlyCancer = FALSE)
              o4 <- oncoSimulPop(4,
                                 cbn1, 
                                 detectionSize = 1e4,
                                 onlyCancer = FALSE,
                                 mc.cores = 2)
              ## many of them are generating warnings, because sampling
              ## with pop size of 0. That is OK.
              expect_message(samplePop(o1),
                             "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o4),
                                           "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "uniform"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, timeSample = "unif",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, timeSample = "last",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
          })


test_that("exercising the sampling code, v2 objects", {
              oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.3, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2") )
              o1 <- oncoSimulIndiv(oi, detectionSize = 1e4,
                                   onlyCancer = TRUE,
                                   max.num.tries = 5000)
              o4 <- oncoSimulPop(2,
                                 oi, 
                                 detectionSize = 1e4,
                                 onlyCancer = TRUE,
                                 max.num.tries = 5000)
              expect_message(samplePop(o1),
                             "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 10 genes")
              expect_message(samplePop(o4),
                                           "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "uniform"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, timeSample = "unif",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
              expect_message(samplePop(o4, timeSample = "last",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 2 subjects and 10 genes")
          })

test_that("exercising the sampling code, v2 objects, more", {
              cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                                child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                                s = 0.1,
                                sh = -0.9,
                                typeDep = "MN")
              cbn1 <- allFitnessEffects(cs)
              o1 <- oncoSimulIndiv(cbn1, detectionSize = 1e4,
                                   onlyCancer = TRUE,
                                   max.num.tries = 5000)
              o4 <- oncoSimulPop(4,
                                 cbn1, 
                                 detectionSize = 1e4,
                                 onlyCancer = TRUE,
                                 mc.cores = 2,
                                 max.num.tries = 5000)
              expect_message(samplePop(o1),
                             "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "single",
                                       timeSample = "unif"),
                            "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "unif"),
                             "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o1, typeSample = "whole",
                                       timeSample = "last"),
                             "Subjects by Genes matrix of 1 subjects and 6 genes")
              expect_message(samplePop(o4),
                                           "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "last"),
                            "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, typeSample = "single",
                                       timeSample = "uniform"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, timeSample = "unif",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
              expect_message(samplePop(o4, timeSample = "last",
                                       typeSample = "whole"),
                             "Subjects by Genes matrix of 4 subjects and 6 genes")
          })
