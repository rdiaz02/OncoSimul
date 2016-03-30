## This takes ~ 22 seconds. Most of them in the very first test.
cat(paste("\n Starting oncoSimulSample-limits tests", date(), "\n"))

RNGkind("Mersenne-Twister")
data(examplePosets)

p701 <- examplePosets[["p701"]]


test_that("oncoSimulSample success with large num tries", {
    ## If nindiv small, sometimes you reach cancer at first try in every
    ## indiv.
    nindiv <- 70 ## this is decreased to increase speed
    p1 <- oncoSimulSample(nindiv, p701,
                          max.num.tries = 5000 * nindiv,
                          sampleEvery = 0.03, ## this to avoid large N all
                                             ## of a sudden
                          onlyCancer = TRUE)
    expect_true(p1$probCancer < 1)
    expect_true(p1$attemptsUsed > nindiv)
    expect_true(nrow(p1$popSummary) == nindiv)
})




test_that("oncoSimulSample success when no onlyCancer", {
    nindiv <- 4
    p4 <- oncoSimulSample(nindiv, p701,
                          sampleEvery = 0.03,
                          max.num.tries = nindiv,
                          onlyCancer = FALSE)
    expect_true(p4$probCancer ==  1)
    expect_true(p4$attemptsUsed ==  nindiv)
    expect_true(nrow(p4$popSummary) == nindiv)
})



test_that("oncoSimulSample exits with minimal num tries", {
    nindiv <- 50
    p5 <- oncoSimulSample(nindiv, p701,
                          initSize = 10,
                          finalTime = 50,
                          sampleEvery = 0.03,                          
                          max.num.tries = nindiv,
                          onlyCancer = TRUE)
    expect_true(p5$HittedMaxTries)
    expect_true(is.na(p5$popSummary))
})


test_that("oncoSimulSample exits with small num tries", {
    nindiv <- 50
    p6 <- oncoSimulSample(nindiv, p701,
                          initSize = 10,
                          finalTime = 50,
                          sampleEvery = 0.03,
                          max.num.tries = nindiv + 2,
                          onlyCancer = TRUE)
    expect_true(p6$HittedMaxTries)
    expect_true(is.na(p6$popSummary))
})


test_that("oncoSimulSample does not use drivers as stopping when there are none", {
    f1 <- allFitnessEffects(noIntGenes = rep(0.1, 5))
    expect_true(all(oncoSimulSample(5, f1,
                                    onlyCancer = FALSE,
                                    sampleEvery = 0.01,
                                    initSize = 1e4, ## no extinction in 2 units
                                    finalTime = 2)$popSummary[, "FinalTime"] == 2))
    oi <- allFitnessEffects(orderEffects =
                                          c("F > D" = -0.03, "D > F" = 0.4),
                                      noIntGenes = rexp(5, 10),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "F" = "f1, f2, f3",
                                            "D" = "d1, d2"),
                            drvNames = character(0))
    expect_true(all(oncoSimulSample(5, oi,
                                    onlyCancer = FALSE,
                                    sampleEvery = 0.01,
                                    initSize = 1e4, ## no extinction in 2 units
                                    finalTime = 2)$popSummary[, "FinalTime"] == 2))
})

cat(paste("\n Ending oncoSimulSample-limits tests", date(), "\n"))
