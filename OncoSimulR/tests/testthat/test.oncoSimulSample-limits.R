data(examplePosets)

p701 <- examplePosets[["p701"]]


test_that("oncoSimulSample success with large num tries", {
    ## If nindiv small, sometimes you reach cancer at first try in every
    ## indiv.
    nindiv <- 70 ## this is decreased to increase speed
    p1 <- oncoSimulSample(nindiv, p701,
                          max.num.tries = 5000 * nindiv,
                          sampleEvery = 0.3, ## this to avoid large N all
                                             ## of a sudden
                          onlyCancer = TRUE)
    expect_true(p1$probCancer < 1)
    expect_true(p1$attemptsUsed > nindiv)
    expect_true(nrow(p1$popSummary) == nindiv)
})




test_that("oncoSimulSample success when no onlyCancer", {
    nindiv <- 4
    p4 <- oncoSimulSample(nindiv, p701,
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
                          max.num.tries = nindiv + 2,
                          onlyCancer = TRUE)
    expect_true(p6$HittedMaxTries)
    expect_true(is.na(p6$popSummary))
})

