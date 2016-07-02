## This test takes about 30 seconds

data(examplePosets)

p701 <- examplePosets[["p701"]]

test_that("oncoSimulSample success with large num tries", {
    ## If nindiv small, sometimes you reach cancer at first try in every
    ## indiv.
    nindiv <- 200 
    p1 <- oncoSimulSample(nindiv, p701,
                          max.num.tries = 5000 * nindiv,
                          sampleEvery = 0.03, ## this to avoid large N all
                                             ## of a sudden
                          onlyCancer = TRUE, detectionProb = NA)
    expect_true(p1$probCancer < 1)
    expect_true(p1$attemptsUsed > nindiv)
    expect_true(nrow(p1$popSummary) == nindiv)
})

