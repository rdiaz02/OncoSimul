data(examplePosets)

p701 <- examplePosets[["p701"]]

nindiv <- 4


test_that("oncoSimulSample success with large num tries", {
    p1 <- oncoSimulSample(nindiv, p701, max.num.tries = 5000 * nindiv,
                          onlyCancer = TRUE)
    expect_true(p1$probCancer < 1)
    expect_true(p1$attemptsUsed > nindiv)
    expect_true(nrow(p1$popSummary) == nindiv)
})



test_that("oncoSimulSample success when no onlyCancer", {
    p4 <- oncoSimulSample(nindiv, p701,
                          max.num.tries = nindiv,
                          onlyCancer = FALSE)
    expect_true(p4$probCancer ==  1)
    expect_true(p4$attemptsUsed ==  nindiv)
    expect_true(nrow(p4$popSummary) == nindiv)
})



test_that("oncoSimulSample exits with minimal num tries", {
    p5 <- oncoSimulSample(nindiv, p701,
                          initSize = 10,
                          finalTime = 50,
                          max.num.tries = nindiv,
                          onlyCancer = TRUE)
    expect_true(p5$HittedMaxTries)
    expect_true(is.na(p5$popSummary))
})


test_that("oncoSimulSample exits with small num tries", {
    p6 <- oncoSimulSample(nindiv, p701,
                          initSize = 10,
                          finalTime = 50,
                          max.num.tries = nindiv + 1,
                          onlyCancer = TRUE)
    expect_true(p6$HittedMaxTries)
    expect_true(is.na(p6$popSummary))
})

