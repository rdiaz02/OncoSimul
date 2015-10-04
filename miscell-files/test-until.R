library(OncoSimulR)
library(testthat)

data(examplePosets)
p701 <- examplePosets[["p701"]]

i <-  as.integer(runif(1, min = 0, max = 1e9))

count <- 0
while(TRUE) {
    count <- count + 1
    cat("\n\n Iteration is ", count, "\n")
    ## i <- i + 1
    i <-  as.integer(runif(1, min = 0, max = 1e9))
    set.seed(i)
    cat("\n\n seed is ", i, "\n")
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
    p1$probCancer
    
    })

    
}
