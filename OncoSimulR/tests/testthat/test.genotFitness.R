test_that("Conversion for data frame", {
    (m4 <- data.frame(G = c("A, B", "A", "WT", "B"), F = c(3, 2, 1, 4)))
    fem4 <- allFitnessEffects(genotFitness = m4)
    expect_true(all(evalAllGenotypes(fem4, addwt = TRUE, order = FALSE) ==
                    m4[c(3, 2, 4, 1), ]))
})


test_that("Conversion for matrix", {

    m5 <- cbind(c(1, 0, 1, 0), c(0, 1, 1, 0), c(2, 3, 5.5, 1))
    fem5 <- allFitnessEffects(genotFitness = m5)
    evalAllGenotypes(fem5, addwt = TRUE, order = FALSE)
    expect_true(all(evalAllGenotypes(fem5, addwt = TRUE, order = FALSE)[, 2] ==
                    m5[c(4, 1, 2, 3), 3]))
    
})


test_that("Conversion for incomplete matrix", {

    m6 <- cbind(c(1, 1), c(1, 0), c(2, 3))
    fem6 <- allFitnessEffects(genotFitness = m6)
    evalAllGenotypes(fem6, addwt = TRUE, order = FALSE)

    expect_true(all(evalAllGenotypes(fem6, addwt = TRUE, order = FALSE)[, 2] ==
                    c(1, m6[2, 3], 1, m6[1, 3])))
    
})


test_that("We fail with order", {
    expect_error(OncoSimulR:::from_genotype_fitness(data.frame(G = "A > B", F = 1)))
})

test_that("We fail with wrong separator, :", {
    expect_error(OncoSimulR:::from_genotype_fitness(data.frame(G = "A : B", F = 1)))
})
