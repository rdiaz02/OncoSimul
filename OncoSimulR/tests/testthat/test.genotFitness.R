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


test_that("Dividing fitness by wt and missing genot is 1", {
    
    m7 <- cbind(c(1, 1, 0), c(1, 0, 0), c(2, 3, 5))
    expect_message(fem7 <- allFitnessEffects(genotFitness = m7),
                   "Fitness of wildtype != 1. Dividing all fitnesses by fitness of wildtype.",
                   fixed = TRUE)
    expect_equivalent(evalAllGenotypes(fem7, order = FALSE, addwt = TRUE)[, 2],
                      c(1, 3/5, 1, 2/5))
})

test_that("The WT is added if absent", {
    m7 <- cbind(c(1, 1, 0), c(1, 0, 0), c(2, 3, 5))
    fem7 <- allFitnessEffects(genotFitness = m7)
    ag <- evalAllGenotypes(fem7, order = FALSE)
    ## internal call
    expect_equivalent(OncoSimulR:::allGenotypes_to_matrix(ag)[, 3],
                      c(1, 3/5, 1, 2/5))
    ## the user visible, which is via plotFitnessLandscape -> to_Fitness_Matrix
    plot(ag)
})


test_that("genotFitness not combined with others", {

    m7 <- cbind(c(1, 1, 0), c(1, 0, 0), c(2, 3, 5))
    
    expect_error(allFitnessEffects(genotFitness = m7,
                                   noIntGenes = c(.1, .2)),
                 "You have a non-null genotFitness.",
                 fixed = TRUE)

    expect_error(allFitnessEffects(genotFitness = m7,
                                   epistasis = c("A" = .3,
                                            "B" = .1,
                                            "C" = .2)),
                 "You have a non-null genotFitness.",
                 fixed = TRUE)
    
    expect_error(allFitnessEffects(genotFitness = m7,
                                   orderEffects = c("G > H" = 1.1)),
                 "You have a non-null genotFitness.",
                 fixed = TRUE)


    p3 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c", "f"),
                  child = c("a", "b", "d", "e", "c", "c", "f", "f", "g", "g"),
                  s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                  sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                  typeDep = c(rep("--", 4), 
                              "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    
    expect_error(allFitnessEffects(rT = p3, genotFitness = m7),
                 "You have a non-null genotFitness.",
                 fixed = TRUE)
   
})


