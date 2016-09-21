test_that("Expect output", {

    expect_output(print(rfitness(4)), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, reference = "max")), "Fitness",
                  fixed = TRUE)
    expect_output(print(rfitness(3, reference = c(1, 0, 1))), "Fitness",
                  fixed = TRUE)
    expect_output(print(rfitness(5, scale = c(1, 7))), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(5, scale = c(1, 4), wt_is_1 = FALSE,
                           log = TRUE)), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, reference = "random2")), "Fitness",
                  fixed = TRUE)
    expect_output(print(rfitness(4, reference = "random2",
                                 min_accessible_genotypes = 6)), "Fitness",
                  fixed = TRUE)
    expect_output(print(rfitness(3, reference = "random2",
                                 min_accessible_genotypes = 6)), "Fitness",
                  fixed = TRUE)
    expect_output(print(rfitness(3, reference = "max",
                                 min_accessible_genotypes = 6)), "Fitness",
                  fixed = TRUE)
})
