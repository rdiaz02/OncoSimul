inittime <- Sys.time()
cat("\n Starting complete_fitness_landscape at", date(), "\n")

test_that("test.complete_fitness_landscape works as expected", {
  set.seed(NULL)
  for (i in 1:10) {

    replace_val <- runif(1, min = 0, max = 1e-3)
    if (i == 1) replace_val <- 0

    ## No removal here
    r1 <- rfitness(3)
    r1[, ncol(r1)] <- 1:nrow(r1)
    r1_trunc <- r1[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE]
    r1_recons <- complete_fitness_landscape(r1_trunc, fill = replace_val)
    expect_equal(r1[, -ncol(r1), drop = FALSE], r1_recons[, -ncol(r1), drop = FALSE])
    expect_equal(r1[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE],
                 r1_recons[r1_recons[, "Birth", drop = FALSE] >= 1, , drop = FALSE])
    expect_true(isTRUE(all(r1_recons[r1[, "Birth", drop = FALSE] < 1, "Birth", drop = FALSE] == replace_val)))



    ##### Removal
    r1 <- rfitness(3, min_accessible_genotypes = 1)
    r1_trunc <- r1[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE]
    r1_recons <- complete_fitness_landscape(r1_trunc, fill = replace_val)
    expect_equal(r1[, -ncol(r1), drop = FALSE], r1_recons[, -ncol(r1), drop = FALSE])
    expect_equal(r1[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE],
                 r1_recons[r1_recons[, "Birth", drop = FALSE] >= 1, , drop = FALSE])
    expect_true(isTRUE(all(r1_recons[r1[, "Birth", drop = FALSE] < 1, "Birth", drop = FALSE] == replace_val)))



    ## Shuffled
    r1_trunc <- r1_trunc[sample(1:nrow(r1_trunc)), , drop = FALSE]
    r1_recons <- complete_fitness_landscape(r1_trunc, fill = replace_val)
    expect_equal(r1[, -ncol(r1), drop = FALSE], r1_recons[, -ncol(r1), drop = FALSE])
    expect_equal(r1[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE],
                 r1_recons[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE])
    expect_true(isTRUE(all(r1_recons[r1[, "Birth", drop = FALSE] < 1, "Birth", drop = FALSE] == replace_val)))


    ## Single row input
    r1 <- rfitness(7, min_accessible_genotypes = 1)
    r1[, ncol(r1)] <- 0
    r1[1, ncol(r1)] <- 1

    colnames(r1)[ncol(r1)] <- "Birth"
    r1_trunc <- r1[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE]
    r1_recons <- complete_fitness_landscape(r1_trunc, fill = replace_val)
    expect_equal(r1[, -ncol(r1), drop = FALSE], r1_recons[, -ncol(r1), drop = FALSE])
    expect_equal(r1[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE],
                 r1_recons[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE])
    expect_true(isTRUE(all(r1_recons[r1[, "Birth", drop = FALSE] < 1, "Birth", drop = FALSE] == replace_val)))

    ## Single row input
    r1 <- rfitness(7, min_accessible_genotypes = 1)
    r1[, ncol(r1)] <- 0
    r1[3, ncol(r1)] <- 1
    colnames(r1)[ncol(r1)] <- "Birth"
    r1_trunc <- r1[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE]
    r1_recons <- complete_fitness_landscape(r1_trunc, fill = replace_val)
    expect_equal(r1[, -ncol(r1), drop = FALSE], r1_recons[, -ncol(r1), drop = FALSE])
    expect_equal(r1[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE],
                 r1_recons[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE])
    expect_true(isTRUE(all(r1_recons[r1[, "Birth", drop = FALSE] < 1, "Birth", drop = FALSE] == replace_val)))

    ## Now, none accessible
    r1 <- rfitness(7, min_accessible_genotypes = 1)
    r1[, ncol(r1)] <- 0
    r1_trunc <- r1[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE]
    r1_recons <- complete_fitness_landscape(r1_trunc, fill = replace_val)
    expect_equal(r1[, -ncol(r1), drop = FALSE], r1_recons[, -ncol(r1), drop = FALSE])
    expect_equal(r1[r1[, "Birth", drop = FALSE] >= 1, , drop = FALSE],
                 r1_recons[r1_recons[, "Birth", drop = FALSE] >= 1, , drop = FALSE])
    expect_true(isTRUE(all(r1_recons[r1[, "Birth", drop = FALSE] < 1, "Birth", drop = FALSE] == replace_val)))





    ####################

    #### Now, Fitness instead of Birth for last column
    ## No removal
    r2 <- rfitness(7, min_accessible_genotypes = 1)
    r2[, ncol(r2)] <- 1:nrow(r2)
    colnames(r2)[ncol(r2)] <- "Fitness"
    r2_trunc <- r2[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE]
    r2_recons <- complete_fitness_landscape(r2_trunc, fill = replace_val)
    expect_equal(r2[, -ncol(r2), drop = FALSE], r2_recons[, -ncol(r2), drop = FALSE])
    expect_equal(r2[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE],
                 r2_recons[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE])
    expect_true(isTRUE(all(r2_recons[r2[, "Fitness", drop = FALSE] < 1, "Fitness", drop = FALSE] == replace_val)))

    ## Removal
    r2 <- rfitness(7, min_accessible_genotypes = 1)
    colnames(r2)[ncol(r2)] <- "Fitness"
    r2_trunc <- r2[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE]
    r2_recons <- complete_fitness_landscape(r2_trunc, fill = replace_val)
    expect_equal(r2[, -ncol(r2), drop = FALSE], r2_recons[, -ncol(r2), drop = FALSE])
    expect_equal(r2[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE],
                 r2_recons[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE])
    expect_true(isTRUE(all(r2_recons[r2[, "Fitness", drop = FALSE] < 1, "Fitness", drop = FALSE] == replace_val)))




    ## Shuffled
    r2_trunc <- r2_trunc[sample(1:nrow(r2_trunc)), , drop = FALSE]
    r2_recons <- complete_fitness_landscape(r2_trunc, fill = replace_val)
    expect_equal(r2[, -ncol(r2), drop = FALSE], r2_recons[, -ncol(r2), drop = FALSE])
    expect_equal(r2[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE],
                 r2_recons[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE])
    expect_true(isTRUE(all(r2_recons[r2[, "Fitness", drop = FALSE] < 1, "Fitness", drop = FALSE] == replace_val)))


    ## Single row input
    r2 <- rfitness(7, min_accessible_genotypes = 1)
    r2[, ncol(r2)] <- 0
    r2[1, ncol(r2)] <- 1

    colnames(r2)[ncol(r2)] <- "Fitness"
    r2_trunc <- r2[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE]
    r2_recons <- complete_fitness_landscape(r2_trunc, fill = replace_val)
    expect_equal(r2[, -ncol(r2), drop = FALSE], r2_recons[, -ncol(r2), drop = FALSE])
    expect_equal(r2[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE],
                 r2_recons[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE])
    expect_true(isTRUE(all(r2_recons[r2[, "Fitness", drop = FALSE] < 1, "Fitness", drop = FALSE] == replace_val)))

    ## Single row input
    r2 <- rfitness(7, min_accessible_genotypes = 1)
    r2[, ncol(r2)] <- 0
    r2[3, ncol(r2)] <- 1
    colnames(r2)[ncol(r2)] <- "Fitness"
    r2_trunc <- r2[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE]
    r2_recons <- complete_fitness_landscape(r2_trunc, fill = replace_val)
    expect_equal(r2[, -ncol(r2), drop = FALSE], r2_recons[, -ncol(r2), drop = FALSE])
    expect_equal(r2[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE],
                 r2_recons[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE])
    expect_true(isTRUE(all(r2_recons[r2[, "Fitness", drop = FALSE] < 1, "Fitness", drop = FALSE] == replace_val)))



    ## Now, none accessible
    r2 <- rfitness(7)
    r2[, ncol(r2)] <- 0
    colnames(r2)[ncol(r2)] <- "Fitness"
    r2_trunc <- r2[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE]
    r2_recons <- complete_fitness_landscape(r2_trunc, fill = replace_val)
    expect_equal(r2[, -ncol(r2), drop = FALSE], r2_recons[, -ncol(r2), drop = FALSE])
    expect_equal(r2[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE],
                 r2_recons[r2[, "Fitness", drop = FALSE] >= 1, , drop = FALSE])
    expect_true(isTRUE(all(r2_recons[r2[, "Fitness", drop = FALSE] < 1, "Fitness", drop = FALSE] == replace_val)))
  }

})


cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
