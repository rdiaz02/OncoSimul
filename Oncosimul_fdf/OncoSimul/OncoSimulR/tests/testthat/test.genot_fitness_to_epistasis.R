inittime <- Sys.time()
cat(paste("\n Starting genot_fitness_to_epistasis at", date(), "\n"))

## A minimal, simple check of this function which is a
## relict from the past
test_that("genot_fitness_to_epistasis works minimally", {

    m6 <- cbind(c(1, 1), c(1, 0), c(2, 3))
    fepi <- OncoSimulR:::genot_fitness_to_epistasis(m6)
    
    fem6 <- allFitnessEffects(epistasis = fepi)
    fem6b <- allFitnessEffects(genotFitness = m6)

    afe1 <- evalAllGenotypes(fem6, addwt = TRUE)
    afe2 <- evalAllGenotypes(fem6b, addwt = TRUE)

    ## recall breaking change!
    expect_true(!identical(afe1[, "Fitness"], afe2[, "Fitness"]))

    expect_true(identical(afe1[c(1, 2, 4), "Fitness"],
                          afe2[c(1, 2, 4), "Fitness"]))

    afe1mm <- afe1
    afe1mm$Fitness[3] <- 0
    
    expect_true(identical(afe1mm[, "Fitness"],
                          afe2[, "Fitness"]))
})

cat(paste("\n Ending genot_fitness_to_epistasis at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
