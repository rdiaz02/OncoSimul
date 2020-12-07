inittime <- Sys.time()
cat(paste("\n Starting FDF-small-fitness-specs", date(), "\n"))


## Testing fitness specs with missing genotypes and with letters too

test_that("We can run and equal with letters", {


    g1 <- data.frame(
        Genotype = c("WT", "A", "A, B"), 
        Fitness = c("1",
                    "1.3 + 1.2 * n_1_2",
                    "1.1"),
        stringsAsFactors = FALSE)
    fg1 <- allFitnessEffects(genotFitness = g1, 
                             frequencyDependentFitness = TRUE)

    g11 <- data.frame(
        Genotype = c("WT", "A", "A, B"), 
        Fitness = c("1",
                    "1.3 + 1.2 * n_A_B",
                    "1.1"),
        stringsAsFactors = FALSE)
    fg11 <- allFitnessEffects(genotFitness = g11, 
                              frequencyDependentFitness = TRUE)

    ## FIXME
    ## So, what does spPopSizes refer to ?? that is ambiguous here. 
    ## It seems it is the Genotypes in the original Genotype spec.
    ## but they are the ones left in fg1$Genotype
    ## Since we cannot now what was in g1, do not allow for this.
    efg1 <- evalAllGenotypes(fg1, spPopSizes = c(9, 2, 6))
    efg11 <- evalAllGenotypes(fg11, spPopSizes = c(9, 2, 6))

    expect_identical(efg1, efg11)
    ## is that correct?
    expect_identical(efg1$Fitness, c(1, 1.3 + 1.2 * 6, 0,  1.1))

    
    g1b <- data.frame(
        Genotype = c("WT", "A", "A, B", "B"), 
        Fitness = c("1",
                    "1.3 + 1.2 * n_1_2",
                    "1.1",
                    "0"),
        stringsAsFactors = FALSE)
    fg1b <- allFitnessEffects(genotFitness = g1b, 
                         frequencyDependentFitness = TRUE)
    ## it thinks n_1_2 is 6. This is right
    (ea1b <- evalAllGenotypes(fg1b, spPopSizes = c(9, 2, 0, 6)))
    
    ## is it correct? Yes
    stopifnot(identical(ea1b[, "Fitness"], c(1, 1.3 + 1.2 * 6, 0, 1.1)))
    
 
    set.seed(1)
    rfg1 <- oncoSimulIndiv(fg1, initMutant = "A", onlyCancer = FALSE, finalTime = 500, seed = NULL)
    set.seed(1)
    rfg1b <- oncoSimulIndiv(fg1b, initMutant = "A", onlyCancer = FALSE, finalTime = 500, seed = NULL)
    
    expect_equivalent(rfg1, rfg1b)
    
    g1c <- data.frame(
        Genotype = c("WT", "A", "A, B", "B"), 
        Fitness = c("1",
                "1.3 + 1.2 * n_A_B",
                "1.1",
                "0"),
        stringsAsFactors = FALSE)
    fg1c <- allFitnessEffects(genotFitness = g1c,
                              frequencyDependentFitness = TRUE)
    set.seed(1)
    rfg1c <- oncoSimulIndiv(fg1c, initMutant = "A", onlyCancer = FALSE, finalTime = 500, seed = NULL)
    expect_equivalent(rfg1, rfg1c)
    
    set.seed(1)
    oncoSimulIndiv(fg1, initMutant = "A, B", onlyCancer = FALSE, finalTime = 500, seed = NULL)
} )


test_that("We can run and equal in different order" , {
    ## Actually, spPopSizes must be order as in the table of fitnessLandscape_df:
    ## spPopSizes are in the same order as given in "Genotype" in the output of
    ## evalAllGenotypes

    g2 <- data.frame(
        Genotype = c("WT", "A", "A, B", "B"), 
        Fitness = c("1",
                    "1 + 2 * n_1_2",
                    "1 + 2 * n_2",
                    "1 + 2 * n_1"),
        stringsAsFactors = FALSE)
    fg2 <- allFitnessEffects(genotFitness = g2, 
                             frequencyDependentFitness = TRUE)

    ## identical, except order of genotypes changes. But the specification is the same
    g2b <- data.frame(
        Genotype = c("WT", "A", "B", "A, B"), 
        Fitness = c("1",
                    "1 + 2 * n_1_2",
                    "1 + 2 * n_1",
                    "1 + 2 * n_2"),
        stringsAsFactors = FALSE)
    fg2b <- allFitnessEffects(genotFitness = g2b, 
                              frequencyDependentFitness = TRUE)

    ## spPopSizes are for genotypes AT, A, B, AB
    ofg2 <- evalAllGenotypes(fg2, spPopSizes = c(9, 2, 6, 3))
    ofg2b <- evalAllGenotypes(fg2b, spPopSizes = c(9, 2, 6, 3))
    ## Are they correct?
    expect_identical(ofg2, ofg2b)
    out_expec_ofg2 <- c(1, 1 + 2 * 3 , 1 + 2 * 2, 1 + 2 * 6)
    expect_identical(ofg2[, "Fitness"], out_expec_ofg2)
})


cat(paste("\n Ending FDF-small-fitness-specs", date(), "\n"))


cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

