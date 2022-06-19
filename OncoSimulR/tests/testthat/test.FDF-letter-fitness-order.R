inittime <- Sys.time()
cat(paste("\n Starting FDF-letter-fitness-order", date(), "\n"))


## Testing fitness specs with missing genotypes and with letters too
## Also tests that are failed by version from October-December 2020 at least.
test_that("We can run and equal with letters", {


    g1 <- data.frame(
        Genotype = c("WT", "A", "A, B"), 
        Fitness = c("1",
                    "1.3 + 1.2 * n_1_2",
                    "1.1"),
        stringsAsFactors = FALSE)
    suppressWarnings(fg1 <- allFitnessEffects(genotFitness = g1, 
                             frequencyDependentFitness = TRUE))

    g11 <- data.frame(
        Genotype = c("WT", "A", "A, B"), 
        Fitness = c("1",
                    "1.3 + 1.2 * n_A_B",
                    "1.1"),
        stringsAsFactors = FALSE)
    suppressWarnings(fg11 <- allFitnessEffects(genotFitness = g11, 
                              frequencyDependentFitness = TRUE))

    ## So, what does spPopSizes refer to ?? that is ambiguous here. 
    ## It seems it is the Genotypes in the original Genotype spec.
    ## but they are the ones left in fg1$Genotype
    ## Since we cannot now what was in g1, emits warnings, that we silence here
    efg1 <- suppressWarnings(evalAllGenotypes(fg1, spPopSizes = c(9, 2, 6)))
    efg11 <- suppressWarnings(evalAllGenotypes(fg11, spPopSizes = c(9, 2, 6)))

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
    suppressWarnings(fg1b <- allFitnessEffects(genotFitness = g1b, 
                         frequencyDependentFitness = TRUE))
    ## it thinks n_1_2 is 6. This is right

    (ea1b <- suppressWarnings(evalAllGenotypes(fg1b,
                                               spPopSizes = c(9, 2, 0, 6))))
    
    
    ## is it correct? Yes
    stopifnot(identical(ea1b[, "Fitness"], c(1, 1.3 + 1.2 * 6, 0, 1.1)))
    
 
    set.seed(1)
    rfg1 <- oncoSimulIndiv(fg1, initMutant = "A", onlyCancer = FALSE,
                           finalTime = 500, seed = NULL)
    set.seed(1)
    rfg1b <- oncoSimulIndiv(fg1b, initMutant = "A", onlyCancer = FALSE,
                            finalTime = 500, seed = NULL)
    
    expect_equivalent(rfg1, rfg1b)
    
    g1c <- data.frame(
        Genotype = c("WT", "A", "A, B", "B"), 
        Fitness = c("1",
                "1.3 + 1.2 * n_A_B",
                "1.1",
                "0"),
        stringsAsFactors = FALSE)
    suppressWarnings(fg1c <- allFitnessEffects(genotFitness = g1c,
                              frequencyDependentFitness = TRUE))
    set.seed(1)
    rfg1c <- oncoSimulIndiv(fg1c, initMutant = "A", onlyCancer = FALSE,
                            finalTime = 500, seed = NULL)
    expect_equivalent(rfg1, rfg1c)
    
    set.seed(1)
    oncoSimulIndiv(fg1, initMutant = "A, B", onlyCancer = FALSE,
                   finalTime = 500, seed = NULL)
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
    suppressWarnings(fg2 <- allFitnessEffects(genotFitness = g2, 
                             frequencyDependentFitness = TRUE))

    ## identical, except order of genotypes changes. But the specification is the same
    g2b <- data.frame(
        Genotype = c("WT", "A", "B", "A, B"), 
        Fitness = c("1",
                    "1 + 2 * n_1_2",
                    "1 + 2 * n_1",
                    "1 + 2 * n_2"),
        stringsAsFactors = FALSE)
    suppressWarnings(fg2b <- allFitnessEffects(genotFitness = g2b, 
                              frequencyDependentFitness = TRUE))

    ## spPopSizes are for genotypes AT, A, B, AB
    ofg2 <- suppressWarnings(evalAllGenotypes(fg2, spPopSizes = c(9, 2, 6, 3)))
    ofg2b <- suppressWarnings(evalAllGenotypes(fg2b, spPopSizes = c(9, 2, 6, 3)))
    ## Are they correct?
    expect_identical(ofg2, ofg2b)
    out_expec_ofg2 <- c(1, 1 + 2 * 3 , 1 + 2 * 2, 1 + 2 * 6)
    expect_identical(ofg2[, "Fitness"], out_expec_ofg2)
})



test_that("Breaks as it should", {
    df1 <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                      Fitness = c("1",
                                  "n_1_3", ## BA
                                  "n_2_3", ## CA
                                  "n_1_2",
                                  "n_1",
                                  "n_2"
                                  ))
    suppressWarnings(suppressMessages(adf1 <- allFitnessEffects(genotFitness = df1,
                              frequencyDependentFitness = TRUE)))
    ## (adf1)
    expect_error(suppressWarnings(evalAllGenotypes(adf1, spPopSizes = 1:6)))
    ## Breaks in old and new: n_2_3
    df1 <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                  Fitness = c("1",
                              "f_1_3", ## BA
                              "f_2_3", ## CA
                              "f_1_2",
                              "f_1",
                              "f_2"
                              ))
    suppressWarnings(suppressMessages(adf1 <- allFitnessEffects(genotFitness = df1,
                              frequencyDependentFitness = TRUE)))
    
    ## (adf1)
    expect_error(suppressWarnings(evalAllGenotypes(adf1, spPopSizes = 1:6))) ## Breaks in old
    ## and new: f_2_3
})


## Good example of how easy it was to screw up
test_that("eval fitness gives wrong answer, as misspecified", {
    df2 <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                  Fitness = c("1",
                              "n_1_3", ## you want BA
                              "n_3", 
                              "n_3",
                              "n_1",
                              "n_2"
                              ))
    suppressWarnings(adf2 <- allFitnessEffects(genotFitness = df2,
                              frequencyDependentFitness = TRUE))

    (adf2)
    suppressWarnings(gg <- evalAllGenotypes(adf2, spPopSizes = 1:6))
    expect_equal(gg[gg$Genotype == "B", "Fitness"], 6) 
    ## Wrong: gives for B fitness of using CA, not BA (both old and new versions)
})


test_that("Can deal with some user misorderig (and naming spPopSizes)" , {
    ## Older one could not
    df3 <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                  Fitness = c("1",
                              "n_B_A", ## you want BA
                              "n_A", 
                              "n_A",
                              "n_B",
                              "n_C"
                              ))
    suppressWarnings(adf3 <- allFitnessEffects(genotFitness = df3,
                              frequencyDependentFitness = TRUE))
    (adf3)
    spPopSizes <- 1:6
    names(spPopSizes) <- c("WT", "B", "C", "A", "B, A", "C, A")
    gg <- evalAllGenotypes(adf3, spPopSizes = spPopSizes) 
    expect_equal(gg[gg$Genotype == "B", "Fitness"], 5) 
})


test_that("Correct values of evalAllGenotypes and substitutions", {
    ## Old one gives wrong result
    df3 <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                      Fitness = c("1",
                                  "f_A_B", ## you want BA
                                  "f_A", 
                                  "f_A",
                                  "f_B",
                                  "f_C"
                                  ))
    suppressWarnings(adf3 <- allFitnessEffects(genotFitness = df3,
                              frequencyDependentFitness = TRUE))
    spPopSizes <- 1:6
    names(spPopSizes) <- c("WT", "B", "C", "A", "B, A", "C, A")
    ## spPopSizes <- c(WT = 1, B = 2, C = 3, A = 4, AB=5, CA = 6)
    
    gg <- evalAllGenotypes(adf3, spPopSizes = spPopSizes) 
    
    ## should give
    ## c(1, c(A = 4, B = 5, C = 4, AB = 2, AC = 3,  BC = 0, ABC = 0)/sum(spPopSizes))
    expect_equal(gg[, "Fitness"],
                 c(1, c(4, 5, 4, 2, 3, 0, 0)/sum(spPopSizes))
                 )

    r <- data.frame(rfitness(3))
    r <- r[, c(2, 3, 1, 4)]
    r <- r[c(1, 3, 4, 2, 5, 6), ]
	
	colnames(r)[which(colnames(r) == "Birth")] <- "Fitness"
    r$Fitness <- c("1",
                   "f_A_B", ## you want BA
                   "f_A", 
                   "f_A_C",
                   "f_B",
                   "f_C"
                   )
    df3m <- r

    df3m2 <- data.frame(rfitness(3))[, c(2, 3, 1, 4)][1:6, ]
	
	colnames(df3m2)[which(colnames(df3m2) == "Birth")] <- "Fitness"
    df3m2$Fitness <- c("1",
                       "f_A_C",
                       "f_A_B", ## you want BA
                       "f_A", 
                       "f_B",
                       "f_C"
                       )
    ## Note the columns are sorted as per message
    suppressWarnings(adf3m <- allFitnessEffects(genotFitness = df3m,
                               frequencyDependentFitness = TRUE))
    suppressWarnings(adf3m2 <- allFitnessEffects(genotFitness = df3m2,
                                frequencyDependentFitness = TRUE))

    o1 <- evalAllGenotypes(adf3m, spPopSizes = spPopSizes)

    o2 <- evalAllGenotypes(adf3m2, spPopSizes = spPopSizes)
    expect_equal(o1, o2)
    
})


test_that("Fail if unequal names spPopSizes", {
    r <- data.frame(rfitness(3))
    r <- r[, c(2, 3, 1, 4)]
    r <- r[c(1, 3, 4, 2, 5, 6), ]
	
	colnames(r)[which(colnames(r) == "Birth")] <- "Fitness"
    r$Fitness <- c("1",
                   "f_A_B", ## you want BA
                   "f_A", 
                   "f_A_C",
                   "f_B",
                   "f_C"
                   )
    df3m <- r

    df3m2 <- data.frame(rfitness(3))[, c(2, 3, 1, 4)][1:6, ]
	colnames(df3m2)[which(colnames(df3m2) == "Birth")] <- "Fitness"
	
    df3m2$Fitness <- c("1",
                       "f_A_C",
                       "f_A_B", ## you want BA
                       "f_A", 
                       "f_B",
                       "f_C"
                       )
    ## Note the columns are sorted as per message
    suppressWarnings(adf3m <- allFitnessEffects(genotFitness = df3m,
                               frequencyDependentFitness = TRUE))
    spPopSizes <- 1:6
    names(spPopSizes) <- c("WT", "B", "C", "A", "B, A", "B, C")
    expect_error(evalAllGenotypes(adf3m, spPopSizes = spPopSizes),
                 "Genotype names in spPopSizes differ ", fixed = TRUE)
})




test_that("Correct values of evalAllGenotypes and substitutions, some math and changed order genes in n expression", {
    ## Old one gives wrong result
    df3 <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                      Fitness = c("1",
                                  "log(n_A_B)", ## you want BA
                                  "sqrt(n_A) + 3 * exp(n_C)", 
                                  "-1 * n_A + 7 * (n_C_A > 2)",
                                  "2 - 1/n_B",
                                  "min(n_C, n_B) - 0.1 * (n_B_A > 2)"))
                                  
    suppressWarnings(adf3 <- allFitnessEffects(genotFitness = df3,
                              frequencyDependentFitness = TRUE))
    spPopSizes <- 1:6
    names(spPopSizes) <- c("WT", "B", "C", "A", "B, A", "C, A")
    
    
    gg <- evalAllGenotypes(adf3, spPopSizes = spPopSizes) 
    
    ## should give
   
    expect_equal(gg[, "Fitness"],
                 c(1, c(-4 + 7,
                        log(5),
                        sqrt(4) + 3 * exp(3),
                        2 - 1/2,
                        2 - 0.1,
                        0,
                        0))
                 )


    df3x <- data.frame(Genotype = c("WT", "B", "C", "A", "B, A", "C, A"),
                      Fitness = c("n_A + 0.1",
                                  "log(n_A_B)", ## you want BA
                                  "sqrt(n_A) + 3 * exp(n_C)", 
                                  "-1 * n_A + 7 * (n_C_A > 2)",
                                  "2 - 0.1/n_B",
                                  "min(n_C, n_B) - 1 * (n_B_A > 2)"))
    suppressWarnings(adf3x <- allFitnessEffects(genotFitness = df3x,
                              frequencyDependentFitness = TRUE))
    spPopSizes <- 1:6
    names(spPopSizes) <- c("WT", "B", "C", "A", "B, A", "C, A")
    ggx <- evalAllGenotypes(adf3x, spPopSizes = spPopSizes) 
    
    expect_equal(ggx[, "Fitness"],
                 c(4.1, c(-4 + 7,
                        log(5),
                        sqrt(4) + 3 * exp(3),
                        2 - 0.1/2,
                        2 - 1,
                        0,
                        0))
                 )


    r <- data.frame(rfitness(3))
    r <- r[, c(2, 3, 1, 4)]
    r <- r[c(1, 3, 4, 2, 5, 6), ]
	
	colnames(r)[which(colnames(r) == "Birth")] <- "Fitness"
    r$Fitness <- df3x$Fitness[c(1, 2, 3, 4, 5, 6)]
    df3m <- r

    df3m2 <- data.frame(rfitness(3))[, c(2, 3, 1, 4)][1:6, ]
	colnames(df3m2)[which(colnames(df3m2) == "Birth")] <- "Fitness"
	
    df3m2$Fitness <- df3x$Fitness[c(1, 4, 2, 3, 5, 6)]
    
    ## Note the columns are sorted as per message
    suppressWarnings(adf3m <- allFitnessEffects(genotFitness = df3m,
                               frequencyDependentFitness = TRUE))
    suppressWarnings(adf3m2 <- allFitnessEffects(genotFitness = df3m2,
                                frequencyDependentFitness = TRUE))

    o1 <- evalAllGenotypes(adf3m, spPopSizes = spPopSizes)

    o2 <- evalAllGenotypes(adf3m2, spPopSizes = spPopSizes)
    expect_equal(o1, o2)
    expect_equal(o1, ggx)
})

test_that("should run", {
    r <- data.frame(Genotype = c("WT", "A", "B", "A, B"), 
                    Fitness = c("10*f_", 
                                "10*f_1", 
                                "50*f_2", 
                                "200*(f_1 + f_2) + 50*f_1_2"))

    suppressWarnings(afe <- allFitnessEffects(genotFitness = r, 
                             frequencyDependentFitness = TRUE,
                             frequencyType = "rel"))

    evalAllGenotypes(afe,
                     spPopSizes = c(WT = 2500, A = 2000,
                                    B = 5500, "A, B" = 700))   


    gffd0 <- data.frame(
        Genotype = c(
            "A", "A, B",
            "C", "C, D", "C, E"),
        Fitness = c(
            "1.3",
            "1.4",
            "1.4",
            "1.1 + 0.7*((f_1 + f_A_B) > 0.3)",
            "1.2 + sqrt(f_A + f_C + f_C_D)"))

    suppressWarnings(afd0 <- allFitnessEffects(genotFitness = gffd0,
                              frequencyDependentFitness = TRUE))

    sp <- 1:5
    names(sp) <- c("A", "C", "A, B", "C, D", "C, E")
    eag0 <- evalAllGenotypes(afd0, spPopSizes = sp)

})


cat(paste("\n Ending FDF-small-fitness-specs", date(), "\n"))


cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)

