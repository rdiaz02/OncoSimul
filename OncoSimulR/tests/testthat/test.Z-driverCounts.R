inittime <- Sys.time()
cat(paste("\n Starting test.Z-driver-counts at", date(), "\n"))
date()
test_that("Assertion is correct when nothing returned",{
    ## Fixing the seed here isn't really needed. These are just a few
    ## specific cases in Linux, but other compilers will give different
    ## runs. None should fail, though.
    RNGkind("L'Ecuyer-CMRG") 
    set.seed(13)
    oi <- allFitnessEffects(orderEffects =
                                c("F > D" = -0.3, "D > F" = 0.4),
                            noIntGenes = runif(5, 0.01, 0.06),
                            geneToModule =
                                c("Root" = "Root",
                                  "F" = "f1, f2, f3",
                                  "D" = "d1, d2"),
                            drvNames = c("d1", "d2", "f1", "f2", "f3"))
    set.seed(13)
    expect_message(ou <- oncoSimulSample(1, 
                                  oi,
                                  sampleEvery = 0.03,
                                  onlyCancer = FALSE,
                                  model = "Bozic",
                                  mutationPropGrowth = TRUE,
                                  seed = NULL),
                  "Subjects by Genes", fixed = TRUE)



    RNGkind("Mersenne-Twister") 
    set.seed(13)
    oi <- allFitnessEffects(orderEffects =
                                c("F > D" = -0.3, "D > F" = 0.4),
                            noIntGenes = runif(5, 0.01, 0.06),
                            geneToModule =
                                c("Root" = "Root",
                                  "F" = "f1, f2, f3",
                                  "D" = "d1, d2"),
                            drvNames = c("d1", "d2", "f1", "f2", "f3"))
    set.seed(13)
    expect_message(ou2 <- oncoSimulSample(1, 
                    oi,
                    sampleEvery = 0.03,
                    onlyCancer = FALSE,
                    model = "Bozic",
                    mutationPropGrowth = TRUE,
                    seed = NULL),
                  "Subjects by Genes", fixed = TRUE)

} )
date()


## Just to check no crashes and to show output
date()
test_that("driverCounts: a run that used to cause crashes", {
    RNGkind("Mersenne-Twister")
    set.seed(1)
    ni <- runif(10, min = -0.01, max = 0.1)
    names(ni) <- paste0("g", 2:11)
    fe31 <- allFitnessEffects(noIntGenes = ni,
                              drvNames = "g2") 
    set.seed(1)
    mue11 <- oncoSimulIndiv(fe31, 
                            mu = 1e-6,
                            initSize = 1e5,
                            model = "McFL",
                            detectionSize = 5e6,
                            finalTime = 5, 
                            onlyCancer = FALSE)
    cat("\n ... output from mue11\n")
    print(mue11)
})

set.seed(NULL)

cat(paste("\n Ending test.Z-driver-counts at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
