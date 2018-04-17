inittime <- Sys.time()
cat(paste("\n Starting test.Z-rfitness-landscape at", date(), "\n"))
test_that("exercise accessible genotypes code", {
    ## Because of compiler differences, this might only test the relevant
    ## code in Linux. But conditions will be true in all systems.
    
    RNGkind("Mersenne-Twister")

    set.seed(3)
    r1 <- rfitness(4, reference = rep(0, 4), c = 2.1,
                   min_accessible_genotypes = 1)
    expect_true(max(r1[-1, "Fitness"]) >= 1.0)
    
    set.seed(3)
    r1 <- rfitness(4, reference = rep(0, 4), c = 2.1,
                   min_accessible_genotypes = 1,
                   accessible_th = 0.06)
    expect_true(max(r1[-1, "Fitness"]) >= 1.06)

    set.seed(3)
    r2 <- rfitness(4, reference = rep(0, 4), c = 2.1,
                   min_accessible_genotypes = 2,
                   accessible_th = 0.06)
    expect_true(max(r2[-1, "Fitness"]) >= 1.06)

    ## Exercising plotting code
    plot(r2)
    plot(r2, only_accessible = TRUE)
    plot(r2, only_accessible = TRUE, accessible_th = 0.06)
    ## this only show WT
    plot(r2, only_accessible = TRUE, accessible_th = 0.2)
   
}) 

set.seed(NULL)

cat(paste("\n Ending test.Z-rfitness-landscape at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
