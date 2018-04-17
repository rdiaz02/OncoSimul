inittime <- Sys.time()
cat(paste("\n Starting wide2long tests", date(), "\n"))

## RNGkind("Mersenne-Twister")
test_that("wide2long, two basic cases, minimal testing",
{
    data(examplePosets)
    ## An object of class oncosimul
    p705 <- examplePosets[["p705"]]
    p1 <- oncoSimulIndiv(p705, sampleEvery = 0.025,
                         keepEvery = 5, onlyCancer = FALSE)
    class(p1)
    lp1 <- OncoSimulWide2Long(p1)
    head(lp1)
    summary(lp1)
    tt1 <- table(table(lp1$Time))
    expect_true(ncol(lp1) == 4)
    ## expect_true(any(is.na(lp1))) in very rare cases we could have full
    ## data
    expect_true(length(tt1) == 1)
    expect_true(tt1 == length(unique(p1$pops.by.time[, 1])))
    expect_true(names(tt1) == as.character(length(unique(lp1$Genotype))))
    
    ## An object of class oncosimul2
    data(examplesFitnessEffects)
    sm <-  oncoSimulIndiv(examplesFitnessEffects$cbn1,
                          model = "McFL", 
                          mu = 5e-7,
                          detectionSize = 1e8, 
                          detectionDrivers = 2,
                          sampleEvery = 0.025,
                          keepEvery = 5,
                          initSize = 2000,
                          onlyCancer = FALSE)
    class(sm)
    lsm <- OncoSimulWide2Long(sm)
    head(lsm)
    summary(lsm)
    tt2 <- table(table(lsm$Time))
    expect_true(ncol(lsm) == 4)
    ## expect_true(any(is.na(lsm))) ## also, rarely we might have full data.
    expect_true(length(tt2) == 1)
    expect_true(tt2 == length(unique(sm$pops.by.time[, 1])))
    expect_true(names(tt2) == as.character(length(unique(lsm$Genotype))))
})

cat(paste("\n Ending wide2long tests", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
