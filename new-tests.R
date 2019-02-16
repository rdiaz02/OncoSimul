library(OncoSimulR)
library(testthat)

df1 <- data.frame(Genotype = c("A", "B, C"), Fitness = c(1.3, 2),
                  stringsAsFactors = FALSE)
expect_warning(OncoSimulR:::allGenotypes_to_matrix(df1),
               "No WT genotype. Setting its fitness to 1.", fixed = TRUE)

df2 <- data.frame(Genotype = c("WT", "A", "B, C"), Fitness = c(5, 1.3, 2))
expect_warning(OncoSimulR:::allGenotypes_to_matrix(df2),
               "First column of genotype fitness is a factor.",
               fixed = TRUE)




df1 <- data.frame(Genotype = c("A", "B, C"), Fitness = c(1.3, 2),
                  stringsAsFactors = FALSE)
expect_warning(OncoSimulR:::to_genotFitness_std(df1),
               "No WT genotype. Setting its fitness to 1.", fixed = TRUE)


df2 <- data.frame(Genotype = c("WT", "A", "B, C"), Fitness = c(5, 1.3, 2))
expect_warning(OncoSimulR:::to_genotFitness_std(df2),
               "First column of genotype fitness is a factor.",
               fixed = TRUE)

df3 <- data.frame(Genotype = c("WT", "A", "B, C"), Fitness = c(1.3, 2, 0),
                  stringsAsFactors = FALSE)

m3 <- rbind(c(0, 0, 0, 1.3), c(1, 0, 0, 2.0))
colnames(m3) <- c("A", "B", "C", "Fitness")
m4 <- rbind(c(0, 0, 0, 1.3), c(1, 0, 0, 2.0), c(0, 1, 1, 0))
colnames(m4) <- c("A", "B", "C", "Fitness")
expect_equal(OncoSimulR:::to_genotFitness_std(df3), m3)
expect_equal(OncoSimulR:::to_genotFitness_std(df3, simplify = FALSE), m4)

for(i in 1:10) {
    rxx <- rfitness(7)
    expect_equal(
        allFitnessEffects(genotFitness = rxx)$fitnessLandscape,
        OncoSimulR:::to_genotFitness_std(rxx))
}


for(i in 1:10) {
    ng <- 7
    rxx <- rfitness(ng)
    cnn <-
        replicate(ng,
                  paste(sample(c(LETTERS,
                                 letters,
                                 0:9), 5, replace = TRUE),
                        collapse = ""))
    if(any(duplicated(cnn))) cnn <- LETTERS[1:ng]
    colnames(rxx)[1:ng] <- cnn
    fex <- allFitnessEffects(genotFitness = rxx)
    gn <- OncoSimulR:::allNamedGenes(fex)
    m1x <- data.frame(Gene = gtools::mixedsort(cnn),
                      GeneNumID = 1:ng,
                      stringsAsFactors = FALSE)
    expect_identical(m1x, gn)
}

## FIXME: to do
## taken an rT, convert to fitness landscape, and verify we get
## same fitnesses

## this should work!
test_that("drv names OK", {
    rxx <- rfitness(5)
    expect_silent(allFitnessEffects(genotFitness = rxx, drvNames = LETTERS[1:4]))
    })

## I think this is already tested
## rxx <- rfitness(3)
## allFitnessEffects(genotFitness = rxx,
##                   geneToModule = c("Root" = "Root",
##                                    "A" = "a1, a2",
##                                    "B" = "b1",
##                                    "C" = "c1"))

## Make sure warning if using Bozic
test_that("Bozic and fitness landscape spec", {
    rxx <- rfitness(7)
    expect_warning(oncoSimulIndiv(
        allFitnessEffects(genotFitness = rxx),
        model = "Bozic", initSize = 5000,
        onlyCancer = FALSE,
        finalTime = 10,
        verbosity = 0),
        "Bozic model passing a fitness landscape will not work for now",
        fixed = TRUE)
    rm(rxx)
})



test_that("fitness evaluation what we expect", {
    for(i in 1:10) {
        rxx <- rfitness(5)
        ## allFitnessEffects(genotFitness = rxx)
        eag <- evalAllGenotypes(allFitnessEffects(genotFitness = rxx),
                                addwt = TRUE)
        rxxf <- rxx[, "Fitness"]
        rxxf[rxxf <= 1e-09] <- 0
        expect_equal(rxxf, eag[, "Fitness"])
    }
})


## set.seed(1)
## rxx <- rfitness(5)
## rxx[2, 6] <- 2
## simul1 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx,
##                                            drvNames = LETTERS[1:5]),
##                          model = "Exp", initSize = 5000,
##                          onlyCancer = FALSE,
##                          finalTime = 300,
##                          verbosity = 3)
## summary(simul1)


## set.seed(1)
## rxx <- rfitness(5)
## rxx[2, 6] <- 2
## simul1 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx,
##                                            drvNames = LETTERS[1:5]),
##                          model = "Exp", initSize = 5000,
##                          onlyCancer = FALSE,
##                          finalTime = 1000,
##                          verbosity = 0)
## summary(simul1)





## rxx <- rfitness(7)
## simul1 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rxx,
##                                            drvNames = LETTERS[1:7]),
##                          model = "Exp", initSize = 5000,
##                          onlyCancer = FALSE,
##                          finalTime = 1000,
##                          verbosity = 0)
## summary(simul1)


## n: number of genes
dag_fitness <- function(n) {
    s1min <- 0.1
    s1max <- 0.7
    dummys1 <- 0.5 ## to put something below, to begin with

    rt <- simOGraph(n, out ="rT", geneNames = LETTERS[1:n],
                    s = dummys1, sh = -99)
    ## nparents = sample(2:5, 1),
    ##                h = sample(2:5, 1))

    ## Make sure we get variation 
    uc <- unique(rt$child)
    s1v <- runif(uc, s1min, s1max)
    names(s1v) <- uc
    rt$s <- s1v[rt$child]
    
    rtf <- evalAllGenotypes(allFitnessEffects(rt), addwt = TRUE)
    fl <- OncoSimulR:::allGenotypes_to_matrix(rtf)
    fl[fl[, "Fitness"] == 0, "Fitness"] <- 1e-9
    return(list(rt = rt, fl = fl))
}


## rtfl <- dag_fitness(5)

## set.seed(2)
## s1 <- oncoSimulIndiv(allFitnessEffects(rtfl$rt))
## set.seed(2)
## s2 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rtfl$fl))
## summary(s1)
## summary(s2)


test_that("rt and fl specifications are the same", {
    ## We test that passing a DAG as a DAG of restrictions or as its
    ## fitness landscape give identical output
    tests <- 10
    ng <- 7
    for(i in 1:tests) {
        rtfl <- dag_fitness(ng)
        is <- round(runif(1) * 100)
        set.seed(is)
        s1 <- oncoSimulIndiv(allFitnessEffects(rtfl$rt))
        set.seed(is)
        s2 <- oncoSimulIndiv(allFitnessEffects(genotFitness = rtfl$fl))
        expect_identical(s1$pops.by.time, s2$pops.by.time)
        print(summary(s1))
        expect_identical(s1[1:length(s1)], s2[1:length(s2)])
        ## they differ in the call attribute, of course
        ## adding a package for this is an overkill
        ## expect_true(compare::compare(s1, s2, ignoreAttrs = TRUE)$result)
    }
})


## FIXME: some tests with mutator, etc?

## NOTE the BREAKING changes!!! missing genotypes set to 0




## catching the label bug
o3 <- allFitnessEffects(orderEffects = c(
                                          "M > D > F" = 0.99,
                                          "D > M > F" = 0.2,
                                          "D > M"     = 0.1,
                                          "M > D"     = 0.9),
                                      noIntGenes = c("u" = 0.01, "z" = 0.01, "w" = 0.02),
                                      geneToModule =
                                          c("Root" = "Root",
                                            "M" = "m",
                                            "F" = "f",
                                            "D" = "d") )
              tmp <- oncoSimulIndiv(o3, model = "McFL",
                      mu = 5e-5, finalTime = 500,
                      detectionDrivers = 3,
                      sampleEvery = 0.03,
                      keepEvery = 1,
                      onlyCancer = FALSE,
                      initSize = 1000,
                      keepPhylog = TRUE
                     , initMutant = c("d > m > w")
                      )
tmp$GenotypesLabels
