seed <- 1

while(TRUE) {
    seed <- seed + 1
    set.seed(seed)
    cat("\n seed = ", seed, "\n")
    cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                                child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                      s = 0.1,
                                sh = -0.9,
                      typeDep = "MN")
    cbn1 <- allFitnessEffects(cs)
    o1 <- oncoSimulIndiv(cbn1, detectionSize = 1e4,
                         onlyCancer = TRUE,
                         max.num.tries = 5000,
                         sampleEvery = 0.03, keepEvery = 1)
    o4 <- oncoSimulPop(4,
                       cbn1, 
                       detectionSize = 1e4,
                       onlyCancer = TRUE,
                       mc.cores = 2,
                       max.num.tries = 5000,
                       sampleEvery = 0.03, keepEvery = 1)
     expect_message(samplePop(o4, typeSample = "single-nowt",
                             popSizeSample = c(9000, 9000, 8500, 9000)),
                   "Subjects by Genes matrix of 4 subjects and 6 genes")
    o41 <- oncoSimulPop(4,
                        cbn1,
                        initSize = 2e3,
                        detectionSize = 1e3,
                        onlyCancer = TRUE,
                        mc.cores = 2,
                        max.num.tries = 5000,
                        sampleEvery = 0.03, keepEvery = 1)
    expect_message(samplePop(o41, typeSample = "single-nowt",
                             popSizeSample = c(900, 800, 850, 900)),
                   "Subjects by Genes matrix of 4 subjects and 6 genes")
    expect_warning(samplePop(o41, typeSample = "single-nowt",
                             popSizeSample = c(900, 800, 850, 900)),
                   "No non-WT clone with required popSize or at required time")
    o91 <- oncoSimulPop(4,
                       cbn1, 
                       detectionSize = 1e5,
                       onlyCancer = TRUE,
                       mc.cores = 2,
                       max.num.tries = 5000,
                       sampleEvery = 0.03, keepEvery = 1)
     expect_message(samplePop(o91, typeSample = "single-nowt",
                             popSizeSample = c(9000, 9000, 8500, 9000)),
                   "Subjects by Genes matrix of 4 subjects and 6 genes")
}
