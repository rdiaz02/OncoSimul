## Root and no root issues in modules
## Works OK
fnme <- allFitnessEffects(epistasis = c("A" = 0.1,
                                        "B" = 0.2),
                          geneToModule = c("A" = "a1, a2",
                                           "B" = "b1"))
evalAllGenotypes(fnme, order = FALSE, addwt = TRUE)



## Clarify the root thing in vignette: here we do not need it. So remove
## from all similar examples. Likewise for help. And say: "need Root in
## modules if it is in the fitness specification. Othewrise, don't need
## it".

sa <- 0.1
sb <- -0.2
sab <- 0.25
sac <- -0.1
sbc <- 0.25
sv2 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                                       "A : -B" = sa,
                                       "A : C" = sac,
                                       "A:B" = sab,
                                       "-A:B:C" = sbc),
                         geneToModule = c(
                             "Root" = "Root",
                             "A" = "a1, a2",
                             "B" = "b",
                             "C" = "c"))
eagr <- evalAllGenotypes(sv2, order = FALSE, addwt = TRUE)


sa <- 0.1
sb <- -0.2
sab <- 0.25
sac <- -0.1
sbc <- 0.25
sv2 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                                       "A : -B" = sa,
                                       "A : C" = sac,
                                       "A:B" = sab,
                                       "-A:B:C" = sbc),
                         geneToModule = c(
                             "A" = "a1, a2",
                             "B" = "b",
                             "C" = "c"))
eagnr <- evalAllGenotypes(sv2, order = FALSE, addwt = TRUE)

expect_identical(eagr, eagnr)


## but here, we need it

s <- 0.2
## OK
m1 <- allFitnessEffects(data.frame(
    parent = c("Root", "A"),
    child = c("A", "B"),
    s = s,
    sh = -1,
    typeDep = "OR"),
    geneToModule = c("Root" = "Root",
                     "A" = "a1, a2",
                     "B" = "b1"))
## But of course this fails
m2 <- allFitnessEffects(data.frame(
    parent = c("Root", "A"),
    child = c("A", "B"),
    s = s,
    sh = -1,
    typeDep = "OR"),
    geneToModule = c("A" = "a1, a2",
                     "B" = "b1"))



