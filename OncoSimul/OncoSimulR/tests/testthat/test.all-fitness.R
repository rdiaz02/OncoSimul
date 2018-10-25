inittime <- Sys.time()
cat(paste("\n Starting all fitness at", date()))
## RNGkind("Mersenne-Twister") ## we have some examples below with random
## genotype generation. We leave the default. We use L'Ecuyer when using
## mclapply and that file could run before this one


test_that("WT named genes give a warning", {
    m1 <- cbind(WT = c(0, 1), B = c(0, 1), Fitness = c(1, 1e-8))
    expect_warning(s1 <- oncoSimulIndiv(allFitnessEffects(genotFitness = m1),
                                       detectionSize = 1, initSize = 100,
                                       keepPhylog = TRUE),
                   "A gene is named WT", fixed = TRUE)
    fee <- allFitnessEffects(epistasis = c("A" = 0.3,
                                           "B" = 0.5),
                             geneToModule = c("Root" = "Root",
                                              "A" = "WT, a2",
                                              "B" = "b1"))
    expect_warning(s1 <- oncoSimulIndiv(fee,
                                       detectionSize = 1, initSize = 100,
                                       keepPhylog = TRUE),
                   "A gene is named WT", fixed = TRUE)
    fee <- allFitnessEffects(epistasis = c("WT" = 0.3,
                                           "B" = 0.5))
    expect_warning(s1 <- oncoSimulIndiv(fee,
                                       detectionSize = 1, initSize = 100,
                                       keepPhylog = TRUE),
                   "A gene is named WT", fixed = TRUE)
})

test_that("Root name in module table or not", {
    expect_silent(fee <- allFitnessEffects(epistasis = c("A" = 0.3,
                                           "B" = 0.5),
                             geneToModule = c("Root" = "Root",
                                              "A" = "a1, a2",
                                              "B" = "b1")))
    expect_silent(fee2 <- allFitnessEffects(epistasis = c("A" = 0.3,
                                           "B" = 0.5),
                             geneToModule = c("A" = "a1, a2",
                                              "B" = "b1")))
    expect_identical(evalAllGenotypes(fee, order = TRUE), evalAllGenotypes(fee2, order = TRUE))
    expect_error(allFitnessEffects(epistasis = c("A" = 0.3,
                                                 "B" = 0.5),
                                   geneToModule = c(
                                       "Root" = "0",
                                       "A" = "a1, a2",
                                       "B" = "b1")),
                 "The name Root is in the module table, but is not of Root",
                 fixed = TRUE)
    expect_error(allFitnessEffects(epistasis = c("A" = 0.3,
                                                 "B" = 0.5),
                                   geneToModule = c(
                                       "A" = "a1, a2",
                                       "B" = "b1",
                                       "Root" = "0")),
                 "If the name Root is in the module table, it must be the first",
                 fixed = TRUE)
     expect_error(allFitnessEffects(epistasis = c("A" = 0.3,
                                                 "B" = 0.5),
                                   geneToModule = c(
                                       "A" = "a1, a2",
                                       "B" = "b1",
                                       "Root" = "Root")),
                 "If Root is in the module table, it must be the first",
                 fixed = TRUE)
    expect_error(allFitnessEffects(epistasis = c("A" = 0.3,
                                                 "B" = 0.5),
                                   geneToModule = c(
                                       "0" = "Root",
                                       "A" = "a1, a2",
                                       "B" = "b1")),
                 "Some values in geneToModule not present in any of",
                 fixed = TRUE)
    expect_error(allFitnessEffects(epistasis = c("A" = 0.3,
                                                 "B" = 0.5),
                                   geneToModule = c(
                                       "B" = "Root",
                                       "A" = "a1, a2",
                                       "Root" = "b1")),
                 "If the name Root is in the module table, it must be the first",
                 fixed = TRUE)
    expect_error(allFitnessEffects(epistasis = c("A" = 0.3,
                                                 "B" = 0.5),
                                   geneToModule = c(
                                       "Root" = "b1",
                                       "A" = "a1, a2",
                                       "B" = "Root")),
                 "If Root is in the module table, it must be the first",
                 fixed = TRUE)
       expect_error(allFitnessEffects(epistasis = c("A" = 0.3,
                                                 "B" = 0.5),
                                   geneToModule = c(
                                       "A" = "a1, a2",
                                       "Root" = "b1",
                                       "B" = "Root")),
                 "If Root is in the module table, it must be the first",
                 fixed = TRUE)
       expect_error(allFitnessEffects(epistasis = c("A" = 0.3,
                                                 "B" = 0.5),
                                   geneToModule = c(
                                       "B" = "Root",
                                       "A" = "a1, a2",
                                       "B" = "c1")),
                 "Root is in the module table, but with a different name",
                 fixed = TRUE)
})




test_that("Bauer example: correct number of fitness classes", {
    sd <- 0.1
    sdp <- 0.15
    sp <- 0.05
    bauer <- data.frame(parent = c("Root", rep("p", 5)),
                        child = c("p", paste0("s", 1:5)),
                        s = c(sd, rep(sdp, 5)),
                        sh = c(0, rep(sp, 5)),
                        typeDep = "MN")
    b1 <- evalAllGenotypes(allFitnessEffects(bauer), order = FALSE)
    b2 <- evalAllGenotypes(allFitnessEffects(bauer), order = TRUE, max = 2000)
    expect_equal(length(unique(b1$Fitness)), 11)
    expect_equal(length(unique(b2$Fitness)), 11)
} )



test_that("Bauer example: identical values fitness classes, unorder and ord", {
    sd <- 0.1
    sdp <- 0.15
    sp <- 0.05
    bauer <- data.frame(parent = c("Root", rep("p", 5)),
                        child = c("p", paste0("s", 1:5)),
                        s = c(sd, rep(sdp, 5)),
                        sh = c(0, rep(sp, 5)),
                        typeDep = "MN")
    b1 <- evalAllGenotypes(allFitnessEffects(bauer), order = FALSE)
    b2 <- evalAllGenotypes(allFitnessEffects(bauer), order = TRUE, max = 2000)
    expect_equal(unique(b1$Fitness), unique(b2$Fitness))
} )


test_that("Bauer example: identical values fitness classes, rename", {
    sd <- 0.1
    sdp <- 0.15
    sp <- 0.05
    bauer <- data.frame(parent = c("Root", rep("p", 5)),
                        child = c("p", paste0("s", 1:5)),
                        s = c(sd, rep(sdp, 5)),
                        sh = c(0, rep(sp, 5)),
                        typeDep = "MN")
    b1 <- evalAllGenotypes(allFitnessEffects(bauer), order = FALSE)
    bauer3 <- data.frame(parent = c("Root", rep("u", 5)),
                         child = c("u", paste0("s", 1:5)),
                         s = c(sd, rep(sdp, 5)),
                         sh = c(0, rep(sp, 5)),
                         typeDep = "MN")
    b3 <- evalAllGenotypes(allFitnessEffects(bauer), order = TRUE, max = 2000)
    expect_equal(unique(b1$Fitness), unique(b3$Fitness))
} )


test_that("Bauer example: identical values fitness classes, diff. order", {
    sd <- 0.1
    sdp <- 0.15
    sp <- 0.05
    bauer <- data.frame(parent = c("Root", rep("p", 5)),
                        child = c("p", paste0("s", 1:5)),
                        s = c(sd, rep(sdp, 5)),
                        sh = c(0, rep(sp, 5)),
                        typeDep = "MN")
    b1 <- evalAllGenotypes(allFitnessEffects(bauer), order = FALSE)
    bauer3 <- data.frame(parent = c(rep("u", 5), "Root"),
                         child = c(paste0("s", 1:5), "u"),
                         s = c(sd, rep(sdp, 5)),
                         sh = c(0, rep(sp, 5)),
                         typeDep = "MN")
    b3 <- evalAllGenotypes(allFitnessEffects(bauer), order = TRUE, max = 2000)
    expect_equal(unique(b1$Fitness), unique(b3$Fitness))
} )


### Order effects


test_that("Order effects, entry of labels and separation", {
    o1 <- evalAllGenotypes(allFitnessEffects(
        orderEffects = c("d>f" = 0.4, "f > d" = -0.3) ),
        order = TRUE)
    o2 <- evalAllGenotypes(allFitnessEffects(
        orderEffects = c("f > d" = -0.3, "d > f" = 0.4) ),
        order = TRUE)
    o1s <- o1[order(o1$Genotype),   ]
    o2s <- o2[order(o2$Genotype),   ]
    expect_equal(o1s, o2s)
})

test_that("Order effects, modules 1", {
    ofe1 <- allFitnessEffects(orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
                              geneToModule =
                                  c("Root" = "Root",
                                    "F" = "f1, f2",
                                    "D" = "d1, d2") )
    ag <- evalAllGenotypes(ofe1, order = TRUE)
    expect_true(all.equal(ag[c(17, 39, 19, 29), "Fitness"], c(1.4, 0.7, 1.4, 0.7)))
    expect_true(all.equal(ag[c(43, 44), "Fitness"], c(1.4, 1.4)))
    expect_true(all(ag[41:52, "Fitness"] == 1.4))
    expect_true(all(ag[53:64, "Fitness"] == 0.7))
})

test_that("Order effects, modules 2", {
    ofe2 <- allFitnessEffects(orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
                              geneToModule =
                                  c("Root" = "Root",
                                    "F" = "f1, f2, f3",
                                    "D" = "d1, d2") )
    ag2 <- evalAllGenotypes(ofe2, max = 326, order = TRUE)
    oe <- c(grep("^f.*d.*", ag2[, 1]), grep("^d.*f.*", ag2[, 1]))
    expect_true(all(ag2[grep("^d.*f.*", ag2[, 1]), "Fitness"] == 1.4))
    expect_true(all(ag2[grep("^f.*d.*", ag2[, 1]), "Fitness"] == 0.7))
    expect_true(all(ag2[-oe, "Fitness"] ==  1))
})



test_that("Order effects, twisted module names", {
    o1 <- evalAllGenotypes(allFitnessEffects(
      orderEffects = c("F > D" = -0.3, "D > F" = 0.4),  
        geneToModule = c("Root" = "Root", "F" = "d", "D" = "f")), order = TRUE)
    o2 <- evalAllGenotypes(allFitnessEffects(
      orderEffects = c("D>F" = 0.4, "F >D" = -0.3),  
        geneToModule = c("Root" = "Root", "F" = "d", "D" = "f")), order = TRUE)
    o1s <- o1[order(o1$Genotype),   ]
    o2s <- o2[order(o2$Genotype),   ]
    expect_equal(o1s, o2s)
    expect_true(all.equal(o1[, "Fitness"], c(1, 1, 0.7, 1.4)))
})



test_that("Order effects, three-gene-orders and modules 1", {
    o3 <- allFitnessEffects(orderEffects = c(
                                  "F > D > M" = -0.3,
                                  "D > F > M" = 0.4,
                                  "D > M > F" = 0.2,
                                  "D > M"     = 0.1,
                                  "M > D"     = 0.5
    ),
                              geneToModule =
                                  c("Root" = "Root",
                                    "M" = "m",
                                    "F" = "f",
                                    "D" = "d") )
    ag <- evalAllGenotypes(o3, order = TRUE)
    expect_true(all.equal(ag[, "Fitness"],
                        c(rep(1, 4),
                          1.1,
                          1, 1,
                          1.5,
                          1,
                          1.4 * 1.1,
                          1.2 * 1.1,
                          0.7 * 1.1,
                          rep(1.5, 3))))
})



################ No interaction


test_that("No interaction genes, 1", {
    ai1 <- evalAllGenotypes(allFitnessEffects(
        noIntGenes = c(0.05, -.2, .1)), order = FALSE)
    
    ai2 <- evalAllGenotypes(allFitnessEffects(
        noIntGenes = c(0.05, -.2, .1)), order = TRUE)

    
    expect_true(all.equal(ai1[, "Fitness"],  c( (1 + .05), (1 - .2), (1 + .1),
       (1 + .05) * (1 - .2),
       (1 + .05) * (1 + .1),
       (1 - .2) * (1 + .1),
       (1 + .05) * (1 - .2) * (1 + .1))))

    expect_true(all.equal(ai2[, "Fitness"],  c((1 + .05), (1 - .2), (1 + .1),
                           1.05 * .8, 1.05 * 1.1, .8 * 1.05, .8 * 1.1,
                           1.05 * 1.1, 1.1 * .8,
                           rep(1.05 * .8 * 1.1, 6) )))
})


## to the above, add gene names

test_that("No interaction genes, 2", {
    
    ai3 <- evalAllGenotypes(allFitnessEffects(
        noIntGenes = c("a" = 0.05, "b" = -.2, "c" = .1)), order = FALSE)
    
    ai4 <- evalAllGenotypes(allFitnessEffects(
        noIntGenes = c("a" = 0.05, "b" = -.2, "c" = .1)), order = TRUE)
    
    expect_true(all.equal(ai3[, "Fitness"],  c( (1 + .05), (1 - .2), (1 + .1),
       (1 + .05) * (1 - .2),
       (1 + .05) * (1 + .1),
       (1 - .2) * (1 + .1),
       (1 + .05) * (1 - .2) * (1 + .1))))

    expect_true(all.equal(ai4[, "Fitness"], c((1 + .05), (1 - .2), (1 + .1),
                           1.05 * .8, 1.05 * 1.1, .8 * 1.05, .8 * 1.1,
                           1.05 * 1.1, 1.1 * .8,
                           rep(1.05 * .8 * 1.1, 6) )))

})


## to the above, resort gene names

test_that("No interaction genes, 3", {
    
    ai3 <- evalAllGenotypes(allFitnessEffects(
        noIntGenes = c("m" = 0.05, "b" = -.2, "f" = .1)), order = FALSE)
    
    ai4 <- evalAllGenotypes(allFitnessEffects(
        noIntGenes = c("m" = 0.05, "b" = -.2, "f" = .1)), order = TRUE)
    
    expect_true(all.equal(ai3[, "Fitness"],  c( (1 + .05), (1 - .2), (1 + .1),
       (1 + .05) * (1 - .2),
       (1 + .05) * (1 + .1),
       (1 - .2) * (1 + .1),
       (1 + .05) * (1 - .2) * (1 + .1))))

    expect_true(all.equal(ai4[, "Fitness"],  c((1 + .05), (1 - .2), (1 + .1),
                           1.05 * .8, 1.05 * 1.1, .8 * 1.05, .8 * 1.1,
                           1.05 * 1.1, 1.1 * .8,
                           rep(1.05 * .8 * 1.1, 6) )))

})




test_that("No interaction genes and order effects, 1", {
    foi1 <- allFitnessEffects(
        orderEffects = c("D>B" = -0.3, "B > D" = 0.3),
        noIntGenes = c("A" = 0.05, "C" = -.2, "E" = .1))
    agoi1 <- evalAllGenotypes(foi1,  max = 325, order = TRUE)
    rn <- 1:nrow(agoi1)
    names(rn) <- agoi1[, 1]
    expect_true(all.equal(agoi1[rn[LETTERS[1:5]], "Fitness"], c(1.05, 1, 0.8, 1, 1.1)))
    ## orders that do not involve all. D > A;   B > C;
    expect_true(all(agoi1[grep("^A > [BD]$", names(rn)), "Fitness"] == 1.05))
    expect_true(all(agoi1[grep("^C > [BD]$", names(rn)), "Fitness"] == 0.8))
    expect_true(all(agoi1[grep("^E > [BD]$", names(rn)), "Fitness"] == 1.1))
    expect_true(all(agoi1[grep("^[BD] > A$", names(rn)), "Fitness"] == 1.05))
    expect_true(all(agoi1[grep("^[BD] > C$", names(rn)), "Fitness"] == 0.8))
    expect_true(all(agoi1[grep("^[BD] > E$", names(rn)), "Fitness"] == 1.1))
    expect_true(all.equal(agoi1[230:253, "Fitness"] ,
                          rep((1 - 0.3) * 1.05 * 0.8 * 1.1, 24)))
    expect_true(all.equal(agoi1[c(260:265, 277, 322, 323, 325), "Fitness"] ,
              rep((1 - 0.3) * 1.05 * 0.8 * 1.1, 10)))
    expect_true(all.equal(agoi1[c(206:229, 254:259, 266:267), "Fitness"] ,
              rep((1 + 0.3) * 1.05 * 0.8 * 1.1, 32)))
    ## some of four, one of which either D or B. 
    expect_true(all.equal(agoi1[c(203:205, 191), "Fitness"] ,
              rep(1.05 * 0.8 * 1.1, 4)))
    ##  a few of three, A, C, and ether D or B
    expect_true(all.equal(agoi1[c(42, 45, 30, 33), "Fitness"] ,
              rep(1.05 * 0.8, 4)))
})


## like above, but change order of names
test_that("No interaction genes and order effects, 2", {
    foi1 <- allFitnessEffects(
        orderEffects = c("D>B" = -0.3, "B > D" = 0.3),
        noIntGenes = c("M" = 0.05, "A" = -.2, "J" = .1))
    agoi1 <- evalAllGenotypes(foi1,  max = 325, order = TRUE)
    rn <- 1:nrow(agoi1)
    names(rn) <- agoi1[, 1]
    expect_true(all.equal(agoi1[rn[c("B", "D", "M", "A", "J")], "Fitness"],
                          c(1, 1, 1.05, 0.8, 1.1)))
    ## orders that do not involve all. D > A;   B > C;
    expect_true(all(agoi1[grep("^M > [BD]$", names(rn)), "Fitness"] == 1.05))
    expect_true(all(agoi1[grep("^A > [BD]$", names(rn)), "Fitness"] == 0.8))
    expect_true(all(agoi1[grep("^J > [BD]$", names(rn)), "Fitness"] == 1.1))
    expect_true(all(agoi1[grep("^[BD] > M$", names(rn)), "Fitness"] == 1.05))
    expect_true(all(agoi1[grep("^[BD] > A$", names(rn)), "Fitness"] == 0.8))
    expect_true(all(agoi1[grep("^[BD] > J$", names(rn)), "Fitness"] == 1.1))
    expect_true(all.equal(agoi1[230:253, "Fitness"] ,
                          rep((1 - 0.3) * 1.05 * 0.8 * 1.1, 24)))
    expect_true(all.equal(agoi1[c(260:265, 277, 322, 323, 325), "Fitness"] ,
              rep((1 - 0.3) * 1.05 * 0.8 * 1.1, 10)))
    expect_true(all.equal(agoi1[c(206:229, 254:259, 266:267), "Fitness"] ,
              rep((1 + 0.3) * 1.05 * 0.8 * 1.1, 32)))
    ## some of four, one of which either D or B. 
    expect_true(all.equal(agoi1[c(203:205, 191), "Fitness"] ,
              rep(1.05 * 0.8 * 1.1, 4)))
    ##  a few of three, A, C, and ether D or B
    expect_true(all.equal(agoi1[c(42, 45, 30, 33), "Fitness"] ,
              rep(1.05 * 0.8, 4)))
})



### synthetic viability and mortality
test_that("synthetic viability, 1", {
    s <- 0.2
    sv <- allFitnessEffects(epistasis = c("-A : B" = -1,
                                "A : -B" = -1,
                                "A:B" = s))
    expect_true( all.equal(
        evalAllGenotypes(sv,
                         order = FALSE, addwt = TRUE)[, "Fitness"] , c(1, 0, 0, 1.2)))
    expect_true( all.equal(
        evalAllGenotypes(sv,
                         order = TRUE, addwt = TRUE)[, "Fitness"], c(1, 0, 0, 1.2, 1.2)))
})


test_that("synthetic viability, with modules, 2", {
    sa <- -0.1
    sb <- -0.2
    sab <- 0.25
    sv2 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                                "A : -B" = sa,
                                 "A:B" = sab),
                             geneToModule = c(
                                 "Root" = "Root",
                                 "A" = "a1, a2",
                                 "B" = "b"))
    expect_true( all.equal(
        evalAllGenotypes(sv2,
                         order = FALSE, addwt = TRUE)[, "Fitness"],
        c(1, 0.9, 0.9, 0.8, 0.9, rep(1.25, 3))))
    expect_true( all.equal(
        evalAllGenotypes(sv2,
                         order = TRUE, addwt = TRUE)[, "Fitness"],
        c(1, 0.9, 0.9, 0.8,
          .90, 1.25, .9, rep(1.25, 9)
          )))
})


test_that("synthetic mortality, 1", {
    sa <- 0.1
    sb <- 0.2
    sab <- -0.8
    sm1 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                                 "A : -B" = sa,
                                 "A:B" = sab))
    expect_true( all.equal(
        evalAllGenotypes(sm1, order = FALSE, addwt = TRUE)[, "Fitness"],
        c(1, 1.1, 1.2, 1 - 0.8)))
    expect_true( all.equal(
        evalAllGenotypes(sm1, order = TRUE, addwt = TRUE)[, "Fitness"], 
            c(1, 1.1, 1.2, 1 - 0.8, 1 - 0.8)))
})


### epistasis

test_that("Epistasis, 1", {
    sa <- 0.2
    sb <- 0.3
    sab <- 0.7
    e2 <- allFitnessEffects(epistasis =
                                c("A: -B" = sa,
                                  "-A:B" = sb,
                                  "A : B" = sab))
    expect_true(all.equal(evalAllGenotypes(e2,
                                     order = FALSE,
                                     addwt = TRUE)[, "Fitness"], 
                    c(1, 1.2, 1.3, 1.7)))
})


test_that("Epistasis, with and without -", {
    sa <- 0.2
    sb <- 0.3
    sab <- 0.7
    s2 <- ((1 + sab)/((1 + sa) * (1 + sb))) - 1
    e2 <- allFitnessEffects(epistasis =
                                c("A: -B" = sa,
                                  "-A:B" = sb,
                                  "A : B" = sab))
    e3 <- allFitnessEffects(epistasis =
                                c("A" = sa,
                                  "B" = sb,
                                  "A : B" = s2))
    expect_true(all.equal(evalAllGenotypes(e2,
                                     order = FALSE,
                                     addwt = TRUE)[, "Fitness"],
                                         c(1, 1.2, 1.3, 1.7)))
    expect_true(all.equal(evalAllGenotypes(e3,
                                     order = FALSE,
                                     addwt = TRUE)[, "Fitness"], 
                                         c(1, 1.2, 1.3, 1.7)))
})


test_that("Epistasis, with and without -, three terms", {
    sa <- 0.1
    sb <- 0.15
    sc <- 0.2
    sab <- 0.3
    sbc <- -0.25
    sabc <- 0.4
    sac <- (1 + sa) * (1 + sc) - 1
    E3 <- allFitnessEffects(epistasis =
                                c("A:-B:-C" = sa,
                                  "-A:B:-C" = sb,
                                  "-A:-B:C" = sc,
                                  "A:B:-C" = sab,
                                  "-A:B:C" = sbc,
                                  "A:-B:C" = sac,
                                  "A : B : C" = sabc)
                            )

    expect_true(all.equal(evalAllGenotypes(E3, order = FALSE, addwt = FALSE)[, "Fitness"], 
                    c(1.1, 1.15, 1.2, 1.3, (1.1 * 1.2), 0.75, 1.4)))

})


test_that("Epistasis, three, with and without -, two alternative specs", {
    sa <- 0.1
    sb <- 0.15
    sc <- 0.2
    sab <- 0.3
    sbc <- -0.25
    sabc <- 0.4
    sac <- (1 + sa) * (1 + sc) - 1
    E3A <- allFitnessEffects(epistasis =
                                 c("A:-B:-C" = sa,
                                   "-A:B:-C" = sb,
                                   "-A:-B:C" = sc,
                                   "A:B:-C" = sab,
                                   "-A:B:C" = sbc,
                                   "A:-B:C" = sac,
                                   "A : B : C" = sabc)
                             )
    expect_true(all.equal(evalAllGenotypes(E3A, order = FALSE, addwt = FALSE)[, "Fitness"], 
                        c(1.1, 1.15, 1.2, 1.3, (1.1 * 1.2), 0.75, 1.4)))
    Sab <- ( (1 + sab)/((1 + sa) * (1 + sb))) - 1
    Sbc <- ( (1 + sbc)/((1 + sb) * (1 + sc))) - 1
    Sabc <- ( (1 + sabc)/( (1 + sa) * (1 + sb) * (1 + sc) * (1 + Sab) * (1 + Sbc) ) ) - 1
    E3B <- allFitnessEffects(epistasis =
                                 c("A" = sa,
                                   "B" = sb,
                                   "C" = sc,
                                   "A:B" = Sab,
                                   "B:C" = Sbc,
                                   ## "A:C" = sac,
                                   "A : B : C" = Sabc)
                             )
    expect_true(all.equal(evalAllGenotypes(E3A, order = FALSE, addwt = FALSE),  
                        evalAllGenotypes(E3B, order = FALSE, addwt = FALSE)))
})




test_that("Epistasis, three, with and without -, two alternative specs, order makes no diff", {
    sa <- 0.1
    sb <- 0.15
    sc <- 0.2
    sab <- 0.3
    sbc <- -0.25
    sabc <- 0.4
    sac <- (1 + sa) * (1 + sc) - 1
    E3A <- allFitnessEffects(epistasis =
                                 c("A:-B:-C" = sa,
                                   "-A:B:-C" = sb,
                                   "-A:-B:C" = sc,
                                   "A:B:-C" = sab,
                                   "-A:B:C" = sbc,
                                   "A:-B:C" = sac,
                                   "A : B : C" = sabc)
                             )
    ge3a <- evalAllGenotypes(E3A, order = FALSE, addwt = FALSE)
    ge3ao <- evalAllGenotypes(E3A, order = TRUE, addwt = FALSE)
    Sab <- ( (1 + sab)/((1 + sa) * (1 + sb))) - 1
    Sbc <- ( (1 + sbc)/((1 + sb) * (1 + sc))) - 1
    Sabc <- ( (1 + sabc)/( (1 + sa) * (1 + sb) * (1 + sc) * (1 + Sab) * (1 + Sbc) ) ) - 1
    E3B <- allFitnessEffects(epistasis =
                                 c("A" = sa,
                                   "B" = sb,
                                   "C" = sc,
                                   "A:B" = Sab,
                                   "B:C" = Sbc,
                                   ## "A:C" = sac,
                                   "A : B : C" = Sabc)
                             )
    ge3b <- evalAllGenotypes(E3A, order = FALSE, addwt = FALSE)
    ge3bo <- evalAllGenotypes(E3A, order = TRUE, addwt = FALSE)
    expect_true(identical(ge3ao, ge3bo))
    ## get names that can be matched against the non-ordered
    nn <- ge3ao[, 1]
    nnn <- unlist(lapply(nn, function(x) paste(sort(unlist(strsplit(x, " > "))),
                                               collapse = ", ")))
    ## Verify all of the same name have same value
    ## Beware this could fail for numerical issues.
    expect_true(all( tapply(ge3ao$Fitness, nnn,
                            function(x) length(unique(x))) == 1))
    ## Is the value identical to the unordered?
    mo <- tapply(ge3ao$Fitness, nnn, mean)
    mu <- tapply(ge3a$Fitness, ge3a[, 1], mean)
    expect_true(all.equal(mo, mu))
})





## posets
test_that("Error if not same sh within child", {
    c1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                     child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                     s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                     sh = c(rep(0, 4), c(-.1, -.2), c(-.05, -.06, -.07)),
                     typeDep = "MN")
    expect_error(allFitnessEffects(c1))
})


### Note for the testing below: I use values far from 0, to easily see
### differences.

## now OK
test_that("Poset, CBN, values and order no changes", {
    c1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                     child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                     s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                     sh = c(rep(0, 4), c(-.1, -.1), rep(-.05, 3)),
                     typeDep = "MN")
    fc1 <- allFitnessEffects(c1)
    gfc1 <- evalAllGenotypes(fc1, order = FALSE)
    gfc1o <- evalAllGenotypes(fc1, order = TRUE, max = 1956)
    expect_true(all.equal(
        gfc1[c(1:21, 22, 28, 41, 44, 56, 63 ) , "Fitness"],
        c(1.01, 1.02, 0.9, 1.03, 1.04, 0.95,
          1.01 * c(1.02, 0.9, 1.03, 1.04, 0.95),
          1.02 * c(0.90, 1.03, 1.04, 0.95),
          0.9 * c(1.03, 1.04, 0.95),
          1.03 * c(1.04, 0.95),
          1.04 * 0.95,
          1.01 * 1.02 * 1.1,
          1.01 * 0.9 * 0.95,
          1.03 * 1.04 * 0.95,
          1.01 * 1.02 * 1.1 * 0.95,
          1.03 * 1.04 * 1.2 * 0.9, ## notice this
          1.01 * 1.02 * 1.03 * 1.04 * 1.1 * 1.2
          )))
    ## Verify all of the same name have same value
    ## Beware this could fail for numerical issues.
    nn <- gfc1o[, 1]
    nnn <- unlist(lapply(nn, function(x) paste(sort(unlist(strsplit(x, " > "))),
                                               collapse = ", ")))
    expect_true(all( tapply(gfc1o$Fitness, nnn,
                            function(x) length(unique(x))) == 1))
    mo <- tapply(gfc1o$Fitness, nnn, mean)
    mu <- tapply(gfc1$Fitness, gfc1[, 1], mean)
    expect_true(all.equal(mo, mu))
    ## type of dep for those from root does not matter
    c1b <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                      child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                      s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                      sh = c(rep(0, 4), c(-.1, -.1), rep(-.05, 3)),
                      typeDep = c("-", "--", "SM", "XMPN", rep("MN",5)))
    fc1b <- allFitnessEffects(c1b)
    gfc1b <- evalAllGenotypes(fc1b, order = FALSE)
    expect_true(all.equal(gfc1, gfc1b))
})


## OR:

test_that("Poset, OR, values and order no changes", {
    s1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                     child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                     s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                     sh = c(rep(0, 4), c(-.1, -.1), rep(-.05, 3)),
                     typeDep = "SM")
    fs1 <- allFitnessEffects(s1)
    gfs1 <- evalAllGenotypes(fs1, order = FALSE)
    gfs1o <- evalAllGenotypes(fs1, order = TRUE, max = 1956)
    expect_true(all.equal(
        gfs1[c(1:21, 22, 28, 41, 44, 56, 63, 39 ) , "Fitness"],
        c(1.01, 1.02, 0.9, 1.03, 1.04, 0.95,
          1.01 * c(1.02, 1.1, 1.03, 1.04, 0.95),
          1.02 * c(1.1, 1.03, 1.04, 0.95),
          0.9 * c(1.03, 1.04, 1.2),
          1.03 * c(1.04, 1.2),
          1.04 * 1.2,
          1.01 * 1.02 * 1.1,
          1.01 * 1.1 * 1.2, ## 28
          1.03 * 1.04 * 1.2, ## 41
          1.01 * 1.02 * 1.1 * 1.2, ## 44
          0.9 * 1.03 * 1.04 * 1.2, ## notice this: 56
          1.01 * 1.02 * 1.03 * 1.04 * 1.1 * 1.2,
          0.9 * 1.03 * 1.2
          )))
    nn <- gfs1o[, 1]
    nnn <- unlist(lapply(nn, function(x) paste(sort(unlist(strsplit(x, " > "))),
                                               collapse = ", ")))
    expect_true(all( tapply(gfs1o$Fitness, nnn,
                            function(x) length(unique(x))) == 1))
    mo <- tapply(gfs1o$Fitness, nnn, mean)
    mu <- tapply(gfs1$Fitness, gfs1[, 1], mean)
    expect_true(all.equal(mo, mu))
    zzz <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                      child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                      s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                      sh = c(rep(0, 4), c(-.1, -.1), rep(-.05, 3)),
                      typeDep = c("-", "--", "SM", "XMPN", rep("SM", 5)))
    zzz <- allFitnessEffects(zzz)
    zzz <- evalAllGenotypes(zzz, order = FALSE)
    expect_true(all.equal(gfs1, zzz))
})




test_that("Poset, XOR, values and order no changes", {
    x1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                     s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                     sh = c(rep(0, 4), c(-.9, -.9), rep(-.95, 3)),
                     typeDep = "XMPN")
    fx1 <- allFitnessEffects(x1)
    gfx1 <- evalAllGenotypes(fx1, order = FALSE)
    gfx1o <- evalAllGenotypes(fx1, order = TRUE, max = 1956)
    expect_true(all.equal(
        gfx1[c(1:21, 22, 28, 41, 44, 56, 63, 39 ) , "Fitness"],
        c(1.01, 1.02, 0.1, 1.03, 1.04, 0.05,
          1.01 * c(1.02, 1.1, 1.03, 1.04, 0.05),
          1.02 * c(1.1, 1.03, 1.04, 0.05),
          0.1 * c(1.03, 1.04, 1.2),
          1.03 * c(1.04, 1.2),
          1.04 * 1.2,
          1.01 * 1.02 * 0.1, ## 22
          1.01 * 1.1 * 1.2, ## 28
          1.03 * 1.04 * 0.05, ## 41
          1.01 * 1.02 * 0.1 * 1.2, ## 44
          0.1 * 1.03 * 1.04 * 0.05, ## notice this: 56
          1.01 * 1.02 * 1.03 * 1.04 * 0.1 * 0.05,
          0.1 * 1.03 * 0.05
          )))
    nn <- gfx1o[, 1]
    nnn <- unlist(lapply(nn, function(x) paste(sort(unlist(strsplit(x, " > "))),
                                               collapse = ", ")))
    expect_true(all( tapply(gfx1o$Fitness, nnn,
                            function(x) length(unique(x))) == 1))
    mo <- tapply(gfx1o$Fitness, nnn, mean)
    mu <- tapply(gfx1$Fitness, gfx1[, 1], mean)
    expect_true(all.equal(mo, mu))
    zzz <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                      child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                      s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                      sh = c(rep(0, 4), c(-.9, -.9), rep(-.95, 3)),
                      typeDep = c("SM", "-", "--", "XMPN", rep("XMPN", 5))) 
    zzz <- allFitnessEffects(zzz)
    zzz <- evalAllGenotypes(zzz, order = FALSE)
    expect_true(all.equal(gfx1, zzz))
})


## this order thing gets boring. A function from now on

fouo <- function(fe) {
    uo <- evalAllGenotypes(fe, order = FALSE, max = 50000)
    oo <- evalAllGenotypes(fe, order = TRUE, max = 50000)
    nn <- oo[, 1]
    nnn <- unlist(lapply(nn, function(x) paste(sort(unlist(strsplit(x, " > "))),
                                               collapse = ", ")))
    expect_true(all( tapply(oo$Fitness, nnn,
                            function(x) length(unique(x))) == 1))
    mo <- tapply(oo$Fitness, nnn, mean)
    mu <- tapply(uo$Fitness, uo[, 1], mean)
    ## expect_true(all.equal(mo, mu))
    return(list(mu, mo))
}

test_that("Poset, all three effects", {
    p3 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c", "f"),
                  child = c("a", "b", "d", "e", "c", "c", "f", "f", "g", "g"),
                  s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                  sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                  typeDep = c(rep("--", 4), 
                      "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    fp3 <- allFitnessEffects(p3)
    gfp3 <- evalAllGenotypes(fp3, order = FALSE)
    expect_true(all.equal(gfp3[c(9, 24, 29, 59, 60, 66, 119, 120, 126, 127),
                               "Fitness"],
                          c(1.01 * 1.1, 1.03 * .05, 1.01 * 1.02 * 0.1, 0.1 * 0.05 * 1.3,
                            1.03 * 1.04 * 1.2, 1.01 * 1.02 * 0.1 * 0.05,
                            0.1 * 1.03 * 1.04 * 1.2 * 1.3,
                            1.01 * 1.02 * 0.1 * 1.03 * 1.04 * 1.2,
                            1.02 * 1.1 * 1.03 * 1.04 * 1.2 * 1.3,
                            1.01 * 1.02 * 1.03 * 1.04 * 0.1 * 1.2 * 1.3)))
    ouo <- fouo(fp3)
    expect_true(all.equal(ouo[[1]], ouo[[2]]))
})

test_that("poset with all effects and modules, 1", {
    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                  child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                  s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                  sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                  typeDep = c(rep("--", 4), 
                      "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    fp4m <- allFitnessEffects(p4,
                              geneToModule = c("Root" = "Root", "A" = "a1",
                              "B" = "b1, b2", "C" = "c1",
                              "D" = "d1, d2", "E" = "e1",
                              "F" = "f1, f2", "G" = "g1"))
    gfp4 <- evalAllGenotypes(fp4m, order = FALSE, max = 1024)
    expect_true(all.equal(gfp4[c(12, 20, 21, 40, 41, 46,
                                 50, 55, 64, 92, 155, 157,
                                 163, 372, 632, 828), "Fitness"],
                          c(1.01 * 1.02, 1.02, 1.02 * 1.1, 0.1 * 1.3, 1.03, 
                            1.03 * 1.04, 1.04 * 0.05, 0.05 * 1.3,  
                            1.01 * 1.02 * 0.1, 1.02 * 1.1, 0.1 * 0.05 * 1.3,
                            1.03 * 0.05, 1.03 * 0.05,
                            1.03 * 1.04 * 1.2, 1.03 * 1.04 * 1.2, 
                            1.02 * 1.1 * 1.03 * 1.04 * 1.2 * 1.3)))
})


test_that("breaks if not all modules", {
    p4 <- data.frame(parent = c(rep("Root", 4), "A", "B", "D", "E", "C", "F"),
                     child = c("A", "B", "D", "E", "C", "C", "F", "F", "G", "G"),
                     s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3),
                     sh = c(rep(0, 4), c(-.9, -.9), c(-.95, -.95), c(-.99, -.99)),
                     typeDep = c(rep("--", 4), 
                         "XMPN", "XMPN", "MN", "MN", "SM", "SM"))
    expect_error(allFitnessEffects(p4,
                                   geneToModule = c("Root" = "Root", "A" = "a1",
                                       "B" = "b1, b2", "C" = "c1",
                                       "D" = "d1, d2", "E" = "e1",
                                       "F" = "f1, f2")))
})


test_that("Bauer example: exercising drvNames", {
    sd <- 0.1
    sdp <- 0.15
    sp <- 0.05
    bauer <- data.frame(parent = c("Root", rep("p", 5)),
                        child = c("p", paste0("s", 1:5)),
                        s = c(sd, rep(sdp, 5)),
                        sh = c(0, rep(sp, 5)),
                        typeDep = "MN")
    b1 <- evalAllGenotypes(allFitnessEffects(bauer, drvNames = c("s1", "s5")),
                           order = FALSE)
    b2 <- evalAllGenotypes(allFitnessEffects(bauer,
                                             drvNames = c("s2", "s3", "s4")),
                           order = TRUE, max = 2000)
    expect_equal(length(unique(b1$Fitness)), 11)
    expect_equal(length(unique(b2$Fitness)), 11)
} )


test_that("non distinct gene names per module caught", {
    expect_error(ofe1 <-
                     allFitnessEffects(
                         orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
                         geneToModule =
                             c("Root" = "Root",
                               "F" = "f1, f2",
                               "D" = "d1, d2, f1") ),
                 "Are there identical gene names in different modules?")
    s <- 0.2
    expect_error(
        m1 <- allFitnessEffects(data.frame(
            parent = c("Root", "A"),
            child = c("A", "B"),
            s = s,
            sh = -1,
            typeDep = "OR"),
            geneToModule = c("Root" = "Root",
                             "A" = "a1, b1, a2",
                             "B" = "b1")),
        "Are there identical gene names in different modules?")
})


test_that("gene by itself and in module", {
    expect_error(ofe1 <- allFitnessEffects(data.frame(parent = c("Root", "Root", "A"),
                                                      child = c("A", "e", "B"),
                                 s = .1,
                                 sh = .1,
                                 typeDep = "OR"),
                      geneToModule =
                          c("Root" = "Root",
                            "e" = "e",
                            "A" = "a1, e, a2",
                            "B" = "b1"))
                ,
                 "Are there identical gene names in different modules?")
    ## I think the "Is a gene part of two ... cannot be reached anymore"
})



test_that("a silly epistasis example", {
    ## make sure we exercise the nrow(df) == 0L
    expect_silent(sv <- allFitnessEffects(
                     orderEffects = c("A > B" = 0.1),
                     epistasis = c("A" = 1)))
    expect_output(print(evalAllGenotypes(sv, order = TRUE)), "Genotype")
})

test_that("can run without keeping input", {
    cs <-  data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                      child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = 0.1,
                 sh = -0.9,
                 typeDep = "MN")
    expect_silent(cbn1 <- allFitnessEffects(cs, keepInput = FALSE))
})



test_that("we are exercising evalGenotype with a comma, echo, and proNeg", {
    ## we do this in the vignette already
    ofe2 <- allFitnessEffects(orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
                              geneToModule =
                                  c("Root" = "Root",
                                    "F" = "f1, f2, f3",
                                    "D" = "d1, d2") )
    expect_equal(evalGenotype("d1 , d2, f3", ofe2, verbose = TRUE, echo = TRUE),
                 1.4)
    expect_equal(evalGenotype("f3 , d1 , d2", ofe2, verbose = TRUE, echo = TRUE),
                 0.7)
    expect_equal(evalGenotype("f3 , d1 , d2", ofe2, verbose = TRUE,
                              echo = TRUE, model = "Bozic"),
                 1.3)
})


test_that("We limit number of genotypes in eval", {
    expect_error(evalAllGenotypes(
        allFitnessEffects(
            noIntGenes = runif(10)), order = FALSE),
        "There are 1024 genotypes. This is larger than max.",
        fixed = TRUE)
})

## how is table geneModule with no ints? are they there?

## check it breaks if same ID

## check breaks if in restriction and no interactions

## Nope, this is not inconsistent

test_that("Bozic limit cases handled consistently", {
    sv <- allFitnessEffects(data.frame(
        parent = c("Root", "Root", "a1", "a2"),
        child = c("a1", "a2", "b", "b"),
        s = 1.2,
        sh = 0.1,
        typeDep = "OR"),
        noIntGenes = c("E" = 0.85, "F" = 1))
    expect_output(print(evalAllGenotypes(sv, order = FALSE, addwt = TRUE,
                                   model = "Bozic")), ## this works
                  "Death_rate", fixed = TRUE, all = FALSE)
    expect_error(oncoSimulIndiv(sv, model = "Bozic"),
                 "You are using a Bozic model with the new restriction specification, and you have at least one s > 1."
                 ) 
    sv2 <- allFitnessEffects(epistasis = c("-A : B" = 1.5,
                                           "A : -B" = 1.5,
                                           "A:B" = 2.5))
    expect_output(print(evalAllGenotypes(sv2, order = FALSE, addwt = TRUE,
                                   model = "Bozic")), ## this works
                  "Death_rate", fixed = TRUE, all = FALSE)
    expect_error(oncoSimulIndiv(sv2, model = "Bozic"),
                 "You are using a Bozic model with the new restriction specification, and you have at least one s > 1."
                 ) 
    sv3 <- allFitnessEffects(epistasis = c("-A : B" = 0.1,
                                           "A : -B" = 1,
                                           "A:B" = 0.1))
    expect_output(print(evalAllGenotypes(sv3, order = FALSE, addwt = TRUE,
                                   model = "Bozic")), ## this works
                  "Death_rate", fixed = TRUE, all = FALSE)
    expect_warning(oncoSimulIndiv(sv3, model = "Bozic"),
                   "You are using a Bozic model with the new restriction specification, and you have at least one s of 1."
                   ) 
    sv4 <- allFitnessEffects(epistasis = c("-A : B" = 0.99,
                                           "A : -B" = 0.95,
                                           "A:B" = 0.5),
                             orderEffects = c("G > H" = 1.1),
                             noIntGenes = c("E" = 0.85, "F" = 1.35))
    expect_output(print(evalAllGenotypes(sv4, order = FALSE, addwt = TRUE,
                                   model = "Bozic")), ## this works
                  "Death_rate", fixed = TRUE, all = FALSE)
    expect_error(oncoSimulIndiv(sv4, model = "Bozic"),
                 "You are using a Bozic model with the new restriction specification, and you have at least one s > 1."
                 ) 
    sv5 <- allFitnessEffects(epistasis = c("-A : B" = 0.99,
                                           "A : -B" = 0.95,
                                           "A:B" = 1),
                             orderEffects = c("G > H" = 1),
                             noIntGenes = c("E" = 0.85, "F" = 1))
    expect_output(print(evalAllGenotypes(sv5, order = FALSE, addwt = TRUE,
                                   model = "Bozic")), ## this works
                  "Death_rate", fixed = TRUE, all = FALSE)
    expect_warning(oncoSimulIndiv(sv5, model = "Bozic"),
                   "You are using a Bozic model with the new restriction specification, and you have at least one s of 1."
                   ) 
    sv4c <- allFitnessEffects(epistasis = c("-A : B" = 0.99,
                                           "A : -B" = 1.5,
                                           "A:B" = 0.5),
                             orderEffects = c("G > H" = 1),
                             noIntGenes = c("E" = 0.85, "F" = 1.35))
    expect_output(print(evalAllGenotypes(sv4c, order = FALSE, addwt = TRUE,
                                   model = "Bozic")), ## this works
                  "Death_rate", fixed = TRUE, all = FALSE)
    expect_error(oncoSimulIndiv(sv4c, model = "Bozic"),
                 "You are using a Bozic model with the new restriction specification, and you have at least one s > 1."
                 ) 
    sv4d <- allFitnessEffects(epistasis = c("-A : B" = 0.99,
                                           "A : -B" = 0.95,
                                           "A:B" = 0.5),
                             orderEffects = c("G > H" = 1.1),
                             noIntGenes = c("E" = 0.85, "F" = 0.35))
    expect_output(print(evalAllGenotypes(sv4d, order = FALSE, addwt = TRUE,
                                   model = "Bozic")), ## this works
                  "Death_rate", fixed = TRUE, all = FALSE)
    expect_error(oncoSimulIndiv(sv4d, model = "Bozic"),
                 "You are using a Bozic model with the new restriction specification, and you have at least one s > 1."
                 ) 
    sv4d <- allFitnessEffects(epistasis = c("-A : B" = 0.99,
                                           "A : -B" = 0.95,
                                           "A:B" = 0.5),
                             orderEffects = c("G > H" = 0.1),
                             noIntGenes = c("E" = 0.85, "F" = 1.3))
    expect_output(print(evalAllGenotypes(sv4d, order = FALSE, addwt = TRUE,
                                   model = "Bozic")), ## this works
                  "Death_rate", fixed = TRUE, all = FALSE)
    expect_error(oncoSimulIndiv(sv4d, model = "Bozic"),
                 "You are using a Bozic model with the new restriction specification, and you have at least one s > 1."
                 )
    svff <- allFitnessEffects(data.frame(
        parent = c("Root", "Root", "a1", "a2"),
        child = c("a1", "a2", "b", "b"),
        s = 0.2,
        sh = -1,
        typeDep = "OR"),
        noIntGenes = c("E" = 0.85, "F" = .1))
    expect_output(print(evalAllGenotypes(svff, order = FALSE, addwt = TRUE,
                                   model = "Bozic")), ## this works
                  "Death_rate", fixed = TRUE, all = FALSE)
    expect_warning(oncoSimulIndiv(svff, model = "Bozic"),
                 "You are using a Bozic model with the new restriction specification, and you have at least one sh <= -1."
                 ) 
    svff2 <- allFitnessEffects(data.frame(
        parent = c("Root", "Root", "a1", "a2"),
        child = c("a1", "a2", "b", "b"),
        s = 0.2,
        sh = -1.5,
        typeDep = "OR"),
        noIntGenes = c("E" = 0.85, "F" = .1))
    expect_output(print(evalAllGenotypes(svff2, order = FALSE, addwt = TRUE,
                                   model = "Bozic")), ## this works
                  "Death_rate", fixed = TRUE, all = FALSE)
    expect_warning(oncoSimulIndiv(svff2, model = "Bozic"),
                 "You are using a Bozic model with the new restriction specification, and you have at least one sh <= -1."
                 ) 
    svff3 <- allFitnessEffects(data.frame(
        parent = c("Root", "Root", "a1", "a2"),
        child = c("a1", "a2", "b", "b"),
        s = 0.2,
        sh = -Inf,
        typeDep = "OR"),
        noIntGenes = c("E" = 0.85, "F" = .1))
    expect_output(print(evalAllGenotypes(svff3, order = FALSE, addwt = TRUE,
                                   model = "Bozic")), ## this works
                  "Death_rate", fixed = TRUE, all = FALSE)
    expect_output(print(oncoSimulIndiv(svff3, model = "Bozic",
                                       sampleEvery = 0.02)),
                 "Individual OncoSimul trajectory with call"
                 ) 
})




test_that("No epistasis, modules", {
    sa <- 0.01
    sb <- 0.02
    sc <- 0.03
    fnme <- allFitnessEffects(epistasis = c("A" = sa,
                                            "B" = sb,
                                            "C" = sc),
                              geneToModule = c("A" = "a1, a2",
                                               "B" = "b1",
                                               "C" = "c1, c2"))
    ea <- evalAllGenotypes(fnme, order = FALSE, addwt = TRUE)
    expect_identical(ea[ea$Genotype == "a1, a2", "Fitness"], 1 + sa)
    expect_identical(ea[ea$Genotype == "a1, b1", "Fitness"],
    (1 + sa) * (1 + sb))
    expect_identical(ea[ea$Genotype == "a2, c1", "Fitness"],
    (1 + sa) * (1 + sc))
    expect_identical(ea[ea$Genotype == "b1, c2", "Fitness"],
    (1 + sb) * (1 + sc))
    expect_identical(ea[ea$Genotype == "a1, a2, c1", "Fitness"],
    (1 + sa) * (1 + sc))
    expect_identical(ea[ea$Genotype == "a1, a2, b1, c1", "Fitness"],
    (1 + sa) * (1 + sb) * (1 + sc))
    expect_identical(ea[ea$Genotype == "a2, b1, c1", "Fitness"],
    (1 + sa) * (1 + sb) * (1 + sc))
})


test_that("noIntGenes, two common errors: character vector and ,>:", {
    expect_error(allFitnessEffects(noIntGenes =
                                       c("a1, a2, b1, b2, b3, c1, c2")),
                 "noIntGenes is a character vector", fixed = TRUE)

    expect_error(allFitnessEffects(noIntGenes =
                                       c("a1", "a2")),
                 "noIntGenes is a character vector", fixed = TRUE)
    u <- 0.3
    names(u) <- "a,b"
    v <- 0.3
    names(v) <- "a > b"
    w <- 0.3
    names(w) <- "a : b"

    expect_error(allFitnessEffects(noIntGenes = u),
                 "The name of some noIntGenes contain a ',' or a '>' or a ':'",
                 fixed = TRUE)
    expect_error(allFitnessEffects(noIntGenes = v),
                 "The name of some noIntGenes contain a ',' or a '>' or a ':'",
                 fixed = TRUE)
    expect_error(allFitnessEffects(noIntGenes = w),
                 "The name of some noIntGenes contain a ',' or a '>' or a ':'",
                 fixed = TRUE)
})


test_that("Some same genes in epistasis and order effects", {
    s1 <- 0.1
    s2 <- 0.2
    s3 <- 0.3
    s4 <- 0.9
    s5 <- 0.7
    s0 <- 0.33
    sh <- -Inf
    o999 <- allFitnessEffects(rT = data.frame(parent = c("Root", "a", "f"),
                                              child  = c("a", "f", "m"),
                                              s = s0,
                                              sh = sh,
                                              typeDep = "MN"),
                              orderEffects = c("a>b" = s1, "b > a" = s2, "b > m" = s3),
                              epistasis = c("a:c" = s4, "b:e" = s5))
    of <- evalAllGenotypes(o999, order = TRUE, max = 1956)
    expect_equal(dplyr::filter(of, Genotype == "b > a > f > c")[, "Fitness"],
    (1 + s2) * (1 + s0) * (1 + s0) * (1 + s4))
    expect_equal(dplyr::filter(of, Genotype == "a > f > c > b > m")[, "Fitness"],
    (1 + s0) * (1 + s0) * (1 + s0) * (1 + s1) * (1 + s4) * (1 + s3))
    expect_equal(dplyr::filter(of, Genotype == "e > a > b")[, "Fitness"],
    (1 + s0) * (1 + s5) * (1 + s1))
    s1 <- -0.2
    s2 <- 0.3
    s3 <- 0.9
    o99 <- allFitnessEffects(
        orderEffects = c("a>b" = s1, "b > a" = s2),
        epistasis = c("a:c" = s3))
    eo99 <- evalAllGenotypes(o99, order = TRUE, addwt = TRUE)
    expect_equal(dplyr::filter(eo99, Genotype == "a > c > b")[, "Fitness"],
                 (1 + s3) * (1 + s1))
    expect_equal(dplyr::filter(eo99, Genotype == "c > a > b")[, "Fitness"],
                 (1 + s3) * (1 + s1))
    expect_equal(dplyr::filter(eo99, Genotype == "c > b > a")[, "Fitness"],
                 (1 + s3) * (1 + s2))
    expect_equal(dplyr::filter(eo99, Genotype == "c > a")[, "Fitness"],
                 (1 + s3))
    expect_equal(dplyr::filter(eo99, Genotype == "a > c")[, "Fitness"],
                 (1 + s3))
})


test_that("noIntGenes errors", {
    expect_error(
        pancrr <- allFitnessEffects(
        data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                              "TP53", "TP53", "MLL3"),
                   child = c("KRAS","SMAD4", "CDNK2A", 
                             "TP53", "MLL3",
                             rep("PXDN", 3), rep("TGFBR2", 2)),
                   s = .2,
                   sh = .3,
                   typeDep = "MN"),
        noIntGenes = c("TP53" = 0.1)),
        "A gene in noIntGenes also present in the other terms",
        fixed = TRUE)
    expect_error(
        nr <- allFitnessEffects(
            noIntGenes = c("A" = 0.1, "B" = 0.2, "A" = 0.05)),
        "Duplicated gene names in geneNoInt",
        fixed = TRUE)
})


test_that("not all genes named", {
    gg <- rep(0.01, 3)
    names(gg) <- letters[1:2]
    expect_error(allFitnessEffects(noIntGenes = gg),
                 "In noIntGenes some genes have names, some don't.",
                 fixed = TRUE)
})


test_that("We can deal with single-gene genotypes and trivial cases" ,{

    ## we get the message
    expect_message(allFitnessEffects(
        genotFitness = data.frame(g = c("A", "B"),
                                  y = c(1, 2))), "All single-gene genotypes",
        fixed = TRUE)

    
    expect_true(identical(
        data.frame(Genotype = c("WT", "A", "B", "A, B"),
                   Fitness = c(1.0, 1.0, 2.0, 0.0), ## 0.0 used to be 1.0
                   stringsAsFactors = FALSE),
        as.data.frame(evalAllGenotypes(
            allFitnessEffects(genotFitness = data.frame(g = c("A", "B"),
                                                        y = c(1, 2))),
            addwt = TRUE))
    ))

    expect_true(identical(
        data.frame(Genotype = c("WT", "A", "B", "A, B"),
                   Fitness = c(1.0, 1.5, 2.9, 0.0),
                   stringsAsFactors = FALSE),
        as.data.frame(evalAllGenotypes(
            allFitnessEffects(genotFitness = data.frame(g = c("A", "B"),
                                                        y = c(1.5, 2.9))),
            addwt = TRUE))
    ))

    expect_true(identical(
        data.frame(Genotype = c("WT", "A", "B", "E", "A, B", "A, E", "B, E", "A, B, E"),
                   Fitness = c(1.0, 1.3, 2.4, 3.2, rep(0, 4)),
                   stringsAsFactors = FALSE),
        as.data.frame(evalAllGenotypes(
            allFitnessEffects(genotFitness = data.frame(g = c("A", "B", "E"),
                                                        y = c(1.3, 2.4, 3.2))),
            addwt = TRUE))
    ))


    expect_true(identical(
        data.frame(Genotype = c("WT", "A"),
                   Fitness = c(1.0, 1.0),
                   stringsAsFactors = FALSE),
        as.data.frame(evalAllGenotypes(
            allFitnessEffects(genotFitness = data.frame(g = c("A"),
                                                        y = c(1))),
            addwt = TRUE))
    ))
    
    expect_true(identical(
        data.frame(Genotype = c("WT", "A"),
                   Fitness = c(1.0, 0.6),
                   stringsAsFactors = FALSE),
        as.data.frame(evalAllGenotypes(
            allFitnessEffects(genotFitness = data.frame(g = c("A"),
                                                        y = c(0.6))),
            addwt = TRUE))
    ))

    expect_true(identical(
        data.frame(Genotype = c("WT", "A", "D", "F", "A, D", "A, F", "D, F",
                                "A, D, F"),
                   Fitness = c(1.0, rep(0, 6), 1.7), ## c(rep(1, 7), 1.7),
                   stringsAsFactors = FALSE),
        as.data.frame(evalAllGenotypes(
            allFitnessEffects(genotFitness = data.frame(g = c("A, D, F"),
                                                        y = c(1.7))),
            addwt = TRUE))
    ))    


    m <- rbind(c(1, 0, 1.2),
               c(0, 1, 2.4))

    expect_true(identical(
        data.frame(Genotype = c("WT", "A", "B", "A, B"),
                   Fitness = c(1.0, 1.2, 2.4, 0.0),
                   stringsAsFactors = FALSE),
        as.data.frame(evalAllGenotypes(
            allFitnessEffects(genotFitness = m),
            addwt = TRUE))
    ))    

    expect_message(evalAllGenotypes(
        allFitnessEffects(genotFitness = m)),
        "No column names", fixed = TRUE)

    mcn <- m
    colnames(mcn) <- c("A", "", "Fitness")
    expect_warning(evalAllGenotypes(
        allFitnessEffects(genotFitness = mcn)),
        "One column named ''", fixed = TRUE)
    rm(mcn)
    
    m2 <- rbind(c(1, 0, 1.2),
               c(0, 1, 2.4))
    colnames(m2) <- c("U", "M", "Fitness")
    expect_true(identical(
        data.frame(Genotype = c("WT", "M", "U", "M, U"),
                   Fitness = c(1.0, 2.4, 1.2, 0),
                   stringsAsFactors = FALSE),
        as.data.frame(evalAllGenotypes(
            allFitnessEffects(genotFitness = m2),
            addwt = TRUE))
    ))

    expect_message(
        evalAllGenotypes(
            allFitnessEffects(genotFitness = m2)),
        "Sorting gene column names", fixed = TRUE
    )

    
    m2df <- data.frame(rbind(c(1, 0, 1.2),
               c(0, 1, 2.4)))
    colnames(m2df) <- c("U", "M", "Fitness")
    expect_true(identical(
        data.frame(Genotype = c("WT", "M", "U", "M, U"),
                   Fitness = c(1.0, 2.4, 1.2, 0),
                   stringsAsFactors = FALSE),
        as.data.frame(evalAllGenotypes(
            allFitnessEffects(genotFitness = m2df),
            addwt = TRUE))
    ))   

    m3 <- matrix(c(1, 1.2), ncol = 2)
    colnames(m3) <- c("U", "Fitness")
    expect_error(
        allFitnessEffects(genotFitness = m3),
        "genotFitness: if two-column must be data frame",
        fixed = TRUE)

    ## Stupid
    m5 <- data.frame(x = 1, y = 2, stringsAsFactors= FALSE)
    expect_error(evalAllGenotypes(allFitnessEffects(genotFitness = m5)),
                 "genotFitness: first column of data frame is numeric.",
                 fixed = TRUE)

    m6 <- matrix(letters[1:4], ncol = 4)
    expect_error(allFitnessEffects(genotFitness = m6),
                 "A genotype fitness matrix/data.frame must be numeric",
                 fixed = TRUE)

    m7 <- as.data.frame(matrix(letters[1:4], ncol = 4))
    expect_error(allFitnessEffects(genotFitness = m7),
                 "A genotype fitness matrix/data.frame must be numeric",
                 fixed = TRUE)

    m8 <- 1:9
    expect_error(allFitnessEffects(genotFitness = m8),
                 "Input must inherit from matrix or data.frame",
                 fixed = TRUE)
})


cat(paste("\n Ending all-fitness at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
