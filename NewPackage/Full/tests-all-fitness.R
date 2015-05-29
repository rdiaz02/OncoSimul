

test_that("Bauer example: correct number of fitness classes", {
    sd <- 0.1
    sdp <- 0.15
    sp <- 0.05
    bauer <- data.frame(parent = c("Root", rep("p", 5)),
                        child = c("p", paste0("s", 1:5)),
                        s = c(sd, rep(sdp, 5)),
                        sh = c(0, rep(sp, 5)),
                        typeDep = "MN",
                        stringsAsFactors = FALSE)
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
                        typeDep = "MN",
                        stringsAsFactors = FALSE)
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
                        typeDep = "MN",
                        stringsAsFactors = FALSE)
    b1 <- evalAllGenotypes(allFitnessEffects(bauer), order = FALSE)
    bauer3 <- data.frame(parent = c("Root", rep("u", 5)),
                         child = c("u", paste0("s", 1:5)),
                         s = c(sd, rep(sdp, 5)),
                         sh = c(0, rep(sp, 5)),
                         typeDep = "MN",
                         stringsAsFactors = FALSE)
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
                        typeDep = "MN",
                        stringsAsFactors = FALSE)
    b1 <- evalAllGenotypes(allFitnessEffects(bauer), order = FALSE)
    bauer3 <- data.frame(parent = c(rep("u", 5), "Root"),
                         child = c(paste0("s", 1:5), "u"),
                         s = c(sd, rep(sdp, 5)),
                         sh = c(0, rep(sp, 5)),
                         typeDep = "MN",
                         stringsAsFactors = FALSE)
    b3 <- evalAllGenotypes(allFitnessEffects(bauer), order = TRUE, max = 2000)
    expect_equal(unique(b1$Fitness), unique(b3$Fitness))
} )


### Order effects


test_that("Order effects, entry of labels and separation", {
    o1 <- evalAllGenotypes(allFitnessEffects(
        orderEffects = c("d>f" = 0.4, "f > d" = -0.3) ))
    o2 <- evalAllGenotypes(allFitnessEffects(
        orderEffects = c("f > d" = -0.3, "d > f" = 0.4) ))
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
    ag <- evalAllGenotypes(ofe1)
    expect_true(all(ag[c(17, 39, 19, 29), "Fitness"] == c(1.4, 0.7, 1.4, 0.7)))
    expect_true(all(ag[c(43, 44), "Fitness"] == c(1.4, 1.4)))
    expect_true(all(ag[41:52, "Fitness"] == 1.4))
    expect_true(all(ag[53:64, "Fitness"] == 0.7))
})

test_that("Order effects, modules 2", {
    ofe2 <- allFitnessEffects(orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
                              geneToModule =
                                  c("Root" = "Root",
                                    "F" = "f1, f2, f3",
                                    "D" = "d1, d2") )
    ag2 <- evalAllGenotypes(ofe2, max = 326)
    oe <- c(grep("^f.*d.*", ag2[, 1]), grep("^d.*f.*", ag2[, 1]))
    expect_true(all(ag2[grep("^d.*f.*", ag2[, 1]), "Fitness"] == 1.4))
    expect_true(all(ag2[grep("^f.*d.*", ag2[, 1]), "Fitness"] == 0.7))
    expect_true(all(ag2[-oe, "Fitness"] == 1))
})



test_that("Order effects, twisted module names", {
    o1 <- evalAllGenotypes(allFitnessEffects(
      orderEffects = c("F > D" = -0.3, "D > F" = 0.4),  
        geneToModule = c("Root" = "Root", "F" = "d", "D" = "f")))
    o2 <- evalAllGenotypes(allFitnessEffects(
      orderEffects = c("D>F" = 0.4, "F >D" = -0.3),  
        geneToModule = c("Root" = "Root", "F" = "d", "D" = "f")))
    o1s <- o1[order(o1$Genotype),   ]
    o2s <- o2[order(o2$Genotype),   ]
    expect_equal(o1s, o2s)
    expect_true(all(o1[, "Fitness"] == c(1, 1, 0.7, 1.4)))
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
    ag <- evalAllGenotypes(o3)
    expect_true(all(ag[, "Fitness"] ==
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

    
    expect_true(all(ai1[, "Fitness"]  == c( (1 + .05), (1 - .2), (1 + .1),
       (1 + .05) * (1 - .2),
       (1 + .05) * (1 + .1),
       (1 - .2) * (1 + .1),
       (1 + .05) * (1 - .2) * (1 + .1))))

    expect_true(all(ai2[, "Fitness"] == c((1 + .05), (1 - .2), (1 + .1),
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
    
    expect_true(all(ai3[, "Fitness"]  == c( (1 + .05), (1 - .2), (1 + .1),
       (1 + .05) * (1 - .2),
       (1 + .05) * (1 + .1),
       (1 - .2) * (1 + .1),
       (1 + .05) * (1 - .2) * (1 + .1))))

    expect_true(all(ai4[, "Fitness"] == c((1 + .05), (1 - .2), (1 + .1),
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
    
    expect_true(all(ai3[, "Fitness"]  == c( (1 + .05), (1 - .2), (1 + .1),
       (1 + .05) * (1 - .2),
       (1 + .05) * (1 + .1),
       (1 - .2) * (1 + .1),
       (1 + .05) * (1 - .2) * (1 + .1))))

    expect_true(all(ai4[, "Fitness"] == c((1 + .05), (1 - .2), (1 + .1),
                           1.05 * .8, 1.05 * 1.1, .8 * 1.05, .8 * 1.1,
                           1.05 * 1.1, 1.1 * .8,
                           rep(1.05 * .8 * 1.1, 6) )))

})




test_that("No interaction genes and order effects, 1", {
    foi1 <- allFitnessEffects(
        orderEffects = c("D>B" = -0.3, "B > D" = 0.3),
        noIntGenes = c("A" = 0.05, "C" = -.2, "E" = .1))
    agoi1 <- evalAllGenotypes(foi1,  max = 325)
    rn <- 1:nrow(agoi1)
    names(rn) <- agoi1[, 1]
    expect_true(all(agoi1[rn[LETTERS[1:5]], "Fitness"] == c(1.05, 1, 0.8, 1, 1.1)))
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
    agoi1 <- evalAllGenotypes(foi1,  max = 325)
    rn <- 1:nrow(agoi1)
    names(rn) <- agoi1[, 1]
    expect_true(all(agoi1[rn[c("B", "D", "M", "A", "J")], "Fitness"] == c(1, 1, 1.05, 0.8, 1.1)))
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
    expect_true( all(
        evalAllGenotypes(sv,
                         order = FALSE, addwt = TRUE)[, "Fitness"] == c(1, 0, 0, 1.2)))
    expect_true( all(
        evalAllGenotypes(sv,
                         order = TRUE, addwt = TRUE)[, "Fitness"] == c(1, 0, 0, 1.2, 1.2)))
})


test_that("synthetic viability, 2", {
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
    expect_true( all(
        evalAllGenotypes(sv2,
                         order = FALSE, addwt = TRUE)[, "Fitness"] ==
                             c(1, 0.9, 0.9, 0.8, 0.9, rep(1.25, 3))))
    expect_true( all(
        evalAllGenotypes(sv2,
                         order = TRUE, addwt = TRUE)[, "Fitness"] ==
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
    expect_true( all(
        evalAllGenotypes(sm1, order = FALSE, addwt = TRUE)[, "Fitness"] ==
            c(1, 1.1, 1.2, 1 - 0.8)))
    expect_true( all(
        evalAllGenotypes(sm1, order = TRUE, addwt = TRUE)[, "Fitness"] ==
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

    expect_true(all(evalAllGenotypes(e2,
                                     order = FALSE,
                                     addwt = TRUE)[, "Fitness"] ==
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
    expect_true(all(evalAllGenotypes(e2,
                                     order = FALSE,
                                     addwt = TRUE)[, "Fitness"] ==
                                         c(1, 1.2, 1.3, 1.7)))
    expect_true(all(evalAllGenotypes(e3,
                                     order = FALSE,
                                     addwt = TRUE)[, "Fitness"] ==
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

    expect_true(all(evalAllGenotypes(E3, order = FALSE, addwt = FALSE)[, "Fitness"] ==
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
    expect_true(all(evalAllGenotypes(E3A, order = FALSE, addwt = FALSE)[, "Fitness"] ==
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
    expect_true(all(evalAllGenotypes(E3A, order = FALSE, addwt = FALSE) == 
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




## this is an error as CBN
c1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                 sh = c(rep(0, 4), c(-.1, -.2), c(-.05, -.06, -.07)),
                 typeDep = "MN",
                 stringsAsFactors = FALSE)

fc1 <- allFitnessEffects(c1)


## now OK
c1 <- data.frame(parent = c(rep("Root", 4), "a", "b", "d", "e", "c"),
                 child = c("a", "b", "d", "e", "c", "c", rep("g", 3)),
                 s = c(0.01, 0.02, 0.03, 0.04, 0.1, 0.1, rep(0.2, 3)),
                 sh = c(rep(0, 4), c(-.1, -.1), rep(-.05, 3)),
                 typeDep = "MN",
                 stringsAsFactors = FALSE)

fc1 <- allFitnessEffects(c1)



## how is table geneModule with no ints? are they there?




## check it breaks if same ID

## check breaks if in restriction and no interactions













