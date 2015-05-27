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

test_that("Bauer example: correct number of fitness classes", {
    expect_equal(length(unique(b1$Fitness)), 11)
    expect_equal(length(unique(b2$Fitness)), 11)
} )

test_that("Bauer example: identical values fitness classes, unorder and ord", {
    expect_equal(unique(b1$Fitness), unique(b2$Fitness))
} )


test_that("Bauer example: identical values fitness classes, rename", {
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
    bauer3 <- data.frame(parent = c(rep("u", 5), "Root"),
                         child = c(paste0("s", 1:5), "u"),
                         s = c(sd, rep(sdp, 5)),
                         sh = c(0, rep(sp, 5)),
                         typeDep = "MN",
                         stringsAsFactors = FALSE)
    b3 <- evalAllGenotypes(allFitnessEffects(bauer), order = TRUE, max = 2000)
    expect_equal(unique(b1$Fitness), unique(b3$Fitness))
} )



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



foi1 <- allFitnessEffects(
    orderEffects = c("D>B" = -0.2, "B > D" = 0.3),
    noIntGenes = c("A" = 0.05, "C" = -.2, "E" = .1))
agoi1 <- evalAllGenotypes(foi1,  max = 325)

rownames(agoi1) <- agoi1[, 1]
agoi1[LETTERS[1:5], ]
all(agoi1[LETTERS[1:5], "Fitness"] == c(1.05, 1, 0.8, 1, 1.1))
## orders that do not involve all. D > A;   B > C;

rn <- rownames(agoi1)
agoi1[grep("^A > [BD]$", rn), ]







af1 <- evalAllGenotypes(f1, max = 325)

agoi1 <- evalAllGenotypes(allFitnessEffects(
    orderEffects = c("A>B" = -0.2, "B > A" = 0.3),
    noIntGenes = c(0.05, -.2, .1)), max = 325)


gm1 <- c("Root" = "Root", "F" = "f1, f2", "D" = "d1, d2, d3")
gm.twisted1 <- c("Root" = "Root", "F" = "d", "D" = "f")
gm.twisted2 <- c("Root" = "Root", "F" = "d1, d2", "D" = "f1, f2, f3")


    ofe2 <- allFitnessEffects(orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
                              geneToModule =
                                  c("Root" = "Root",
                                    "F" = "f1, f2, f3",
                                    "D" = "d1, d2") )



## to the above, add genes with no interactions and see they have no effect

## add three gene order restrictions


## add differences in module names or something similar


## how is table geneModule with no ints? are they there?




## check it breaks if same ID















