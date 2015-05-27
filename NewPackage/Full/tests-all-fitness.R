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


## add differences in module names or something similar


## how is table geneModule with no ints? are they there?




















