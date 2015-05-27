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


gm1 <- c("Root" = "Root", "F" = "f1, f2", "D" = "d1, d2, d3")
gm.twisted1 <- c("Root" = "Root", "F" = "d", "D" = "f")
gm.twisted2 <- c("Root" = "Root", "F" = "d1, d2", "D" = "f1, f2, f3")



ofe1 <- allFitnessEffects(orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2",
                                "D" = "d1, d2") )

ag <- evalAllGenotypes(ofe1)




ofe1 <- allFitnessEffects(orderEffects = c("F > D" = -0.3, "D > F" = 0.4),
                          geneToModule =
                              c("Root" = "Root",
                                "F" = "f1, f2",
                                "D" = "d1, d2") )
ag <- evalAllGenotypes(ofe1)
ag[c(17, 39, 19, 29), "Fitness"] == c(1.4, 0.7, 1.4, 0.7)
ag[c(43, 44), "Fitness"] == c(1.4, 1.4)
all(ag[41:52, "Fitness"] == 1.4)

## to the above, add genes with no interactions and see they have no effect

## add three gene order restrictions


## add differences in module names or something similar


## how is table geneModule with no ints? are they there?




## check it breaks if same ID















