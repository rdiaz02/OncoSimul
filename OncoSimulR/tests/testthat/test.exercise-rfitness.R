## This does test some options, but mainly exercises others
inittime <- Sys.time()
cat(paste("\n Starting exercise-rfitness at", date(), "\n"))
test_that("Expect output", {

    expect_output(print(rfitness(4)), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, reference = "max")), "Fitness",
                  fixed = TRUE)
    expect_output(print(rfitness(3, reference = c(1, 0, 1))), "Fitness",
                  fixed = TRUE)
    expect_output(print(rfitness(5, scale = c(1, 7))), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(5, scale = c(1, 4), wt_is_1 = "no",
                           log = TRUE)), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, reference = "random2")), "Fitness",
                  fixed = TRUE)
    expect_output(print(rfitness(4, reference = "random2",
                                 min_accessible_genotypes = 6)), "Fitness",
                  fixed = TRUE)
    expect_output(print(rfitness(3, reference = "random2",
                                 min_accessible_genotypes = 6)), "Fitness",
                  fixed = TRUE)
    expect_output(print(rfitness(3, reference = "max",
                                 min_accessible_genotypes = 6)), "Fitness",
                  fixed = TRUE)
})


test_that("Minimal tests of generate_matrix_genotypes", {
    ## By induction, if it works for the few first, should work for all,
    ## unless memory issues. And if we go beyond, say, 10 or 12, it can
    ## take long in slow machines.
    for(i in 1:13) {
        tmp <- OncoSimulR:::generate_matrix_genotypes(i)
        expect_true(nrow(tmp) == (2^i))
        expect_true(ncol(tmp) == i)
        cstmp <- colSums(tmp)
        lucstmp <- unique(cstmp)
        expect_true(length(lucstmp) == 1)
        expect_true(lucstmp[1] == ((2^i)/2)) ## yes, 2^(i - 1) but do full
        ## simple logic
        rm(tmp)
        rm(cstmp)
        rm(lucstmp)
    }
})

test_that("Warnings if scale out of scale", {
    expect_warning(rfitness(4, wt_is_1 = "force", scale = c(0, 0.5)),
                   "Using wt_is_1 = force", fixed = TRUE)
})

test_that("Error if wrong arguments", {
    expect_error(rfitness(NA),
                 "Number of genes argument (g) is null or NA",
                 fixed = TRUE)
    expect_error(rfitness(NULL),
                 "Number of genes argument (g) is null or NA",
                 fixed = TRUE)
})

test_that("Additive, Ising, Eggbox, Full exercising", {

    expect_output(print(rfitness(4, model = "Ising")), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Eggbox")), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Full")), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Additive")), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Ising", i = 0.0002, I = 0.5,
                                 circular = TRUE)), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Ising", i = 2, I = 0,
                                 circular = TRUE)), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Ising", i = 2, I = -2,
                                 circular = TRUE)), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Ising", i = 2, I = 0.5,
                                 circular = FALSE)), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Eggbox", e = 2, E = 0)),
                  "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Eggbox", e = 2, E = 2.4)),
                  "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Eggbox", e = 0, E = 2.4)),
                  "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Full", i = 0.002, I = 2,
                                 K = 2, r = TRUE,
                                 p = 0.2, P = 0.3, o = 0.3, O = 1)), "Fitness", fixed = TRUE)
    expect_error(rfitness(4, model = "Full", K = 5), "It makes no sense",
                 fixed = TRUE)
    expect_output(print(rfitness(4, model = "Full", i = 0.002, I = 2,
                                 K = 2, r = TRUE,
                                 p = 0.2, P = 0.3, o = 0.3, O = 1,
                                 s = 0.5, S=0, d = 0.1,
                                 e = 1, E = 1.5,
                                 H = 0.5)), "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Additive", mu = 0, sd = 1)),
                  "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Additive", mu = 1, sd = 0)),
                  "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Additive", mu = 3, sd = 1)),
                  "Fitness", fixed = TRUE)
    expect_output(print(rfitness(4, model = "Additive", mu = -4, sd = 10)),
                  "Fitness", fixed = TRUE)

})
 
test_that("Testing the three-element scale argument", {
    for(j in 1:5) {
        ## First, WT not equal to 1
        ## Don't use too stupid arguments
        scalev <- sort(runif(3, -2, 5), decreasing = TRUE)
        scalev <- scalev[c(1, 3, 2)]
        i <- round(100 * runif(1))
        set.seed(i)
        ra <- rfitness(7, truncate_at_0 = FALSE, wt_is_1 = "no")
        set.seed(i)
        rb <- rfitness(7, scale = scalev, truncate_at_0 = FALSE)

        rap <- ra[, "Fitness"][ra[, "Fitness"] >= ra[1, "Fitness"]]
        ran <- ra[, "Fitness"][ra[, "Fitness"] < ra[1, "Fitness"]]
        rbp <- rb[, "Fitness"][rb[, "Fitness"] >= scalev[3]]
        rbn <- rb[, "Fitness"][rb[, "Fitness"] < scalev[3]]

        expect_equivalent(cor(rap, rbp), 1)
        expect_equivalent(cor(ran, rbn), 1)

        expect_equivalent(max(rb[, "Fitness"]), scalev[1])
        expect_equivalent(min(rb[, "Fitness"]), scalev[2])
        expect_equivalent(rb[1, "Fitness"], scalev[3])

        ## X, Y, Z: X: original, Y: after scaling
        ## Just dealing with the positive
        ## Do not assume WT to be 1

        ## Write as original * slope + intercept
        ## so you can compare against
        ## lm(rbp ~ rap)
        xM <- max(ra[, "Fitness"])
        xW <- ra[1, "Fitness"]
        y <- rap * (scalev[1] - scalev[3])/(xM - xW) +
            (scalev[3] - xW * ((scalev[1] - scalev[3])/(xM - xW)))

        expect_equivalent(rbp, y)
        
        set.seed(i)
        rc <- rfitness(7, scale = exp(scalev), truncate_at_0 = FALSE, log = TRUE)
        
        rcp <- rc[, "Fitness"][rc[, "Fitness"] >= scalev[3]]
        rcn <- rc[, "Fitness"][rc[, "Fitness"] < scalev[3]]

        expect_equivalent(cor(rap, exp(rcp)), 1)
        expect_equivalent(cor(ran, exp(rcn)), 1)

        ## X, U, Z: X: original, U: after scaling; Z: after log
        ## Just dealing with the positive
        ## Do not assume WT to be 1
        ## can compare against lm(exp(rcp) ~ rap)

        u <- rap * (exp(scalev[1]) - exp(scalev[3]))/(xM - xW) +
            (exp(scalev[3]) - xW * ((exp(scalev[1]) - exp(scalev[3]))/(xM - xW)))
        expect_equivalent(u, exp(rcp))

        ## This shows, by the way, that the log-transformed ain't linear function of
        ## log of original
        z <- log( (rap - xW) * ((exp(scalev[1]) - exp(scalev[3]))/(xM - xW)) + exp(scalev[3]))
        expect_equivalent(z, rcp)

        ############################################
        ## Repeat setting WT = 1
        set.seed(i)
        rax <- rfitness(7, truncate_at_0 = FALSE, wt_is_1 = "subtract")
        rapx <- rax[, "Fitness"][rax[, "Fitness"] >= 1]
        ranx <- rax[, "Fitness"][rax[, "Fitness"] < 1]

        expect_equivalent(cor(rapx, rbp), 1)
        expect_equivalent(cor(ranx, rbn), 1)
        
        xMx <- max(rax[, "Fitness"])
        xWx <- rax[1, "Fitness"]
        yx <- rapx * (scalev[1] - scalev[3])/(xMx - xWx) +
            (scalev[3] - xWx * ((scalev[1] - scalev[3])/(xMx - xWx)))

        expect_equivalent(rbp, yx)

        ## log
        expect_equivalent(cor(rapx, exp(rcp)), 1)
        expect_equivalent(cor(ranx, exp(rcn)), 1)

        ux <- rapx * (exp(scalev[1]) - exp(scalev[3]))/(xMx - xWx) +
            (exp(scalev[3]) - xWx * ((exp(scalev[1]) - exp(scalev[3]))/(xMx - xWx)))
        expect_equivalent(ux, exp(rcp))

        ## This shows, by the way, that the log-transformed ain't linear function of
        ## log of original
        zx <- log( (rapx - xWx) * ((exp(scalev[1]) - exp(scalev[3]))/(xMx - xWx)) + exp(scalev[3]))
        expect_equivalent(zx, rcp)
    }
})

cat(paste("\n Ending exercise-rfitness at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
