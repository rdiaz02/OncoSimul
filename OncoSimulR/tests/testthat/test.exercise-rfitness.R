## This does test some options, but mainly exercises others
inittime <- Sys.time()
cat(paste("\n Starting exercise-rfitness at", date(), "\n"))
test_that("Expect output", {

    expect_output(print(rfitness(4)), "Birth", fixed = TRUE)
    expect_output(print(rfitness(4, reference = "max")), "Birth",
                  fixed = TRUE)
    expect_output(print(rfitness(3, reference = c(1, 0, 1))), "Birth",
                  fixed = TRUE)
    expect_output(print(rfitness(5, scale = c(1, 7))), "Birth", fixed = TRUE)
    expect_output(print(rfitness(5, scale = c(1, 4), wt_is_1 = "no",
                           log = TRUE)), "Birth", fixed = TRUE)
    expect_output(print(rfitness(4, reference = "random2")), "Birth",
                  fixed = TRUE)
    expect_output(print(rfitness(4, reference = "random2",
                                 min_accessible_genotypes = 6)), "Birth",
                  fixed = TRUE)
    expect_output(print(rfitness(3, reference = "random2",
                                 min_accessible_genotypes = 6)), "Birth",
                  fixed = TRUE)
    expect_output(print(rfitness(3, reference = "max",
                                 min_accessible_genotypes = 6)), "Birth",
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


    expect_output(print(rfitness(4, model = "Ising")), "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Eggbox")), "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Full")), "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Additive")), "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Ising", i = 0.0002, I = 0.5,
                             circular = TRUE)), "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Ising", i = 2, I = 0,
                             circular = TRUE)), "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Ising", i = 2, I = -2,
                             circular = TRUE)), "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Ising", i = 2, I = 0.5,
                             circular = FALSE)), "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Eggbox", e = 2, E = 0)),
              "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Eggbox", e = 2, E = 2.4)),
              "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Eggbox", e = 0, E = 2.4)),
              "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Full", i = 0.002, I = 2,
                K = 2, r = TRUE,
                p = 0.2, P = 0.3, o = 0.3, O = 1)), "Birth", fixed = TRUE)
expect_error(rfitness(4, model = "Full", K = 5), "It makes no sense",
              fixed = TRUE)
expect_output(print(rfitness(4, model = "Full", i = 0.002, I = 2,
                K = 2, r = TRUE,
                p = 0.2, P = 0.3, o = 0.3, O = 1,
                s = 0.5, S=0, d = 0.1,
                e = 1, E = 1.5,
                H = 0.5)), "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Additive", mu = 0, sd = 1)),
              "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Additive", mu = 1, sd = 0)),
              "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Additive", mu = 3, sd = 1)),
              "Birth", fixed = TRUE)
expect_output(print(rfitness(4, model = "Additive", mu = -4, sd = 10)),
              "Birth", fixed = TRUE)

})
 
test_that("Testing the three-element scale argument", {
    jj <- 0
    for(j in 1:100) {
    ## while(TRUE) {
        ## jj <- jj + 1; cat(paste("\n ######### Doing jj = ", jj, "\n"))
        ## First, WT not equal to 1
        ## Don't use too stupid arguments
        scalev <- sort(runif(3, -2, 5), decreasing = TRUE)
        scalev <- scalev[c(1, 3, 2)]
        
        wt_min <- wt_max <- FALSE

        i <- round(runif(1, 1, .Machine$integer.max - 1))
        ## cat(paste("    i = "), i, "\n")

        set.seed(i)
        ra <- rfitness(7, truncate_at_0 = FALSE, wt_is_1 = "no")
        
        if(max(ra[, "Birth"]) == ra[1, "Birth"]) {
            set.seed(i)
            expect_warning(rb <- rfitness(7, scale = scalev,
                                          truncate_at_0 = FALSE),
                           "WT has maximum fitness. Range will be from scale[2] to scale[3]",
                           fixed = TRUE)
            wt_max <- TRUE
            wt_min <- FALSE
        } else if(min(ra[, "Birth"]) == ra[1, "Birth"]) {
            set.seed(i)
            expect_warning(rb <- rfitness(7, scale = scalev,
                                          truncate_at_0 = FALSE),
                           "WT has minimum fitness. Range will be from scale[3] to scale[1]",
                           fixed = TRUE)
            wt_min <- TRUE
            wt_max <- FALSE
        } else {
            set.seed(i)
            rb <- rfitness(7, scale = scalev,
                       truncate_at_0 = FALSE)
        }
        
        rap <- ra[-1, "Birth"][ra[-1, "Birth"] > ra[1, "Birth"]]
        ran <- ra[-1, "Birth"][ra[-1, "Birth"] < ra[1, "Birth"]]
        rbp <- rb[-1, "Birth"][rb[-1, "Birth"] > scalev[3]]
        rbn <- rb[-1, "Birth"][rb[-1, "Birth"] < scalev[3]]

        ## Numerical issues with equality comparisons
        rap <- c(ra[1, "Birth"], rap)
        rbp <- c(scalev[3], rbp)
        ran <- c(ra[1, "Birth"], ran)
        rbn <- c(scalev[3], rbn)        

        expect_true(abs(length(rap) - length(rbp)) < 2)
        expect_true(abs(length(ran) - length(rbn)) < 2)

        if((length(rap) > 1) && (length(rap) == length(rbp)))
            expect_equivalent(cor(rap, rbp), 1)
        if((length(ran) > 1) && (length(ran) == length(rbn)))
            expect_equivalent(cor(ran, rbn), 1)

        expect_true( (length(ran) + length(rap)) == (2^7 + 1) )
        expect_true( (length(rbn) + length(rbp)) == (2^7 + 1) )
                    
        
        if(!wt_max) {
            expect_equivalent(max(rb[, "Birth"]), scalev[1])
        } else {
            expect_equivalent(max(rb[, "Birth"]), scalev[3])
        }
        if(!wt_min) {
            expect_equivalent(min(rb[, "Birth"]), scalev[2])
        } else {
            expect_equivalent(min(rb[, "Birth"]), scalev[3])            
        }

        expect_equivalent(rb[1, "Birth"], scalev[3])

        ## X, Y, Z: X: original, Y: after scaling
        ## Just dealing with the positive
        ## Do not assume WT to be 1
        ## Write as original * slope + intercept
        ## so you can compare against
        ## lm(rbp ~ rap)

        xM <- max(ra[, "Birth"])
        xW <- ra[1, "Birth"]
        if(!wt_max) {    
            y <- rap * (scalev[1] - scalev[3])/(xM - xW) +
                (scalev[3] - xW * ((scalev[1] - scalev[3])/(xM - xW)))
            expect_equivalent(rbp, y)
            rm(y)
        } else {
            expect_equivalent(scalev[3], rbp)
        }

        xm <- min(ra[, "Birth"])
        if(!wt_min) {
            y <- ran * (scalev[3] - scalev[2])/(xW - xm) +
                (scalev[3] - xW * ((scalev[3] - scalev[2])/(xW - xm)))
            expect_equivalent(rbn, y)
            rm(y)
        } else {
            expect_equivalent(scalev[3], rbn)
        }

        ## log transformation
        set.seed(i)
        suppressWarnings(
            rc <- rfitness(7, scale = exp(scalev), truncate_at_0 = FALSE,
                           log = TRUE))
        rcp <- rc[-1, "Birth"][rc[-1, "Birth"] > scalev[3]]
        rcn <- rc[-1, "Birth"][rc[-1, "Birth"] < scalev[3]]
        
        ## Numerical issues with equality comparisons
        rcp <- c(scalev[3], rcp)
        rcn <- c(scalev[3], rcn)
        
        expect_true(abs(length(rap) - length(rcp)) < 2)
        expect_true(abs(length(ran) - length(rcn)) < 2)
        
        if((length(rap) > 1) && (length(rap) == length(rcp)))
            expect_equivalent(cor(rap, exp(rcp)), 1)
        if((length(ran) > 1) && (length(ran) == length(rcn)))
            expect_equivalent(cor(ran, exp(rcn)), 1)
        
        ## X, U, Z: X: original, U: after scaling; Z: after log
        ## Just dealing with the positive
        ## Do not assume WT to be 1
        ## can compare against lm(exp(rcp) ~ rap)
        ## if(!isTRUE(all.equal(xM, xW, check.attributes = FALSE))) {
        if(!wt_max) {
            u <- rap * (exp(scalev[1]) - exp(scalev[3]))/(xM - xW) +
                (exp(scalev[3]) - xW * ((exp(scalev[1]) - exp(scalev[3]))/(xM - xW)))
            expect_equivalent(u, exp(rcp))
        } else {
            expect_equivalent(scalev[3], rcp)
        }

        ## if(!isTRUE(all.equal(xm, xW, check.attributes = FALSE))) {
        if(!wt_min) {
            v <- ran * (exp(scalev[3]) - exp(scalev[2]))/(xW - xm) +
                (exp(scalev[3]) - xW * ((exp(scalev[3]) - exp(scalev[2]))/(xW - xm)))
            expect_equivalent(exp(rcn), v)
        } else {
            expect_equivalent(scalev[3], rcn)
        }
        
        ## This shows, by the way, that the log-transformed ain't linear function of
        ## log of original
        if(!isTRUE(all.equal(xM, xW, check.attributes = FALSE))) {
            z <- log( (rap - xW) * ((exp(scalev[1]) - exp(scalev[3]))/(xM - xW)) + exp(scalev[3]))
            expect_equivalent(z, rcp)
        }
        
        ############################################
        ## Repeat setting WT = 1
        set.seed(i)
        rax <- rfitness(7, truncate_at_0 = FALSE, wt_is_1 = "subtract")
        rapx <- rax[-1, "Birth"][rax[-1, "Birth"] > 1]
        ranx <- rax[-1, "Birth"][rax[-1, "Birth"] < 1]
        rapx <- c(rax[1, "Birth"], rapx)
        ranx <- c(rax[1, "Birth"], ranx)
        
        expect_true(abs(length(rapx) - length(rbp)) < 2)
        expect_true(abs(length(ranx) - length(rbn)) < 2)
        
        if((length(rapx) > 1) && (length(rapx) == length(rbp)))
            expect_equivalent(cor(rapx, rbp), 1)
        if((length(ranx) > 1) && (length(ranx) == length(rbn)))
            expect_equivalent(cor(ranx, rbn), 1)

        if(length(rapx) < 2) expect_true(length(ranx) > 2)
        
        xMx <- max(rax[, "Birth"])
        xWx <- rax[1, "Birth"]
        
        if(!isTRUE(all.equal(xMx, xWx, check.attributes = FALSE))) {
            yx <- rapx * (scalev[1] - scalev[3])/(xMx - xWx) +
                (scalev[3] - xWx * ((scalev[1] - scalev[3])/(xMx - xWx)))
            expect_equivalent(rbp, yx)
            rm(yx)
        } else {
            expect_equivalent(scalev[3], rbp)
        }

        xmx <- min(rax[, "Birth"])
        if(!isTRUE(all.equal(xmx, xWx, check.attributes = FALSE))) {
            yx <- ran * (scalev[3] - scalev[2])/(xWx - xmx) +
                (scalev[3] - xW * ((scalev[3] - scalev[2])/(xWx - xmx)))
            expect_equivalent(rbn, yx)
            rm(yx)
        } else {
            expect_equivalent(scalev[3], rbn)
        }
        
        expect_true(abs(length(rapx) - length(rcp)) < 2)
        expect_true(abs(length(ranx) - length(rcn)) < 2)
        
        ## log
        if((length(rapx) > 1) && (length(rapx) == length(rcp)))
            expect_equivalent(cor(rapx, exp(rcp)), 1)
        if((length(ranx) > 1) && (length(ranx) == length(rcn)))
            expect_equivalent(cor(ranx, exp(rcn)), 1)
        if(!isTRUE(all.equal(xMx, xWx, check.attributes = FALSE))) {
            ux <- rapx * (exp(scalev[1]) - exp(scalev[3]))/(xMx - xWx) +
                (exp(scalev[3]) - xWx * ((exp(scalev[1]) - exp(scalev[3]))/(xMx - xWx)))
            expect_equivalent(ux, exp(rcp))
        } else {
            expect_equivalent(scalev[3], rcp)
        }

         if(!isTRUE(all.equal(xm, xW, check.attributes = FALSE))) {
            vx <- ran * (exp(scalev[3]) - exp(scalev[2]))/(xW - xm) +
                (exp(scalev[3]) - xW * ((exp(scalev[3]) - exp(scalev[2]))/(xW - xm)))
            expect_equivalent(exp(rcn), vx)
        } else {
            expect_equivalent(scalev[3], rcn)
        }
        
        ## This shows, by the way, that the log-transformed ain't linear function of
        ## log of original
        if(!isTRUE(all.equal(xMx, xWx, check.attributes = FALSE))) {
            zx <- log( (rapx - xWx) *
                       ((exp(scalev[1]) - exp(scalev[3]))/(xMx - xWx)) +
                       exp(scalev[3]))
            expect_equivalent(zx, rcp)
        } else {
            expect_equivalent(scalev[3], rcp)
        }
        
        suppressWarnings(
            rm(scalev, ra, rb, rap, ran, rbp, rbn, rc, rcp, rcn, u, v, y, z,
           rapx, ranx, rax, yx, ux, zx, vx, xmx, xWx, xMx, xW, xM, xm))
    }

})

cat(paste("\n Ending exercise-rfitness at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
