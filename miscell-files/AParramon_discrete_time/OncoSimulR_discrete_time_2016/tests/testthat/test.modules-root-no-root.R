## Starting with v. 2.1.2 modules do not require Root (unless, of course,
## you have a node named Root, as when using a DAG). They can use it, though.

test_that("Root not needed, but OK", {
    sa <- 0.1
    sb <- -0.2
    sab <- 0.25
    sac <- -0.1
    sbc <- 0.25
    ## OK that you use it
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
    ## But not needed
    sv3 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                                           "A : -B" = sa,
                                           "A : C" = sac,
                                           "A:B" = sab,
                                           "-A:B:C" = sbc),
                             geneToModule = c(
                                 "A" = "a1, a2",
                                 "B" = "b",
                                 "C" = "c"))
    ## Of course identical
    eagr <- evalAllGenotypes(sv2, order = FALSE, addwt = TRUE)
    eagnr <- evalAllGenotypes(sv3, order = FALSE, addwt = TRUE)
    expect_identical(eagr, eagnr)
    ## of course order does not matter
    eagr2 <- evalAllGenotypes(sv2, order = TRUE, addwt = TRUE)
    eagnr2 <- evalAllGenotypes(sv3, order = TRUE, addwt = TRUE)
    expect_identical(eagr2, eagnr2)
})

test_that("Root needed with DAGs", {
    s <- 0.2
    ## OK
    expect_silent(m1 <- allFitnessEffects(data.frame(
                      parent = c("Root", "A"),
                      child = c("A", "B"),
                      s = s,
                      sh = -1,
                      typeDep = "OR"),
                      geneToModule = c("Root" = "Root",
                                       "A" = "a1, a2",
                                       "B" = "b1")))
## But of course this fails
    expect_error(m2 <- allFitnessEffects(data.frame(
                     parent = c("Root", "A"),
                     child = c("A", "B"),
                     s = s,
                     sh = -1,
                     typeDep = "OR"),
                     geneToModule = c("A" = "a1, a2",
                                      "B" = "b1")),
                 "Some values in rT, epistasis,  or order effects not in geneToModule",
                 fixed = TRUE)
})
