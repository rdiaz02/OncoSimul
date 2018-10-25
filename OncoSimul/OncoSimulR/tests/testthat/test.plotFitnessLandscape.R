inittime <- Sys.time()
## This actually tests much more than plotFitnessLandscape
cat(paste("\n Starting plotFitnessLandscape at", date()))
test_that("Exercise plotting and dealing with different matrix input", {
    r1 <- rfitness(4)
    expect_silent(plot(r1))
    expect_silent(plot(r1, log = TRUE))
    expect_silent(plot(r1, log = TRUE, use_ggrepel = TRUE))
    expect_silent(plot(r1, log = TRUE, show_labels = FALSE))
    
    
    ## Specify fitness in a matrix, and plot it
    m5 <- cbind(A = c(0, 1, 0, 1), B = c(0, 0, 1, 1), F = c(1, 2, 3, 5.5))
    expect_silent(plotFitnessLandscape(m5))

    m6 <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1), c(1, 2, 3, 5.5))
    expect_message(plotFitnessLandscape(m6),
                   "No column names:", fixed = TRUE)

    ## the next are so ill formed that they should not be accepted
    m7 <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1), F = c(1, 2, 3, 5.5))
    expect_error(plotFitnessLandscape(m7),
                    "duplicated column names", fixed = TRUE)

    ## zz: why isn't this working?
    m8 <- cbind(A = c(0, 1, 0, 1), c(0, 0, 1, 1), F = c(1, 2, 3, 5.5))
    expect_warning(plotFitnessLandscape(m8),
                   "One column named ''", fixed = TRUE)

    m88 <- cbind(B = c(0, 1, 0, 1), c(0, 0, 1, 1), F = c(1, 2, 3, 5.5))
    expect_identical(as.data.frame(
        evalAllGenotypes(allFitnessEffects(genotFitness = m88),
                                      addwt = TRUE)),
                     data.frame(Genotype = c("WT", "A", "B", "A, B"),
                                Fitness = c(1, 3, 2, 5.5),
                                stringsAsFactors = FALSE))
    expect_warning(plotFitnessLandscape(m88),
                   "One column named ''", fixed = TRUE)

    
    ## Specify fitness with allFitnessEffects, and plot it
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1))
    
    expect_silent(plot(evalAllGenotypes(fe, order = FALSE)))

    ## same as
    expect_silent(plotFitnessLandscape(evalAllGenotypes(fe, order = FALSE)))
    ## more ggrepel
    expect_silent(plot(evalAllGenotypes(fe, order = FALSE), use_ggrepel = TRUE))

    m98 <- cbind(B = c(2, 1, 0, 1), c(0, 0, 1, 1), F = c(1, 2, 3, 5.5))
    expect_error(allFitnessEffects(genotFitness = m98),
                 "First ncol - 1 entries not in ",
                 fixed = TRUE)
})


test_that("to_FitnessMatrix stops as it should", {
    x1 <- data.frame(a = 1:2, b = 1:2)
    expect_error(OncoSimulR:::to_Fitness_Matrix(x1, 2000),
                 "We cannot guess what you are passing",
                 fixed = TRUE)
    x2 <- list(a = 12, b = 13)
    expect_error(OncoSimulR:::to_Fitness_Matrix(x2, 2000),
                 "We cannot guess what you are passing",
                 fixed = TRUE)
    ## This is done above
    ## g <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1))
    ## s1 <- c(1, 1.4, 1.2, 1.5)
    ## expect_error(OncoSimulR:::to_Fitness_Matrix(cbind(g, s1), 2000),
    ##              "Matrix x must have column names",
    ##              fixed = TRUE)
    ## expect_message(plotFitnessLandscape(cbind(g, s1)),
    ##              "Matrix x must have column names",
    ##              fixed = TRUE)
    ## expect_message(plotFitnessLandscape(cbind(g, A = c(1, 2))),
    ##              "Matrix x must have column names",
    ##              fixed = TRUE)
})



test_that("to_FitnessMatrix can deal with df", {
    m4 <- data.frame(G = c("A, B", "A", "WT", "B"),
                     Fitness = c(3, 2, 1, 4))
    expect_message(OncoSimulR:::to_Fitness_Matrix(m4, 2000),
                   "Column names of object", fixed = TRUE)
    m5 <- data.frame(G = c("A, B", "B"),
                     Fitness = c(3, 2))
    expect_message(OncoSimulR:::to_Fitness_Matrix(m5, 2000),
                   "Column names of object", fixed = TRUE)
    x1 <- data.frame(a = c("A, B"), Fitness = 2)
    expect_message(OncoSimulR:::to_Fitness_Matrix(x1, 2000),
                   "Column names of object", fixed = TRUE)
    x2 <- data.frame(a = c("A, B", "B"), Fitness = c(2, 3))
    expect_message(OncoSimulR:::to_Fitness_Matrix(x2, 2000),
                   "Column names of object", fixed = TRUE)
    x3 <- data.frame(a = c("A, B", "C"), Fitness = c(2, 3))
    expect_message(OncoSimulR:::to_Fitness_Matrix(x3, 2000),
                   "Column names of object", fixed = TRUE)
    ## Now, the user code
    expect_message(plotFitnessLandscape(x1))
    expect_message(plotFitnessLandscape(x2))
    expect_message(plotFitnessLandscape(x3))
    expect_message(plotFitnessLandscape(m5))
    expect_message(plotFitnessLandscape(m4))
})


test_that("internal peak valley functions", {
    
    x <- matrix(NA, 14, 14)
    x[1, 3] <- -2
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 0
    x[4, 5] <- 0
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    
    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 6, 9, 13), pv$peak)
    expect_equal(c(2, 7, 8, 14), pv$valley)

    x <- matrix(NA, 15, 15)
    x[1, 3] <- -2
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 0
    x[4, 5] <- 0
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    x[14, 15] <- 2
    
    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 6, 9, 13, 15), pv$peak)
    expect_equal(c(2, 7, 8, 14), pv$valley)


    x <- matrix(NA, 15, 15)
    x[1, 3] <- -2
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 3
    x[4, 5] <- 0
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    x[14, 15] <- 2

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 6, 9, 13, 15), pv$peak)
    expect_equal(c(2, 7, 8, 14), pv$valley)



    x <- matrix(NA, 15, 15)
    x[1, 3] <- -2
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 3
    x[4, 5] <- -1
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    x[14, 15] <- 2

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 4, 6, 9, 13, 15), pv$peak)
    expect_equal(c(2, 5, 7, 8, 14), pv$valley)


    x <- matrix(NA, 15, 15)
    x[1, 3] <- -2
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 3
    x[4, 5] <- -1
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    x[14, 15] <- 2
    x[2, 7] <- 1

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 4, 6, 9, 13, 15), pv$peak)
    expect_equal(c(2, 5, 14), pv$valley)



    x <- matrix(NA, 15, 15)
    x[1, 3] <- -2
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 0
    x[4, 5] <- -1
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    x[14, 15] <- 2
    x[2, 7] <- 1 ## hummm.. 3 and 4 should be a peak?Nope, from 1

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 6, 9, 13, 15), pv$peak)
    expect_equal(c(2, 5, 14), pv$valley)


    x <- matrix(NA, 15, 15)
    x[1, 3] <- 1
    x[1, 2] <- -4
    x[2, 3] <- 5
    x[3, 4] <- 0
    x[4, 5] <- -1
    x[5, 6] <- 4
    x[3, 7] <- -3
    x[7, 8] <- 0
    x[3, 10] <- -4
    x[10, 11] <- -4
    x[11, 12] <- 0
    x[12, 13] <- 3
    x[8, 9] <- 5
    x[12, 14] <- -5
    x[14, 15] <- 2
    x[2, 7] <- 1

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(3, 4, 6, 9, 13, 15), pv$peak)
    expect_equal(c(2, 5, 14), pv$valley)



    x <- matrix(NA, 5, 5)
    x[1, 3] <- -2
    x[2, 3] <- 4
    x[3, 4] <- 0
    x[4, 5] <- 6

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 5), pv$peak)
    expect_equal(c(2), pv$valley)


    x <- matrix(NA, 5, 5)
    x[1, 3] <- -2
    x[2, 3] <- -5
    x[3, 4] <- 0
    x[4, 5] <- 6

    (pv <-  OncoSimulR:::peak_valley(x))
    expect_equal(c(1, 2, 5), pv$peak)
    expect_equal(c(3, 4), pv$valley)

    
})


## Beware that using peak_valley on only_accessible makes a difference
test_that("internal peak valley functions w/wo inaccessible filter", {
    ## A is accessible, a peak
    ## AB is a peak if only forward. But there is no
    ## reciprocal sign epistasis here!

    ## We want peaks in general, not just
    ## under assumption of "no back mutation"?

    ## Well, no, that is not obvious with cancer progression models if we
    ## do not allow back mutations.
    
    ## We get a different result when we restrict to accessible
    ## because all < 0 in adjacency are turned to NAs.

    ## Thinking in terms of adjacency matrix, AB is not a peak if it has a
    ## positive and a negative entry in its column, because the negative
    ## entry means there is an ancestor with larger fitness.
    ## But see below for why plainly using the adjacency matrix can give bad results.
    
    ## The next matrices are all fitness matrix. Last column is fitness.
    mf1 <- rbind(
        c(0, 0, 1),
        c(1, 0, 4),
        c(0, 1, 2),
        c(1, 1, 3)
    )

    plotFitnessLandscape(mf1)
    
    expect_equal(
        OncoSimulR:::peak_valley(
                                OncoSimulR:::genot_to_adj_mat(mf1))$peak, 2)
    
    expect_equal(
            OncoSimulR:::peak_valley(
                             OncoSimulR:::filter_inaccessible(
                                              OncoSimulR:::genot_to_adj_mat(mf1), 0))$peak,
        c(2, 4))

    expect_equal(
        OncoSimulR:::fast_peaks(mf1, 0),
        c(2, 4))


    ## reorder the rows of the matrix. Affects fast_peaks, as it should
    mf1 <- rbind(
        c(1, 0, 4),
        c(0, 0, 1),
        c(1, 1, 3),
        c(0, 1, 2)
    )
    
    plotFitnessLandscape(mf1)
    ## this is not affected, since it uses, by construction, the ordered matrix
    expect_equal(
        OncoSimulR:::peak_valley(
                                OncoSimulR:::genot_to_adj_mat(mf1))$peak, 2)
    ## ditto
    expect_equal(
            OncoSimulR:::peak_valley(
                             OncoSimulR:::filter_inaccessible(
                                              OncoSimulR:::genot_to_adj_mat(mf1), 0))$peak,
        c(2, 4))
    expect_equal(
        OncoSimulR:::fast_peaks(mf1, 0),
        c(1, 3))


    
    ## filtering by inaccessible also likely gets rid of all
    ## peaks in the non-accessible part of the fitness landscape.
    ## But of course those cannot be peaks, since they are inaccessible

    mf3 <- rbind(
        c(0, 0, 0, 1),
        c(1, 0, 0, 2),
        c(0, 1, 0, 0.1),
        c(0, 0, 1, 0.3),
        c(1, 1, 0, 3),
        c(1, 0, 1, 4),
        c(0, 1, 1, 0.4),
        c(1, 1, 1, 0.2)
    )

    ## plotFitnessLandscape(mf3)
    ## BC is detected as a peak, the seventh entry
    expect_equal(OncoSimulR:::peak_valley(OncoSimulR:::genot_to_adj_mat(mf3))$peak,
                 c(5, 6, 7))

    ## recall this gives the columns of the reduced matrix, which are the former
    ## 5 and 6
    expect_equal(OncoSimulR:::peak_valley(
                                  OncoSimulR:::filter_inaccessible(
                                                   OncoSimulR:::genot_to_adj_mat(mf3), 0))$peak,
                 c(3, 4))

    ## correct indices from original matrix
    expect_equal(
        OncoSimulR:::fast_peaks(mf3, 0),
        c(5, 6))

    ## works under reorder?
    expect_equal(
        OncoSimulR:::fast_peaks(mf3[c(5, 1, 2, 3, 7, 4, 6), ], 0),
        c(1, 7))

    
  

    mf4 <- rbind(
        c(0, 0, 0, 1),
        c(1, 0, 0, 2),
        c(0, 1, 0, 0.1),
        c(0, 0, 1, 0.3),
        c(1, 1, 0, 3),
        c(1, 0, 1, 4),
        c(0, 1, 1, 0.4),
        c(1, 1, 1, 1.2)
    )

    ## plotFitnessLandscape(mf4)
    
    ## ABC is not detected as a peak, because it is not.
    ## Issue is not its accessibility, but that AC and AB have larger fitness
    ## see example with mf5
    expect_equal(OncoSimulR:::peak_valley(OncoSimulR:::genot_to_adj_mat(mf4))$peak,
                 c(5, 6))

    ## recall this gives the columns of the reduced matrix, which are the former
    ## 5 and 6
    expect_equal(OncoSimulR:::peak_valley(
                                  OncoSimulR:::filter_inaccessible(
                                                   OncoSimulR:::genot_to_adj_mat(mf4), 0))$peak,
                 c(3, 4))

    expect_equal(
        OncoSimulR:::fast_peaks(mf4, 0),
        c(5, 6))

    
    ## Now ABC is accessible
      mf5 <- rbind(
        c(0, 0, 0, 1),
        c(1, 0, 0, 2),
        c(0, 1, 0, 0.1),
        c(0, 0, 1, 0.3),
        c(1, 1, 0, 3),
        c(1, 0, 1, 4),
        c(0, 1, 1, 0.4),
        c(1, 1, 1, 3.5)
    )
 
      ## plotFitnessLandscape(mf5)
      ## plotFitnessLandscape(mf5, only_accessible = TRUE)
      
      ## But only AC is the peak, correctly
      expect_equal(OncoSimulR:::peak_valley(OncoSimulR:::genot_to_adj_mat(mf5))$peak,
                   c(6))

      ## Now, both AC and ABC are peaks
      ## columns 4 and 5 correspond to genotypes 6 and 8
      expect_equal(OncoSimulR:::peak_valley(
                                    OncoSimulR:::filter_inaccessible(
                                                     OncoSimulR:::genot_to_adj_mat(mf5), 0))$peak,
                   c(4, 5))

      expect_equal(
          OncoSimulR:::fast_peaks(mf5, 0),
          c(6, 8))

    ## AC and ABC same max fitness
      mf6 <- rbind(
        c(0, 0, 0, 1),
        c(1, 0, 0, 2),
        c(0, 1, 0, 0.1),
        c(0, 0, 1, 0.3),
        c(1, 1, 0, 3),
        c(1, 0, 1, 4),
        c(0, 1, 1, 0.4),
        c(1, 1, 1, 4)
    )
 
      ## plotFitnessLandscape(mf6)
      ## Both AC and ABC are peaks. Correctly
      expect_equal(OncoSimulR:::peak_valley(OncoSimulR:::genot_to_adj_mat(mf6))$peak,
                   c(6, 8))

      ## fast peaks should refuse to run
      expect_error(
          OncoSimulR:::fast_peaks(mf6, 0),
          "There could be several connected maxima",
          fixed = TRUE)



      ## A and AC
      mf7 <- rbind(
        c(0, 0, 0, 1),
        c(1, 0, 0, 4),
        c(0, 1, 0, 0.1),
        c(0, 0, 1, 0.3),
        c(1, 1, 0, 3),
        c(1, 0, 1, 4),
        c(0, 1, 1, 0.4),
        c(1, 1, 1, 3.4)
    )
      ## plotFitnessLandscape(mf7)
      ## Both A and AC are peaks. Correctly
      expect_equal(OncoSimulR:::peak_valley(OncoSimulR:::genot_to_adj_mat(mf7))$peak,
                   c(2, 6))

      ## fast peaks should refuse to run
      expect_error(
          OncoSimulR:::fast_peaks(mf7, 0),
          "There could be several connected maxima",
          fixed = TRUE)
      

      
      ## A, AC, ABC same max fitness
      mf8 <- rbind(
        c(0, 0, 0, 1),
        c(1, 0, 0, 4),
        c(0, 1, 0, 0.1),
        c(0, 0, 1, 0.3),
        c(1, 1, 0, 3),
        c(1, 0, 1, 4),
        c(0, 1, 1, 0.4),
        c(1, 1, 1, 4)
    )
      ## plotFitnessLandscape(mf8)
      ## Both A and AC are peaks. Correctly
      expect_equal(OncoSimulR:::peak_valley(OncoSimulR:::genot_to_adj_mat(mf8))$peak,
                   c(2, 6, 8))

      
      ## fast peaks should refuse to run
      expect_error(
          OncoSimulR:::fast_peaks(mf8, 0),
          "There could be several connected maxima",
          fixed = TRUE)
      
      ## A, AC, AB same max fitness
      mf9 <- rbind(
        c(0, 0, 0, 1),
        c(1, 0, 0, 4),
        c(0, 1, 0, 0.1),
        c(0, 0, 1, 0.3),
        c(1, 1, 0, 4),
        c(1, 0, 1, 4),
        c(0, 1, 1, 0.4),
        c(1, 1, 1, 2.4)
    )
      ## plotFitnessLandscape(mf9, use_ggrepel = TRUE)
      ## Both A and AC are peaks. Correctly
      expect_equal(OncoSimulR:::peak_valley(OncoSimulR:::genot_to_adj_mat(mf9))$peak,
                   c(2, 5, 6))

      
      
      ## This illustrates that the "filter_inaccessible" is not just "do
      ## not take into account inaccessible genotypes" but, properly, do
      ## not take into account, do not allow any travelling through
      ## inaccessible paths.

      ## Thus, filter_inaccessible is the way to go if we want to exclude
      ## backmutation. In no bakcmutation, it is not possible to go from
      ## m+1 to m mutations.
      
      ## It also shows that naively looking at the adjacency matrix can
      ## fail. Two reasons:

      ## a) the last row will never have any entries and yet it need not
      ## be a peak.

      ## b) simply looking at adjacency matrix is not the correct
      ## procedure when some fitnesses can be equal. That is what the
      ## function peak_valley works hard to get right :-)
      
      
    cp2 <- structure(c(0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 
1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 
1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 
0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 
1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 
1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 
1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 
1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 
0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 
0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 
0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 
0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 
1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 
1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 
1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 
0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 
0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 
1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 
1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 
1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 
1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 
0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 
1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 
1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 
0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 
0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 
1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 
0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 
0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 
1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 
0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 
1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 
1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0.852873407703003, 
1.51520969989942, 1.09934414414554, 1.08362391548151, 1.06352377058758, 
0.875558455823467, 1.69351291065104, 2.92492684398312, 1.02057836095586, 
0.994559647972076, 1.01807462848707, 0.782398502758159, 0.755318352952028, 
1.81553780643735, 1.7427209966816, 1.00116069406198, 0.790243245268257, 
3.38168029927183, 1.18573953796889, 1.24679706264807, 0.944183929293486, 
1.04153712305771, 1.20232261789798, 1.0345783487807, 1.04678440594199, 
0.993244793867836, 0.97914067773803, 0.79321495112376, 0.868101325153957, 
0.866235177920767, 4.1155779007473, 3.163209721772, 4.34977195536485, 
1.09932137400121, 1.08612305022998, 0.916953742980573, 0.850115441923501, 
1.06277833622263, 0.865087563773651, 0.928169473201598, 0.904902930158639, 
0.897493717866434, 0.71149600120298, 1.06538015204221, 1.07859259299858, 
0.858803230350538, 2.25551012930227, 1.09241633274047, 0.870425423271033, 
2.17687545546796, 0.84459090869647, 4.58149975106353, 3.85245245151455, 
1.28342034151899, 1.08529050597462, 1.02256835452167, 1.04982916832593, 
1.0457848642841, 0.90107628754529, 1.08969768294891, 1.05766476796899, 
0.902394628842996, 0.888348932462492, 1.01037474862489, 0.954093541062801, 
0.807820459139572, 2.74832174163312, 1.01318977068049, 0.854004033396404, 
0.842034005421367, 0.800544915243185, 5.31108977064723, 5.31423066433053, 
1.16539625099584, 0.983449927610599, 0.996320237843515, 0.9794158873742, 
1.02038748073625, 0.808875731463122, 0.964868528161141, 0.966566509486774, 
0.860373057266184, 0.81168825662344, 1.19978481918247, 0.98157798351476, 
0.999463234369357, 0.98711106267367, 0.961995700808845, 4.79391503400402, 
0.998909701750288, 0.996465768481649, 0.785688019266101, 0.778917380394268, 
1.17230915723272, 1.19911647477422, 0.961939861987872, 0.981542927739855, 
0.999822362533057, 1.15236749698624, 0.919688401637553, 0.876733220798505, 
0.92069327916386, 0.958801043337062, 0.670589798279379, 0.84152795885645, 
5.93895353544503, 0.723329951949942, 0.733188455582477, 1.07557023464861, 
1.09180382079188, 0.923957719945906, 0.93313538716072, 0.896562810368268, 
1.09769821865825, 1.10615389985864, 0.94426955155254, 0.898545873061366, 
0.876269943340891, 1.11556411094416, 0.94930544641744, 1.02495854041569, 
0.794907983845338, 0.847332095413669, 0.776896984008625, 0.928896557877041, 
0.945135371172636, 0.892100531723894), .Dim = c(128L, 8L), .Dimnames = list(
    NULL, c("CDKN2A", "KRAS", "MLL3", "PXDN", "SMAD4", "TGFBR2", 
    "TP53", "")))

    expect_equal(length(
        OncoSimulR:::peak_valley(OncoSimulR:::genot_to_adj_mat(cp2))$peak), 4)

    expect_equal(length(
        OncoSimulR:::peak_valley(
                         OncoSimulR:::filter_inaccessible(
                                          OncoSimulR:::genot_to_adj_mat(cp2), 0))$peak), 6)


    expect_equal(
          OncoSimulR:::fast_peaks(cp2, 0),
        c(51, 55, 68, 74, 90, 107))

    ## Nope, since filter inaccessible removes genotypes
    expect_false(all(
        OncoSimulR:::peak_valley(
                         OncoSimulR:::filter_inaccessible(
                                          OncoSimulR:::genot_to_adj_mat(cp2), 0))$peak ==
        OncoSimulR:::fast_peaks(cp2, 0)))
    

    ## compare with the probl
    gnn <- OncoSimulR:::to_Fitness_Matrix(cp2, 1000)$afe[, "Genotype"]

    plotFitnessLandscape(cp2, use_ggrepel = TRUE, only_accessible = TRUE)
    
    expect_equal(
        gnn[OncoSimulR:::fast_peaks(cp2, 0)],
        c("KRAS, PXDN, TP53",
          "MLL3, PXDN, SMAD4",
          "CDKN2A, KRAS, MLL3, TP53",
          "CDKN2A, KRAS, TGFBR2, TP53",
          "KRAS, MLL3, TGFBR2, TP53",
          "CDKN2A, KRAS, PXDN, SMAD4, TP53"))

    ## can also check by removing the inacessible genotypes so the indices are the same
    agg <- OncoSimulR:::wrap_accessibleGenotypes(cp2, 0)
    cp3 <- cp2[agg, ]

    ## This is NOT correct: we have removed the inacessible,
    ## but we allow backmutation
    ## OncoSimulR:::peak_valley(
    ##                  OncoSimulR:::genot_to_adj_mat(cp3))$peak
    
    expect_equal(OncoSimulR:::peak_valley(
                     OncoSimulR:::filter_inaccessible(
                                      OncoSimulR:::genot_to_adj_mat(cp3), 0))$peak,
                 OncoSimulR:::fast_peaks(cp3, 0))

    gnn3 <- gnn[agg]

    expect_equal(
        gnn3[OncoSimulR:::fast_peaks(cp3, 0)],
        c("KRAS, PXDN, TP53",
          "MLL3, PXDN, SMAD4",
          "CDKN2A, KRAS, MLL3, TP53",
          "CDKN2A, KRAS, TGFBR2, TP53",
          "KRAS, MLL3, TGFBR2, TP53",
          "CDKN2A, KRAS, PXDN, SMAD4, TP53"))

    
})



test_that("Some random checks of the fast peaks function", {
    niter <- 50
    for(i in 1:niter) {
        for(ng in 2:6) {
            rtmp <- rfitness(ng)
            p1 <- OncoSimulR:::peak_valley(
                                   OncoSimulR:::filter_inaccessible(
                                                    OncoSimulR:::genot_to_adj_mat(rtmp), 0))$peak
            expect_equal(length(p1),
                         length(OncoSimulR:::fast_peaks(rtmp, 0)))
            agg <- OncoSimulR:::wrap_accessibleGenotypes(rtmp, 0)
            if(length(agg) >= 2) {
                ## cat(".")
                p2 <- OncoSimulR:::peak_valley(
                                       OncoSimulR:::filter_inaccessible(
                                                        OncoSimulR:::genot_to_adj_mat(rtmp[agg, , drop = FALSE]), 0))$peak
                expect_equal(p2, OncoSimulR:::fast_peaks(rtmp[agg, , drop = FALSE], 0))
                expect_equal(length(p2), length(p1))
            }
        }
    }
})
cat(paste("\n Ending plotFitnessLandscape at", date()), "\n")
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
