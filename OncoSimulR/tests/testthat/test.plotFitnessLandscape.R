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
                   "Setting/resetting gene names because", fixed = TRUE)

    m7 <- cbind(c(0, 1, 0, 1), c(0, 0, 1, 1), F = c(1, 2, 3, 5.5))
    expect_message(plotFitnessLandscape(m7),
                   "Setting/resetting gene names because", fixed = TRUE)

    m8 <- cbind(A = c(0, 1, 0, 1), c(0, 0, 1, 1), F = c(1, 2, 3, 5.5))
    expect_message(plotFitnessLandscape(m8),
                   "Setting/resetting gene names because", fixed = TRUE)

    
    ## Specify fitness with allFitnessEffects, and plot it
    fe <- allFitnessEffects(epistasis = c("a : b" = 0.3,
                                          "b : c" = 0.5),
                            noIntGenes = c("e" = 0.1))
    
    expect_silent(plot(evalAllGenotypes(fe, order = FALSE)))

    ## same as
    expect_silent(plotFitnessLandscape(evalAllGenotypes(fe, order = FALSE)))
    ## more ggrepel
    expect_silent(plot(evalAllGenotypes(fe, order = FALSE), use_ggrepel = TRUE))
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
    ## A is inaccessible, a peak
    ## AB is a peak if only forward. But there is no
    ## reciprocal sign epistasis here!

    ## We want peaks in general, not just
    ## under assumption of "no back mutation"
    ## We get a different result when we restrict to accessible
    ## because all < 0 in adjacency are turned to NAs.
    
    mf1 <- rbind(c(0, 0, 1),
             c(1, 0, 4),
             c(0, 1, 2),
             c(1, 1, 3))

    expect_equal(
        length(OncoSimulR:::peak_valley(
                                OncoSimulR:::genot_to_adj_mat(mf1))$peak), 1)
    
    expect_equal(
        length(
            OncoSimulR:::peak_valley(
                             OncoSimulR:::filter_inaccessible(
                                              OncoSimulR:::genot_to_adj_mat(mf1), 0))$peak), 2)


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
   
})


