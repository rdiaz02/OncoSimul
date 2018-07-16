## Copyright 2013, 2014, 2015, 2016 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## plot.evalAllGenotypes <- plot.evalAllGenotypesMut <-
##     plot.genotype_fitness_matrix <- plotFitnessLandscape 

## FIXME: show only accessible paths? 
## FIXME: when show_labels = FALSE we still show the boxes
##        and some of the labels.!!!
## FIXME: if using only_accessible, maybe we
## can try to use fast_peaks, and use the slower
## approach as fallback (if identical fitness)
plotFitnessLandscape <- function(x, show_labels = TRUE,
                                 col = c("green4", "red", "yellow"),
                                 lty = c(1, 2, 3), use_ggrepel = FALSE,
                                 log = FALSE, max_num_genotypes = 2000,
                                 only_accessible = FALSE,
                                 accessible_th = 0,
                                 ...) {

    ## FIXME future:

    ## - Allow passing order effects. Change "allGenotypes_to_matrix"
    ##   below? Probably not, as we cannot put order effects as a
    ##   matrix. Do something else, like allow only order effects if from
    ##   and allFitnessEffects object.

    ## - Allow selection only some genotypes or alleles

    ## Allow selecting only paths that involve a particular genotype in
    ## adjacency matrix of genotype i go at row i and column i.  Follow back
    ## all entries in row i, and their previous, and forward all column i.

    ## gfm: genotype fitness matrix
    ## afe: all fitness effects

    ## We seem to go back and forth, but we need to ensure the afe and the
    ## gfm are coherent. Since there are many ways to make them fail
    ## (e.g., a user passing the wrong order, or incomplete matrices, etc)
    ## we do it the following way.

    ## Yes, we set WT to 1 as we call from_genotype_fitness on
    ## matrices and genotype_fitness_matrix objects.
    ## We do that because we call evalAllGenotypes to
    ## get the string representation, etc. And this is for use
    ## with OncoSimul.

  
    tfm <- to_Fitness_Matrix(x, max_num_genotypes = max_num_genotypes)

    mutated <- rowSums(tfm$gfm[, -ncol(tfm$gfm)])
    gaj <- genot_to_adj_mat(tfm$gfm)
    if(only_accessible) {
        gaj <- filter_inaccessible(gaj, accessible_th)
        remaining <- as.numeric(colnames(gaj))
        mutated <- mutated[remaining]
        tfm$afe <- tfm$afe[remaining, , drop = FALSE]
    }
    vv <- which(!is.na(gaj), arr.ind = TRUE)
    
    ## plot(x = mutated, y = e1$Fitness, ylab = "Fitness",
    ##      xlab = "", type = "n", axes = FALSE)
    ## box()
    ## axis(2)
    ## text(x = mutated, y = x$Fitness, labels = x$Genotype)

    ## The R CMD CHEKC notes about no visible binding for global variable

    x_from <- y_from <- x_to <- y_to <- Change <- muts <-
        label <- fitness <- Type <- NULL
                
                
    dd <- data.frame(muts = mutated,
                     fitness = tfm$afe$Fitness,
                     label = tfm$afe$Genotype,
                     stringsAsFactors = TRUE)
    cl <- gaj[vv]
    sg <- data.frame(x_from = mutated[vv[, 1]],
                     y_from = tfm$afe$Fitness[vv[, 1]],
                     x_to = mutated[vv[, 2]],
                     y_to = tfm$afe$Fitness[vv[, 2]],
                     Change = factor(ifelse(cl == 0, "Neutral",
                                     ifelse(cl > 0, "Gain", "Loss")),
                                     levels = c("Gain", "Loss", "Neutral")),
                     stringsAsFactors = TRUE)
    ## From http://stackoverflow.com/a/17257422
    number_ticks <- function(n) {function(limits) pretty(limits, n)}
    
    p <- ggplot() +
        xlab("") + ylab("Fitness") + 
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              panel.grid.minor.x = element_blank()) +
        geom_segment(data = sg,
                       aes(x = x_from, y = y_from,
                           xend = x_to, yend = y_to,
                           colour = Change,
                           lty = Change)) + scale_color_manual("Change",
                                                               values = col,
                                                               drop = FALSE) +
        scale_linetype_manual("Change", values = lty, drop = FALSE) +
        scale_x_continuous(breaks = seq(0, max(mutated), 1))
    
    if(log) {
        p <- p + scale_y_log10(breaks = number_ticks(5))
    } else {
        p <- p + scale_y_continuous(breaks = number_ticks(5))
    }
    peaks_valleys <- peak_valley(gaj)
    maxF <- peaks_valleys$peak
    minF <- peaks_valleys$valley

    ddMM <- dd
    ddMM$Type <- NA
    ddMM$Type[minF] <- "Sink"
    ddMM$Type[maxF] <- "Peak"
    ddMM <- ddMM[ddMM$Type %in% c("Sink", "Peak"), ]
    
    if(show_labels && use_ggrepel) {
        p <- p + geom_text_repel(data = dd[-c(minF, maxF), ],
                                 aes(x = muts, y = fitness, label = label)) +
            geom_label_repel(data = ddMM,
                             aes(x = muts, y = fitness, label = label, fill = Type),
                             color = "black",
                             fontface = "bold")
            ## geom_label_repel(data = dd[maxF, ], aes(x = muts, y = fitness, label = label),
            ##                  color = "black",
            ##                  fontface = "bold", fill = col[1]) +
            ## geom_label_repel(data = dd[minF, ], aes(x = muts, y = fitness, label = label),
            ##                  color = "black",
            ##                  fontface = "bold", fill = col[2]) 
        
        }
    else {
        if(!show_labels) {
            ## Use same code, but make label empty
            ddMM$label <- ""
        }
            
        p <- p + geom_text(data = dd[-c(minF, maxF), ],
                           aes(x = muts, y = fitness, label = label),
                           vjust = -0.2, hjust = "inward") +
            geom_label(data = ddMM,
                       aes(x = muts, y = fitness, label = label, fill = Type),
                       vjust = -0.2, hjust = "inward", color = "black",
                       fontface = "bold")
            ## geom_label(data = dd[maxF, ], aes(x = muts, y = fitness, label = label),
            ##            vjust = -0.2, hjust = "inward", color = "black",
            ##            fontface = "bold", fill = col[1]) +
            ## geom_label(data = dd[minF, ], aes(x = muts, y = fitness, label = label),
            ##            vjust = -0.2, hjust = "inward", color = "black",
            ##            fontface = "bold", fill = col[2])
    }
    p <- p + scale_fill_manual("Local\nmax/min",  values = col)
    p
}


plot.evalAllGenotypes <- plot.evalAllGenotypesMut <-
    plot.genotype_fitness_matrix <-
        plot_fitness_landscape <-
            plotFitnessLandscape


######################################################################
######################################################################
######################################################################
####
####            Internal functions
####
######################################################################
######################################################################
######################################################################


## wrap_filter_inaccessible <- function(x, max_num_genotypes, accessible_th) {
##     ## wrap it, for my consumption
##     tfm <- to_Fitness_Matrix(x, max_num_genotypes = max_num_genotypes)
##     mutated <- rowSums(tfm$gfm[, -ncol(tfm$gfm)])
##     gaj <- genot_to_adj_mat(tfm$gfm)
##     gaj <- filter_inaccessible(gaj, accessible_th)
##     remaining <- as.numeric(colnames(gaj))
##     mutated <- mutated[remaining]
##     tfm$afe <- tfm$afe[remaining, , drop = FALSE]
##     return(list(remaining = remaining,
##                 mutated = mutated,
##                 tfm = tfm))
## }

## No longer being used. Used to be in rfitness
## count_accessible_g <- function(gfm, accessible_th) {
##     gaj <- genot_to_adj_mat(gfm)
##     gaj <- filter_inaccessible(gaj, accessible_th)
##     return(ncol(gaj) - 1)
## }


## There is now C++ code to get just the locations/positions of the
## accessible genotypes
filter_inaccessible <- function(x, accessible_th) {
    ## Return an adjacency matrix with only accessible paths. x is the gaj
    ## matrix created in the plots. The difference between genotypes
    ## separated by a hamming distance of 1

    ## FIXME: could do the x[, -1] before loop and not add the 1
    ## inside while, and do that at end
    colnames(x) <- rownames(x) <- 1:ncol(x)
    while(TRUE) {
        ## remove first column
        ## We use fact that all(na.omit(x) < u) is true if all are NA
        ## so inaccessible rows are removed and thus destination columns
        wrm <- which(apply(x[, -1, drop = FALSE], 2,
                           function(y) {all(na.omit(y) < accessible_th)})) + 1
        if(length(wrm) == 0) break;
        x <- x[-wrm, -wrm, drop = FALSE]
    }
    x[x < 0] <- NA
    return(x)
}


## wrapper to the C++ code
fast_peaks <- function(x, th) {
    ## x is the fitness matrix, not adjacency matrix

    ## Only works if no connected genotypes that form maxima. I.e., no
    ## identical fitness. Do a sufficient check for it (too inclusive, though)
    ## And only under no back mutation

    original_pos <- 1:nrow(x)
    numMut <- rowSums(x[, -ncol(x)])
    o_numMut <- order(numMut)
    x <- x[o_numMut, ]
    numMut <- numMut[o_numMut]
    f <- x[, ncol(x)]
    ## Two assumptions
    stopifnot(numMut[1] == 0)
    ## make sure no repeated in those that could be maxima
    if(any(duplicated(f[f >= f[1]])))
        stop("There could be several connected maxima genotypes.",
             " This function is not safe to use.")
    
    y <- x[, -ncol(x)]
    storage.mode(y) <- "integer"
    original_pos <- original_pos[o_numMut]
    return(sort(original_pos[peaksLandscape(y, f,
                          as.integer(numMut),
                          th)]))
}


## wrapper to the C++ code
wrap_accessibleGenotypes <- function(x, th) {
    ## x is the fitness matrix, not adjacency matrix
    original_pos <- 1:nrow(x)
    numMut <- rowSums(x[, -ncol(x)])
    o_numMut <- order(numMut)
    x <- x[o_numMut, ]
    numMut <- numMut[o_numMut]
    original_pos <- original_pos[o_numMut]
    
    y <- x[, -ncol(x)]
    storage.mode(y) <- "integer"

    acc <- accessibleGenotypes(y, x[, ncol(x)],
                               as.integer(numMut),
                               th)
    ## return(acc[acc > 0])
    return(sort(original_pos[acc[acc > 0]]))
}


## wrapper to the C++ code; the former one, only for testing. Remove
## eventually FIXME
wrap_accessibleGenotypes_former <- function(x, th) {
    ## x is the fitness matrix, not adjacency matrix
    original_pos <- 1:nrow(x)
    numMut <- rowSums(x[, -ncol(x)])
    o_numMut <- order(numMut)
    x <- x[o_numMut, ]
    numMut <- numMut[o_numMut]
    original_pos <- original_pos[o_numMut]
    
    y <- x[, -ncol(x)]
    storage.mode(y) <- "integer"

    acc <- accessibleGenotypes_former(y, x[, ncol(x)],
                                      as.integer(numMut),
                                      th)
    return(sort(original_pos[acc[acc > 0]]))
}

## A transitional function
faster_accessible_genotypes_R <- function(x, th) {
   rs0 <- rowSums(x[, -ncol(x)])
    x <- x[order(rs0), ]
    rm(rs0)
    
    y <- x[, -ncol(x)]
    f <- x[, ncol(x)]
    rs <- rowSums(y)

   ## If 0, not accessible
   ## adm <- slam::simple_triplet_zero_matrix(nrow = length(rs), ncol = length(rs),
   ##                                          mode = "integer")
   
   adm <- matrix(0, nrow = length(rs), ncol = length(rs))
   storage.mode(adm) <- "integer"
   
   ## Most time is gone here
    for(i in 1:length(rs)) { ## i is the current genotype
        candidates <- which(rs == (rs[i] + 1))
        for(j in candidates) {
            if( (sum(abs(y[j, ] - y[i, ])) == 1) &&
                ( (f[j] - f[i]) >= th ) ) {
                ## actually, this is the largest time sink using slam
                adm[i, j] <- 1L
                }
        }
    }

    colnames(adm) <- rownames(adm) <- 1:ncol(adm)
    admtmp <- adm[, -1, drop = FALSE] ## we do not want the root column.
    while(TRUE) {
        ## We remove inaccessible cols (genotypes) and the corresponding
        ## rows repeatedly until nothing left to remove; any column left
        ## is therefore accessible throw at least one path.

        ## inacc_col <- which(slam::colapply_simple_triplet_matrix(admtmp, FUN = sum) == 0L)
        inacc_col <- which(colSums(admtmp) == 0L)
        if(length(inacc_col) == 0) break;
        inacc_row <- inacc_col + 1 ## recall root row is left
        admtmp <- admtmp[-inacc_row, -inacc_col, drop = FALSE]
    }
    return(as.numeric(c(colnames(adm)[1], colnames(admtmp))))

}


## ## This uses slam, but that is actually slower because
## ## of the assignment
## faster_accessible_genots_slam <- function(x, th = 0) {

##     ## Given a genotype matrix, return the genotypes that are accessible
##     ## via creating a directed adjacency matrix between genotypes
##     ## connected (i.e., those that differ by gaining one mutation). 0
##     ## means not connected, 1 means connected.
    
##     ## There is a more general function in OncoSimulR that will give the
##     ## fitness difference. But not doing the difference is faster than
##     ## just setting a value, say 1, if all we want is to keep track of
##     ## accessible ones. And by using only 0/1 we can store only an
##     ## integer. And no na.omits, etc. Is too restricted? Yes. But for
##     ## simulations and computing just accessible genotypes, probably a
##     ## hell of a lot faster.

##     ## Well, this is not incredibly fast either.
    
##     ## Make sure sorted, so ancestors always before descendants
##     rs0 <- rowSums(x[, -ncol(x)])
##     x <- x[order(rs0), ]
##     rm(rs0)
    
##     y <- x[, -ncol(x)]
##     f <- x[, ncol(x)]
##     rs <- rowSums(y)

##     ## If 0, not accessible
##     adm <- slam::simple_triplet_zero_matrix(nrow = length(rs), ncol = length(rs),
##                                       mode = "integer")
##     for(i in 1:length(rs)) { ## i is the current genotype
##         candidates <- which(rs == (rs[i] + 1))
##         for(j in candidates) {
##             ## sumdiff <- sum(abs(y[j, ] - y[i, ]))
##             ## if(sumdiff == 1)
##             if( (sum(abs(y[j, ] - y[i, ])) == 1) &&
##                 ( (f[j] - f[i]) >= th ) )
##                 adm[i, j] <- 1L
##         }
##     }

##     colnames(adm) <- rownames(adm) <- 1:ncol(adm)
##     admtmp <- adm[, -1, drop = FALSE] ## we do not want the root column.
##     while(TRUE) {
##         ## We remove inaccessible cols (genotypes) and the corresponding
##         ## rows repeatedly until nothing left to remove; any column left
##         ## is therefore accessible throw at least one path.

##         ## inacc_col <- which(slam::colapply_simple_triplet_matrix(admtmp, FUN = sum) == 0L)
##         inacc_col <- which(slam::col_sums(admtmp) == 0L)
##         if(length(inacc_col) == 0) break;
##         inacc_row <- inacc_col + 1 ## recall root row is left
##         admtmp <- admtmp[-inacc_row, -inacc_col, drop = FALSE]
##     }
##     return(as.numeric(c(colnames(adm)[1], colnames(admtmp))))
## }


generate_matrix_genotypes <- function(g) {
    ## FIXME future: do this for order too? Only if rfitness for order.
    ## Given a number of genes, generate all possible genotypes.
    
    if(g > 20) stop("This would generate more than one million genotypes")
    ## g: number of genes
    f1 <- function(n) {
        lapply(seq.int(n), function(x) combinations(n = n, r = x))
    }
    genotNums <- f1(g)
    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}),
                      recursive = FALSE),
               function(m) m[[1]])
    }
    genotNums <- list.of.vectors(genotNums)
    
    v <- rep(0, g)
    mat <- matrix(unlist(lapply(genotNums,
                                function(x) {
                                    v[x] <- 1
                                    return(v)
                                }
                                )), ncol = g, byrow = TRUE)
    mat <- rbind(rep(0, g), mat)
    colnames(mat) <- LETTERS[1:g]
    return(mat)
}


## The R version. See also the C++ one
genot_to_adj_mat_R <- function(x) {
    ## x is the fitness matrix
    
    ## FIXME this can take about 23% of the time of the ggplot call.
    ## But them, we are quickly constructing a 2000*2000 matrix
    ## Given a genotype matrix, as given by allGenotypes_to_matrix, produce a
    ## directed adjacency matrix between genotypes connected (i.e., those
    ## that differ by gaining one mutation) with value being the
    ## difference in fitness between destination and origin

    ## FIXME: code is now in place to do all of this in C++
    
    ## Make sure sorted, so ancestors always before descendants
    original_pos <- 1:nrow(x)
    numMut <- rowSums(x[, -ncol(x)])
    o_numMut <- order(numMut)
    x <- x[o_numMut, ]
    original_pos <- original_pos[o_numMut]
    rm(numMut)
    
    y <- x[, -ncol(x)]
    f <- x[, ncol(x)]
    rs <- rowSums(y) ## redo for paranoia; could have ordered numMut

    ## Move this to C++?
    adm <- matrix(NA, nrow = length(rs), ncol = length(rs))
    for(i in 1:length(rs)) { ## i is the current genotype
        candidates <- which(rs == (rs[i] + 1))
        for(j in candidates) {
            sumdiff <- sum(abs(y[j, ] - y[i, ]))
            ## if(sumdiff < 0) stop("eh?")
            if(sumdiff == 1) adm[i, j] <- (f[j] - f[i])
        }
    }
    colnames(adm) <- rownames(adm) <- original_pos
    return(adm)
}

genot_to_adj_mat <- function(x) {
    ## x is the fitness matrix

    ## adding column and row names should rarely be necessary
    ## as these are internal functions, etc. But just in case
    original_pos <- 1:nrow(x)
    numMut <- rowSums(x[, -ncol(x)])
    o_numMut <- order(numMut)
    x <- x[o_numMut, ]
    numMut <- numMut[o_numMut]
    original_pos <- original_pos[o_numMut]
    
    y <- x[, -ncol(x)]
    storage.mode(y) <- "integer"

    adm <- genot2AdjMat(y, x[, ncol(x)],
                        as.integer(numMut))
    colnames(adm) <- rownames(adm) <- original_pos
    return(adm)
}


## ## to move above to C++ note that loop can be
## for(i in 1:length(rs)) { ## i is the current genotype
##     for(j in (i:length(rs))) {
##         if(rs[j] > (rs[i] + 1)) break;
##         else if(rs[j] == (rs[i] + 1)) {
##             ## and use here my HammingDistance function
##             ## sumdiff <- sum(abs(y[j, ] - y[i, ]))
##             ## if(sumdiff == 1) adm[i, j] <- (f[j] - f[i])
##             if(HammingDistance(y[j, ], y[i, ]) == 1) adm[i, j] = (f[j] - f[i]);
##             }
##     }
## }

## actually, all that is already in accessibleGenotypes except for the
## filling up of adm.





peak_valley <- function(x) {
    ## FIXME: when there are no identical entries, all this
    ## can be simplified a lot. Actually, there could be a way
    ## to only use the slow code on paths with 0.
    ## But this does not seem to be the CPU hog.
    ## x is an adjacency matrix where i,j is fitness_j - fitness_i.  Thus,
    ## negative values means j has less fitness than the incoming.  And if
    ## all not NA entries of a column are negative, then column j is is
    ## sink candidate
    ## Thus we locate local maxima and minima

    ## Some of this complicated as we need to detect paths like -3, 0, 0,
    ## 3.  The -3 and the two 0 are local minima with identical values.

    ## We assume crucially that ancestors always have smaller indexes in
    ## the matrix than descendants.

    ## Valleys
    bad_back <- vector("integer", nrow(x))
    for(i in ncol(x):1) {
        if(any(x[i, ] < 0, na.rm = TRUE) || bad_back[i]) {
            ## this node is bad. Any ancestor with fitness >= is bad
            bad_back[i] <- 1
            reach_b <- which(x[, i] <= 0)
            bad_back[reach_b] <- 1
        }
    }
    bad_fwd <- vector("integer", nrow(x))
    for(i in 1:nrow(x)) {
        if( any(x[, i] > 0, na.rm = TRUE) || bad_fwd[i] ) {
            ## this node is bad. Any descendant with fitness >= is bad
            bad_fwd[i] <- 1
            reach_f <- which(x[i, ] >= 0)
            bad_fwd[reach_f] <- 1
        }
    }
    bad <- union(which(bad_back == 1), which(bad_fwd == 1))
    candidate <- which(apply(x, 2, function(z) all(z <= 0, na.rm = TRUE)))
    valley <- setdiff(candidate, bad)

    null <- suppressWarnings(try({rm(bad_back, bad_fwd, reach_b, reach_f, candidate)},
        silent = TRUE))
    ## Peaks
    bad_back <- vector("integer", nrow(x))
    for(i in ncol(x):1) {
        if(any(x[i, ] > 0, na.rm = TRUE) || bad_back[i]) {
            ## this node is bad. Any ancestor with fitness >= is bad
            bad_back[i] <- 1
            reach_b <- which(x[, i] >= 0)
            bad_back[reach_b] <- 1
        }
    }
    bad_fwd <- vector("integer", nrow(x))
    for(i in 1:nrow(x)) {
        ## Eh, why any? All.
        ## Nope, any: we want peaks in general, not just
        ## under assumption of "no back mutation"
        ## We get a different result when we restrict to accessible
        ## because all < 0 in adjacency are turned to NAs.
        if( any(x[, i] < 0, na.rm = TRUE) || bad_fwd[i] ) {
        ## if( all(x[, i] < 0, na.rm = TRUE) ) {
            ## this node is bad. Any descendant with fitness >= is bad
            bad_fwd[i] <- 1
            reach_f <- which(x[i, ] <= 0)
            bad_fwd[reach_f] <- 1
        }
    }
    bad <- union(which(bad_back == 1), which(bad_fwd == 1))
    candidate <- which(apply(x, 2, function(z) all(z >= 0, na.rm = TRUE)))
    peak <- setdiff(candidate, bad)
    return(list(peak = peak, valley = valley))
}



## For the future
## ## data.frame (two columns: genotype with "," and Fitness) -> fitness graph (DAG)
## ## Return an adj matrix of the fitness graph from a fitness
## ## landscape
## ## Based on code in plotFitnessLandscape
## flandscape_to_fgraph <- function(afe) {
##     gfm <- OncoSimulR:::allGenotypes_to_matrix(afe)
##     ## mutated <- rowSums(gfm[, -ncol(gfm)])
##     gaj <- OncoSimulR:::genot_to_adj_mat(gfm)
##     gaj2 <- OncoSimulR:::filter_inaccessible(gaj, 0)
##     Eh! that is assuming no genotype is inaccessible!
##     stopifnot(all(na.omit(as.vector(gaj == gaj2))))
##     remaining <- as.numeric(colnames(gaj2))
##     ## mutated <- mutated[remaining]
##     afe <- afe[remaining, , drop = FALSE]
##     ## vv <- which(!is.na(gaj2), arr.ind = TRUE)

##     gaj2 <- gaj2
##     gaj2[is.na(gaj2)] <- 0
##     gaj2[gaj2 > 0] <- 1
##     colnames(gaj2) <- rownames(gaj2) <- afe[, "Genotype"]
##     return(gaj2)
## }
## ## This could be done easily in C++, taking care of row/colnames at end,
## ## without moving around the full adjacency matrix.
## ## Skeleton for C++
## ## a call to accessibleGenotypesPeaksLandscape
## ## (with another argument or changing the returnpeaks by a three value thing)
## ## after done with first loop, 
## ## return the matrix adm[accessible > 0, accessible >0]
## ## only need care with row/colnames
