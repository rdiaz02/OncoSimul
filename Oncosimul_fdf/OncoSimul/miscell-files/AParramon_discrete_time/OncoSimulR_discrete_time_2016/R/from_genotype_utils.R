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



## Functions that allow passing a matrix or data frame of mappings
## genotype -> fitness so this is taken as input in fitnessEffects.
## and some related functions



to_Magellan <- function(x, file,
                        max_num_genotypes = 2000) {
    if(is.null(file)) {
        file <- tempfile()
        cat("\n Using file ", file, "\n")
    }
    gfm <- to_Fitness_Matrix(x, max_num_genotypes = max_num_genotypes)$gfm
    write(rep(2, ncol(gfm) - 1), file = file, ncolumns = ncol(gfm) - 1)
    write.table(gfm, file = file, append = TRUE,
                row.names = FALSE, col.names = FALSE, sep = " ")
}

to_Fitness_Matrix <- function(x, max_num_genotypes) {
    ## A general converter. Ready to be used by plotFitnessLandscape and
    ## Magellan exporter.
    
    if( (inherits(x, "genotype_fitness_matrix")) ||
        ( (is.matrix(x) || is.data.frame(x)) && (ncol(x) > 2) ) ) {
        ## Why this? We go back and forth twice. We need both things. We
        ## could construct the afe below by appropriately pasting the
        ## columns names
        ## if( (is.null(colnames(x))) || any(grepl("^$", colnames(x))))
        ##    stop("Matrix x must have column names")
        afe <- evalAllGenotypes(allFitnessEffects(
            epistasis = from_genotype_fitness(x)),
            order = FALSE, addwt = TRUE, max = max_num_genotypes)
        ## Might not be needed with the proper gfm object (so gmf <- x)
        ## but is needed if arbitrary matrices.
        gfm <- allGenotypes_to_matrix(afe) 
    } else if(inherits(x, "fitnessEffects")) {
        if(!is.null(x$orderEffects) )
            stop("We cannot yet deal with order effects")
        afe <- evalAllGenotypes(x,
                                order = FALSE,
                                addwt = TRUE, max = max_num_genotypes)
        gfm <- allGenotypes_to_matrix(afe)
    } else if( (inherits(x, "evalAllGenotypes")) ||
               (inherits(x, "evalAllGenotypesMut"))) {
        if(any(grepl(">", x[, 1], fixed = TRUE)))
            stop("We cannot deal with order effects yet.")
        x <- x[, c(1, 2)]
        if(x[1, "Genotype"] != "WT") {
            ## Yes, because we expect this present below
            x <- rbind(data.frame(Genotype = "WT",
                                  Fitness = 1,
                                  stringsAsFactors = FALSE),
                       x)
        }
        afe <- x
        ## in case we pass an evalAllgenotypesfitandmut
        gfm <- allGenotypes_to_matrix(afe)
    } else if(is.data.frame(x)) {
        ## Assume a two-column data frame of genotypes as character
        ## vectors and fitness
        if(colnames(x)[2] != "Fitness")
            stop("We cannot guess what you are passing here") 
        afe <- evalAllGenotypes(allFitnessEffects(genotFitness = x),
                                order = FALSE, addwt = TRUE,
                                max = max_num_genotypes)
        gfm <- allGenotypes_to_matrix(afe)
    } else {
        stop("We cannot guess what you are passing here") 
    }
    return(list(gfm = gfm, afe = afe))
}   


from_genotype_fitness <- function(x) {
    ## Would break with output from allFitnessEffects and
    ## output from allGenotypeAndMut
    if(ncol(x) > 2) {
        ## Of course, can ONLY work with epistastis, NOT order
        return(genot_fitness_to_epistasis(x))
    } else {
        ## Make sure no factors
        if(is.factor(x[, 1])) x[, 1] <- as.character(x[, 1])
        omarker <- any(grepl(">", x[, 1], fixed = TRUE))
        emarker <- any(grepl(",", x[, 1], fixed = TRUE))
        nogoodepi <- any(grepl(":", x[, 1], fixed = TRUE))
        ## if(omarker && emarker) stop("Specify only epistasis or order, not both.")
        if(nogoodepi && emarker) stop("Specify the genotypes separated by a ',', not ':'.")
        if(nogoodepi && !emarker) stop("Specify the genotypes separated by a ',', not ':'.")
        ## if(nogoodepi && omarker) stop("If you want order, use '>' and if epistasis ','.")
        ## if(!omarker && !emarker) stop("You specified neither epistasis nor order")
        if(omarker) {
            ## do something. To be completed
            stop("This code not yet ready")
            ## You can pass to allFitnessEffects genotype -> fitness mappings that
            ## involve epistasis and order. But they must have different
            ## genes. Otherwise, it is not manageable.
        }
        if(emarker) {
            x <- x[, c(1, 2)]
            if(!all(colnames(x) == c("Genotype", "Fitness"))) {
                message("Column names of object not Genotype and Fitness.",
                        " Renaming them assuming that is what you wanted")
                colnames(x) <- c("Genotype", "Fitness")
            }
            ## Yes, we need to do this to  scale the fitness and put the "-"
            return(genot_fitness_to_epistasis(allGenotypes_to_matrix(x)))
        }
        
    }
}



genot_fitness_to_epistasis <- function(x) {
    ## FIXME future:

    ## - Nope, order cannot be dealt with here. Not a matrix of 0 and 1.

    ## - modify "fitnessEffects" so we can take a component that is
    ## - "genot_fitness"; so this would never be exposed to the user


    ## Why we should not combine this specification with other terms? If
    ## you use this is because you say "this is the mapping genotype ->
    ## fitness. Period." so we should not combine other terms (or other
    ## terms that involve these genes)
    
    nr <- nrow(x)
    if(nr < (2^(ncol(x) - 1)))
        message("Number of genotypes less than 2^L.",
                " Missing genotype will be set to 1")
    ## This is specific if only epistasis, not order
    if(nr > (2^(ncol(x) - 1)))
        stop("Number of genotypes > 2^L. Repeated entries?")
    f <- x[, ncol(x)]
    ## Why should I stop?
    if(any(f < 0))
        message("Negative fitnesses. Watch out if you divide by the wildtype")
    x <- x[, -ncol(x)]
    wt <- which(rowSums(x) == 0)
    fwt <- 1
    if(length(wt) == 1)
        fwt <- f[wt]
    if(!isTRUE(all.equal(fwt, 1))) {
        message("Fitness of wildtype != 1. ",
                "Dividing all fitnesses by fitness of wildtype.")
        f <- f/fwt
    }
    
    if(is.null(colnames(x)) || any(grepl("^$", colnames(x))) ) {
        message("Setting/resetting gene names because one or more are missing.",
                " If this is not what you want, pass a matrix",
                " with all columns named.")
        if(ncol(x) <= 26)
            colnames(x) <- LETTERS[1:ncol(x)]
        else
            colnames(x) <- paste0("G", seq.int(ncol(x)))
    }
    cn <- colnames(x)
    x2 <- matrix("", nrow = nrow(x), ncol = ncol(x))
    x2[x == 0] <- "-"
    epin <- apply(x2, 1, function(z) paste(paste0(z, cn), collapse = ":"))
    if(anyDuplicated(epin))
        stop("Non unique names")
    s <- f - 1
    names(s) <- epin
    return(s)
}




allGenotypes_to_matrix <- function(x) { 
    ## Makes no sense to allow passing order: the matrix would have
    ## repeated rows. A > B and B > A both have exactly A and B
    
    ## Take output of evalAllGenotypes or identical data frame and return
    ## a matrix with 0/1 in a column for each gene and a final column of
    ## Fitness

    ## A WT can be specified with string "WT"
    anywt <- which(x[, 1] == "WT")
    if(length(anywt) > 1) stop("More than 1 WT")
    if(length(anywt) == 1) {
        fwt <- x[anywt, 2]
        x <- x[-anywt, ]
    } else {
        fwt <- 1
    }
    splitted_genots <- lapply(x$Genotype,
                             function(z) nice.vector.eo(z, ","))

    all_genes <- unique(unlist(splitted_genots))

    m <- matrix(0, nrow = length(splitted_genots), ncol = length(all_genes))
    the_match <- lapply(splitted_genots,
                        function(z) match(z, all_genes))
    ## A lot simpler with a loop
    for(i in 1:nrow(m)) {
        m[i, the_match[[i]]] <- 1
    }
    m <- cbind(m, x[, 2])
    colnames(m) <- c(all_genes, "Fitness")
    m <- rbind(c(rep(0, length(all_genes)), fwt),
               m)
    ## Ensure sorted
    m <- data.frame(m)
    rs <- rowSums(m[, -ncol(m)])
    m <- m[order(rs), ]
    ## m <- m[do.call(order, as.list(cbind(rs, m[, -ncol(m)]))), ]
    return(m)
}


## Magellan_stats and Magellan_draw cannot be tested
## routinely, as they depend on external software
Magellan_stats <- function(x, max_num_genotypes = 2000,
                           verbose = FALSE,
                           use_log = TRUE,
                           short = TRUE,
                           fl_statistics = "fl_statistics") { # nocov start
    ## I always use 
    ## if(!is.null(x) && is.null(file))
    ##     stop("one of object or file name")
    ## if(is.null(file))
    fn <- tempfile()
    fnret <- tempfile()
    if(verbose)
        cat("\n Using input file", fn, " and output file ", fnret, "\n")
    
    if(use_log) {
        logarg <- "-l"
    } else {
        logarg <- NULL
    }
    if(short) {
        shortarg <- "-s"
    } else {
        shortarg <- NULL
    }
    
    to_Magellan(x, fn, max_num_genotypes = max_num_genotypes)
    call_M <- system2(fl_statistics, args = paste(fn, shortarg, logarg, "-o", fnret))
    if(short) {
        tmp <- as.vector(read.table(fnret, skip = 1, header = TRUE)[-1])
        ## Make names more explicit, but check we have what we think we have
        stopifnot(length(tmp) == 23)
        stopifnot(identical(names(tmp),
                            c("ngeno", "npeaks", "nsinks", "gamma", "gamma.", "r.s",
                              "nchains", "nsteps", "nori", "depth", "magn", "sign",
                              "rsign", "w.1.", "w.2.", "w.3..", "mode_w", "s.1.",
                              "s.2.", "s.3..", "mode_s", "pairs_s", "outD_v")))
        names(tmp) <- c("n_genotypes", "n_peaks", "n_sinks", "gamma", "gamma_star",
                        "r/s","n_chains", "n_steps", "n_origins", "max_depth",
                        "epist_magn", "epist_sign", "epist_rsign",
                        "walsh_coef_1", "walsh_coef_2", "walsh_coef_3", "mode_walsh",
                        "synerg_coef_1", "synerg_coef_2", "synerg_coef_3", "mode_synerg",
                        "std_dev_pairs", "var_outdegree")
    } else {
        tmp <- readLines(fnret)
    }
    return(tmp)
} # nocov end

Magellan_draw <- function(x, max_num_genotypes = 2000,
                          verbose = FALSE,
                          args = "-f",
                          fl_draw = "fl_draw",
                          svg_open = "xdg-open",
                          file_name = NULL) { # nocov start
    ## It always works by appending the name so file_name is without the .svg
    if(is.null(file_name)) {
        fn <- tempfile()
    } else {
        fn <- file_name
    }
    fn_out <- paste0(fn, ".svg")
    if(verbose)
        cat("\n Using input file", fn, " and output file ", fn_out, "\n")
       
    to_Magellan(x, fn, max_num_genotypes = max_num_genotypes)
    call_M <- system2(fl_draw, args = paste(fn, args), wait = FALSE)
    call_view <- system2(svg_open, args = fn_out, wait = FALSE,
                         stdout = ifelse(verbose, "", FALSE),
                         stderr = ifelse(verbose, "", FALSE))
    
    invisible() 
} # nocov end



## ## Example of Bozic issues in conversions of fitnes
## m1 <- cbind(c(0, 1), c(1, 0), c(2, 3))
## m2 <- cbind(c(1, 1), c(1, 0), c(2, 3))
## m3 <- data.frame(G = c("A, B", "A"), F = c(1, 2))
## m4 <- data.frame(G = c("A, B", "A", "WT", "B"), F = c(3, 2, 1, 4))
## m5 <- data.frame(G = c("A, B", "A", "WT", "B"), F = c(3, 2, 1, 0))
## m6 <- data.frame(G = c("A, B", "A", "WT", "B"), F = c(3, 2.5, 2, 0))

## And no, it makes no sense to use any of this for mutator: in mutator I
## directly have the multiplication factor of each gene. Which is likely
## what people want anyway. Add it later if needed by using a ratio
## instead of a "-"
