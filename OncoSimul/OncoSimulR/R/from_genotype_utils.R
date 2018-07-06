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
    ## Go directly if you have, as entry, an object from
    ## rfitness!! to_Fitness_Matrix can be very slow.
    if(is.null(file)) {
        file <- tempfile()
        cat("\n Using file ", file, "\n")
    }
    if(inherits(x, "genotype_fitness_matrix")) {
        write(rep(2, ncol(x) - 1), file = file, ncolumns = ncol(x) - 1)
        write.table(x, file = file, append = TRUE,
                    row.names = FALSE, col.names = FALSE, sep = " ")
    } else {
        gfm <- to_Fitness_Matrix(x, max_num_genotypes = max_num_genotypes)$gfm
        write(rep(2, ncol(gfm) - 1), file = file, ncolumns = ncol(gfm) - 1)
        write.table(gfm, file = file, append = TRUE,
                    row.names = FALSE, col.names = FALSE, sep = " ")
    }
}



## ## genotype_fitness_matrix -> fitness landscape as data frame
## fmatrix_to_afe <- function(x) {
##     stopifnot(inherits(x, "genotype_fitness_matrix"))
##     y <- x[, -ncol(x)]
##     nn <- apply(y, 1,
##                 function(u) paste(sort(colnames(y)[as.logical(u)]),
##                                   collapse = ", "))
##     nn[nn == ""] <- "WT"
##     return(data.frame(Genotype = nn, Fitness = x[, ncol(x)],
##            stringsAsFactors = FALSE))
## }

to_Fitness_Matrix <- function(x, max_num_genotypes) {
    ## A general converter. Ready to be used by plotFitnessLandscape and
    ## Magellan exporter.

    ## FIXME: really, some of this is inefficient. Very. Fix it.
    if( (inherits(x, "genotype_fitness_matrix")) ||
        ( (is.matrix(x) || is.data.frame(x)) && (ncol(x) > 2) ) ) {
        ## Why this? We go back and forth twice. We need both things. We
        ## could construct the afe below by appropriately pasting the
        ## columns names
        ## if( (is.null(colnames(x))) || any(grepl("^$", colnames(x))))
        ##    stop("Matrix x must have column names")

        ## This could use fmatrix_to_afe, above!!!
        ## Major change as of flfast: no longer using from_genotype_fitness
        afe <- evalAllGenotypes(allFitnessEffects(
            genotFitness = x
            ##, epistasis = from_genotype_fitness(x)
        ),
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

## Based on from_genotype_fitness
## but if we are passed a fitness landscapes as produced by
## rfitness, do nothing

##New function
from_letters_to_vector_genotype <- function(letters_genotype) {

  if(!all(unlist(strsplit(letters_genotype, ', ' ) ) %in% c(LETTERS,  letters, "")))
    stop("All elements must be letters separated by ', '")

  letters_genotype <-
    paste(sort(unlist(strsplit(toupper(letters_genotype), ',') ) ),
          collapse = ',')

  if (letters_genotype == "") {
    vector_genotype <- c()
  }else{
    vector_genotype <- sapply(unlist(strsplit(letters_genotype, ', ') ),
                              function(x) {which(x == LETTERS ) },
                              USE.NAMES = FALSE  )
  }
  return (as.vector(vector_genotype))
}

##New function
replaceWithNumbers <- function(str){
  locs <- gregexpr("f\\(([aA-zZ],? ?)*\\)", str)
  regmatches(str, locs) <- list(sapply(regmatches(str, locs)[[1]], function(x) {
    paste0("f(",
           toString(from_letters_to_vector_genotype(sub(".*\\((.*)\\).*",
                                                        "\\1",
                                                        x))),
           ")")
  }))
  return(str)
}

##Modified
to_genotFitness_std <- function(x,
                                frequencyDependentFitness,
                                simplify = TRUE,
                                min_filter_fitness = 1e-9,
                                sort_gene_names = TRUE) {
  ## Would break with output from allFitnessEffects and
  ## output from allGenotypeAndMut

  ## For the very special and weird case of
  ## a matrix but only a single gene so with a 0 and 1
  ## No, this is a silly and meaningless case.
  ## if( ( ncol(x) == 2 ) && (nrow(x) == 1) && (x[1, 1] == 1) ) {

  ## } else  blabla:

  if(! (inherits(x, "matrix") || inherits(x, "data.frame")) )
    stop("Input must inherit from matrix or data.frame.")

  ## if((ncol(x) > 2) && !(inherits(x, "matrix"))
  ##     stop(paste0("Genotype fitness input either two-column data frame",
  ##          " or a numeric matrix with > 2 columns."))
  ## if( (ncol(x) > 2) && (nrow(x) == 1) )
  ##     stop(paste0("It looks like you have a matrix for a single genotype",
  ##                 " of a single gene. For this degenerate cases use",
  ##                 " a data frame specification."))

  if (!frequencyDependentFitness){
    if(ncol(x) > 2) {
      if(inherits(x, "matrix")) {
        if(!is.numeric(x))
          stop("A genotype fitness matrix/data.frame must be numeric.")
      } else if(inherits(x, "data.frame")) {
        if(!all(unlist(lapply(x, is.numeric))))
          stop("A genotype fitness matrix/data.frame must be numeric.")
      }

      ## We are expecting here a matrix of 0/1 where columns are genes
      ## except for the last column, that is Fitness
      ## Of course, can ONLY work with epistastis, NOT order
      ## return(genot_fitness_to_epistasis(x))
      if(any(duplicated(colnames(x))))
        stop("duplicated column names")

      cnfl <- which(colnames(x)[-ncol(x)] == "")
      if(length(cnfl)) {
        freeletter <- setdiff(LETTERS, colnames(x))[1]
        if(length(freeletter) == 0) stop("Renaiming failed")
        warning("One column named ''. Renaming to ", freeletter)
        colnames(x)[cnfl] <- freeletter
      }
      if(!is.null(colnames(x)) && sort_gene_names) {
        ncx <- ncol(x)
        cnx <- colnames(x)[-ncx]
        ocnx <- gtools::mixedorder(cnx)
        if(!(identical(cnx[ocnx], cnx))) {
          message("Sorting gene column names alphabetically")
          x <- cbind(x[, ocnx, drop = FALSE], Fitness = x[, (ncx)])
        }
      }

      if(is.null(colnames(x))) {
        ncx <- (ncol(x) - 1)
        message("No column names: assigning gene names from LETTERS")
        if(ncx > length(LETTERS))
          stop("More genes than LETTERS; please give gene names",
               " as you see fit.")
        colnames(x) <- c(LETTERS[1:ncx], "Fitness")
      }

      if(!all(as.matrix(x[, -ncol(x)]) %in% c(0, 1) ))
        stop("First ncol - 1 entries not in {0, 1}.")
    } else {
      if(!inherits(x, "data.frame"))
        stop("genotFitness: if two-column must be data frame")
      ## Make sure no factors
      if(is.factor(x[, 1])) {
        warning("First column of genotype fitness is a factor. ",
                "Converting to character.")
        x[, 1] <- as.character(x[, 1])
      }
      ## Make sure no numbers
      if(any(is.numeric(x[, 1])))
        stop(paste0("genotFitness: first column of data frame is numeric.",
                    " Ambiguous and suggests possible error. If sure,",
                    " enter that column as character"))

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
      if( emarker || ( (!omarker) && (!emarker) && (!nogoodepi)) ) {
        ## the second case above corresponds to passing just single letter genotypes
        ## as there is not a single marker
        x <- x[, c(1, 2), drop = FALSE]
        if(!all(colnames(x) == c("Genotype", "Fitness"))) {
          message("Column names of object not Genotype and Fitness.",
                  " Renaming them assuming that is what you wanted")
          colnames(x) <- c("Genotype", "Fitness")
        }
        if((!omarker) && (!emarker) && (!nogoodepi)) {
          message("All single-gene genotypes as input to to_genotFitness_std")
        }
        ## Yes, we need to do this to  scale the fitness and put the "-"
        x <- allGenotypes_to_matrix(x)
      }
    }
    ## And, yes, scale all fitnesses by that of the WT
    whichroot <- which(rowSums(x[, -ncol(x), drop = FALSE]) == 0)
    if(length(whichroot) == 0) {
      warning("No wildtype in the fitness landscape!!! Adding it with fitness 1.")
      x <- rbind(c(rep(0, ncol(x) - 1), 1), x)
    } else if(x[whichroot, ncol(x)] != 1) {
      warning("Fitness of wildtype != 1.",
              " Dividing all fitnesses by fitness of wildtype.")
      vwt <- x[whichroot, ncol(x)]
      x[, ncol(x)] <- x[, ncol(x)]/vwt
    }
    if(any(is.na(x)))
      stop("NAs in fitness matrix")
    if(simplify) {
      return(x[x[, ncol(x)] > min_filter_fitness, , drop = FALSE])
    } else {
      return(x)
    }
  }

  if (frequencyDependentFitness){

    if(inherits(x, "matrix")) {
      if(!is.numeric(x[-ncol(x)]))
        stop("All columns except the last one must be numeric.")
    } else if(inherits(x, "data.frame")) {
      if(!all(unlist(lapply(x[-ncol(x)], is.numeric))))
        stop("All columns except the last one must be numeric.")
    }

    if(!all(unlist(lapply(x[ncol(x)], is.character))))
      stop("All elements in last column must be character.")

    if(any(duplicated(colnames(x))))
      stop("duplicated column names")

    cnfl <- which(colnames(x)[-ncol(x)] == "")
    if(length(cnfl)) {
      freeletter <- setdiff(LETTERS, colnames(x))[1]
      if(length(freeletter) == 0) stop("Renaiming failed")
      warning("One column named ''. Renaming to ", freeletter)
      colnames(x)[cnfl] <- freeletter
    }

    if(!is.null(colnames(x)) && sort_gene_names) {
      ncx <- ncol(x)
      cnx <- colnames(x)[-ncx]
      ocnx <- gtools::mixedorder(cnx)
      if(!(identical(cnx[ocnx], cnx))) {
        message("Sorting gene column names alphabetically")
        x <- cbind(x[, ocnx, drop = FALSE], Fitness = x[, (ncx)])
      }
    }

    if(is.null(colnames(x))) {
      ncx <- (ncol(x) - 1)
      message("No column names: assigning gene names from LETTERS")
      if(ncx > length(LETTERS))
        stop("More genes than LETTERS. Not supported by now.")
      colnames(x) <- c(LETTERS[1:ncx], "Fitness")
    }

    if(!all(as.matrix(x[, -ncol(x)]) %in% c(0, 1) ))
      stop("First ncol - 1 entries not in {0, 1}.")

    if(any(is.na(x)))
      stop("NAs in fitness matrix")

    locsLetters <- gregexpr("f\\(([aA-zZ]?,? ?)*\\)", x[, ncol(x)])
    locsNumbers <- gregexpr("f\\(([0-9]?,? ?)*\\)", x[, ncol(x)])
    genesVectorLetters <-
      unique(unlist(strsplit(unlist(sapply(regmatches(x[, ncol(x)], locsLetters),
                                                               function(x){sub(".*\\((.*)\\).*",
                                                                               "\\1",
                                                                               x)})),
                                                 ", ")))
    genesVectorNumbers <-
      unique(unlist(strsplit(unlist(sapply(regmatches(x[, ncol(x)], locsNumbers),
                                                               function(x){sub(".*\\((.*)\\).*"
                                                                               ,"\\1",
                                                                               x)})),
                                                 ", ")))

    if(all(genesVectorLetters %in% c(LETTERS,  letters, ""))) {
      genotypesBy = "letters"
    }else if(is.numeric(locsNumbers)) {
      genotypesBy = "numbers"
    }else{
      stop("All genotypes must be letters or integers separated by ', '.")
    }

    if(genotypesBy == "letters") {
      locs <- gregexpr("f\\(([aA-zZ],? ?)*\\)", x[, ncol(x)])
      lettersGenesVector <-
        sort(unlist(strsplit(gsub(".*\\((.*)\\).*",
                                  "\\1",
                                  unlist(regmatches(x[, ncol(x)], locs))),
                             ", ")))
      lastLetter <- tolower(lettersGenesVector[length(lettersGenesVector)])
      n_max <- match(lastLetter, letters)

      if(n_max > ncol(x) - 1)
        stop("Gene's letters in fitness column must match with columns number")

      fitnessColumn <- sapply(x[, ncol(x)], function(x) replaceWithNumbers(x))
      x[, ncol(x)] <- fitnessColumn
    }
    return(x)
  }
}

## Deprecated after flfast
## to_genotFitness_std is faster and has better error checking
## and is very similar and does not use
## the genot_fitness_to_epistasis, which is not reasonable anymore.

## from_genotype_fitness <- function(x) {
##     ## Would break with output from allFitnessEffects and
##     ## output from allGenotypeAndMut

##     ## For the very special and weird case of
##     ## a matrix but only a single gene so with a 0 and 1
##     ## No, this is a silly and meaningless case.
##     ## if( ( ncol(x) == 2 ) && (nrow(x) == 1) && (x[1, 1] == 1) ) {

##     ## } else  blabla:

##     if(! (inherits(x, "matrix") || inherits(x, "data.frame")) )
##         stop("Input must inherit from matrix or data.frame.")

##     ## if((ncol(x) > 2) && !(inherits(x, "matrix"))
##     ##     stop(paste0("Genotype fitness input either two-column data frame",
##     ##          " or a numeric matrix with > 2 columns."))
##     ## if( (ncol(x) > 2) && (nrow(x) == 1) )
##     ##     stop(paste0("It looks like you have a matrix for a single genotype",
##     ##                 " of a single gene. For this degenerate cases use",
##     ##                 " a data frame specification."))

##     if(ncol(x) > 2) {
##         if(inherits(x, "matrix")) {
##             if(!is.numeric(x))
##                 stop("A genotype fitness matrix/data.frame must be numeric.")
##         } else if(inherits(x, "data.frame")) {
##             if(!all(unlist(lapply(x, is.numeric))))
##                 stop("A genotype fitness matrix/data.frame must be numeric.")
##         }

##         ## We are expecting here a matrix of 0/1 where columns are genes
##         ## except for the last column, that is Fitness
##         ## Of course, can ONLY work with epistastis, NOT order
##         return(genot_fitness_to_epistasis(x))
##     } else {
##         if(!inherits(x, "data.frame"))
##             stop("genotFitness: if two-column must be data frame")
##         ## Make sure no factors
##         if(is.factor(x[, 1])) x[, 1] <- as.character(x[, 1])
##         ## Make sure no numbers
##         if(any(is.numeric(x[, 1])))
##             stop(paste0("genotFitness: first column of data frame is numeric.",
##                         " Ambiguous and suggests possible error. If sure,",
##                         " enter that column as character"))

##         omarker <- any(grepl(">", x[, 1], fixed = TRUE))
##         emarker <- any(grepl(",", x[, 1], fixed = TRUE))
##         nogoodepi <- any(grepl(":", x[, 1], fixed = TRUE))
##         ## if(omarker && emarker) stop("Specify only epistasis or order, not both.")
##         if(nogoodepi && emarker) stop("Specify the genotypes separated by a ',', not ':'.")
##         if(nogoodepi && !emarker) stop("Specify the genotypes separated by a ',', not ':'.")
##         ## if(nogoodepi && omarker) stop("If you want order, use '>' and if epistasis ','.")
##         ## if(!omarker && !emarker) stop("You specified neither epistasis nor order")
##         if(omarker) {
##             ## do something. To be completed
##             stop("This code not yet ready")
##             ## You can pass to allFitnessEffects genotype -> fitness mappings that
##             ## involve epistasis and order. But they must have different
##             ## genes. Otherwise, it is not manageable.
##         }
##         if( emarker || ( (!omarker) && (!emarker) && (!nogoodepi)) ) {
##             ## the second case above corresponds to passing just single letter genotypes
##             ## as there is not a single marker
##             x <- x[, c(1, 2), drop = FALSE]
##             if(!all(colnames(x) == c("Genotype", "Fitness"))) {
##                 message("Column names of object not Genotype and Fitness.",
##                         " Renaming them assuming that is what you wanted")
##                 colnames(x) <- c("Genotype", "Fitness")
##             }
##             if((!omarker) && (!emarker) && (!nogoodepi)) {
##                 message("All single-gene genotypes as input to from_genotype_fitness")
##             }
##             ## Yes, we need to do this to  scale the fitness and put the "-"
##             return(genot_fitness_to_epistasis(allGenotypes_to_matrix(x)))
##         }
##     }
## }





## No longer used for real
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
    x <- x[, -ncol(x), drop = FALSE]
    wt <- which(rowSums(x) == 0)
    fwt <- 1
    if(length(wt) == 1)
        fwt <- f[wt]
    ## No longer being used when we pass fitness landscapse: flfast
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

    if(is.factor(x[, 1])) {
        warning("First column of genotype fitness is a factor. ",
                "Converting to character.")
        x[, 1] <- as.character(x[, 1])
    }
    ## A WT can be specified with string "WT"
    anywt <- which(x[, 1] == "WT")
    if(length(anywt) > 1) stop("More than 1 WT")
    if(length(anywt) == 1) {
        fwt <- x[anywt, 2]
        x <- x[-anywt, ]
        ## Trivial case of passing just a WT?
    } else {
        warning("No WT genotype. Setting its fitness to 1.")
        fwt <- 1
    }
    splitted_genots <- lapply(x$Genotype,
                             function(z) nice.vector.eo(z, ","))

    all_genes <- sort(unique(unlist(splitted_genots)))

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
    ## m <- data.frame(m)
    rs <- rowSums(m[, -ncol(m), drop = FALSE])
    m <- m[order(rs), , drop = FALSE]
    ## m <- m[do.call(order, as.list(cbind(rs, m[, -ncol(m)]))), ]
    return(m)
}


## Magellan_stats and Magellan_draw cannot be tested
## routinely, as they depend on external software
Magellan_stats <- function(x, max_num_genotypes = 2000,
                           verbose = FALSE,
                           use_log = TRUE,
                           short = TRUE,
                           fl_statistics = "fl_statistics",
                           replace_missing = FALSE) { # nocov start
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

    if(replace_missing) {
        zarg <- "-z"
    } else {
        zarg <- NULL
    }

    to_Magellan(x, fn, max_num_genotypes = max_num_genotypes)
    ## newer versions of Magellan provide some extra values to standard output
    call_M <- system2(fl_statistics,
                      args = paste(fn, shortarg, logarg, zarg, "-o", fnret),
                      stdout = NULL)
    if(short) {
        ## tmp <- as.vector(read.table(fnret, skip = 1, header = TRUE)[-1])

        tmp <- as.vector(read.table(fnret, skip = 1, header = TRUE)[c(-1)])
        ## Make names more explicit, but check we have what we think we have
        ## New versions of Magellan produce different output apparently of variable length
        stopifnot(length(tmp) >= 23) ## 23) ## variable length
        stopifnot(identical(names(tmp)[1:13], ## only some
                            c("ngeno", "npeaks", "nsinks", "gamma", "gamma.", "r.s",
                              "nchains", "nsteps", "nori", "depth", "magn", "sign",
                              "rsign"))) ## , "w.1.", "w.2.", "w.3..", "mode_w", "s.1.",
        ## "s.2.", "s.3..", "mode_s", "pairs_s", "outD_v")))
        if(length(tmp) >= 24) ## the new version
            stopifnot(identical(names(tmp)[c(20, 24)],
                                c("steps_m", "mProbOpt_0")))
        ## steps_m: the mean number of steps over the entire landscape to
        ## reach the global optimum
        ## mProbOpt_0: The mean probability to
        ## reach that optimum (again averaged over the entire landscape).
        ## but if there are multiple optima, there are many of these
        names(tmp)[1:13] <- c("n_genotypes", "n_peaks", "n_sinks", "gamma", "gamma_star",
                        "r/s","n_chains", "n_steps", "n_origins", "max_depth",
                        "epist_magn", "epist_sign", "epist_rsign")## ,
                        ## "walsh_coef_1", "walsh_coef_2", "walsh_coef_3", "mode_walsh",
                        ## "synerg_coef_1", "synerg_coef_2", "synerg_coef_3", "mode_synerg",
        ## "std_dev_pairs", "var_outdegree")
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
