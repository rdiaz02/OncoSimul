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
    ## But often, when we use allFitnessEffects, we
    ## obtain the fitness landscape as genotype_fitness_matrix
    ## FIXME: we could use that fact here.
    ## or maybe in to_Fitness_Matrix
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



to_Fitness_Matrix <- function(x, max_num_genotypes) {
    ## A general converter. Ready to be used by plotFitnessLandscape and
    ## Magellan exporter.
    
    ## Very bad, but there is no other way
    
    
    ## FIXME: really, some of this is inefficient. Very. Fix it.
    if( (inherits(x, "genotype_fitness_matrix")) ||
        ( (is.matrix(x) || is.data.frame(x)) && (ncol(x) > 2) ) ) {
        
        if("Fitness" %in% colnames(x)) {
            afe <- evalAllGenotypes(allFitnessEffects(
                genotFitness = x, frequencyDependentFitness = FALSE
                ##, epistasis = from_genotype_fitness(x)
            ),
            order = FALSE, addwt = TRUE, max = max_num_genotypes)
        } else {
            afe <- evalAllGenotypes(allFitnessEffects(
                genotFitness = x
                ##, epistasis = from_genotype_fitness(x)
            ),
            order = FALSE, addwt = TRUE, max = max_num_genotypes)
        }
        ## Why this? We go back and forth twice. We need both things. We
        ## could construct the afe below by appropriately pasting the
        ## columns names
        ## if( (is.null(colnames(x))) || any(grepl("^$", colnames(x))))
        ##    stop("Matrix x must have column names")

        ## This could use fmatrix_to_afe, above!!!
        ## Major change as of flfast: no longer using from_genotype_fitness
        

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
            if ("Death" %in% colnames(x)) {
                x <- rbind(data.frame(Genotype = "WT",
                                      Birth = 1, Death = 1,
                                      stringsAsFactors = FALSE),
                           x)
            } else {
                
                if("Fitness" %in% colnames(x)) {
                    x <- rbind(data.frame(Genotype = "WT",
                                          Fitness = 1,
                                          stringsAsFactors = FALSE),
                               x)
                } else {
                    x <- rbind(data.frame(Genotype = "WT",
                                          Birth = 1,
                                          stringsAsFactors = FALSE),
                               x)
                }
            }
            
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
## rfitness, do nothing. Well, it actually does something.

to_genotFitness_std <- function(x,
                                frequencyDependentBirth = FALSE,
                                frequencyDependentDeath = FALSE,
                                frequencyDependentFitness = NULL,
                                frequencyType = NA,
                                deathSpec = FALSE,
                                simplify = TRUE,
                                min_filter_birth_death = 1e-9,
                                sort_gene_names = TRUE) {
    ## Would break with output from allFitnessEffects and
    ## output from allGenotypeAndMut

    if(! (inherits(x, "matrix") || inherits(x, "data.frame")) )
        stop("Input must inherit from matrix or data.frame.")
    
    if (frequencyDependentDeath && !deathSpec) {
        deathSpec = TRUE
        warning("Assumming death is being specified. Setting deathSpec to TRUE.")
    }
    # Has to be 0's and 1's specification without freq_dep_birth and no death
    # or without freq_dep_birth and without freq_dep_death.
    
    if((ncol(x) > 2 && !deathSpec) || ncol(x) > 3) {
        if (!frequencyDependentBirth && !frequencyDependentDeath){
            if(inherits(x, "matrix")) {
                if(!is.numeric(x))
                    stop("A genotype fitness matrix/data.frame must be numeric.")
            } else if(inherits(x, "data.frame")) {
                if(!all(unlist(lapply(x, is.numeric))))
                    stop("A genotype fitness matrix/data.frame must be numeric.")
            }
        }else{
            if(!inherits(x, "data.frame"))
                stop("Input must inherit from data.frame.")
            if(ncol(x) == 0){
                stop("You have an empty data.frame")
            }
            if (frequencyDependentBirth && frequencyDependentDeath) {
                if(!all(unlist(lapply(x[,-c((ncol(x)-1):ncol(x))], is.numeric)))){
                    stop("All columns except from the last two must be numeric.")
                }
                
                if(is.factor(x[, ncol(x)])) {
                    warning("Last column of genotype birth-death is a factor. ",
                            "Converting to character.")
                    x[, ncol(x)] <- as.character(x[, ncol(x)])
                }
                
                if(is.factor(x[, ncol(x)-1])) {
                    warning("Second last column of genotype birth-death is a factor. ",
                            "Converting to character.")
                    x[, ncol(x)-1] <- as.character(x[, ncol(x)-1])
                }
                if(!all(unlist(lapply(x[, ncol(x)], is.character)))){
                    stop("All elements in last column must be character.")
                }
                
                if(!all(unlist(lapply(x[, ncol(x)-1], is.character)))){
                    stop("All elements in second last column must be character.")
                }
            }
            
            else if(frequencyDependentBirth && deathSpec) {
                if(!all(unlist(lapply(x[,-(ncol(x)-1)], is.numeric)))){
                    stop("All columns except from the second last must be numeric.")
                }
                if(is.factor(x[, ncol(x)-1])) {
                    warning("Second last column of genotype birth-death is a factor. ",
                            "Converting to character.")
                    x[, ncol(x)-1] <- as.character(x[, ncol(x)-1])
                }
                if(!all(unlist(lapply(x[, ncol(x)-1], is.character)))){
                    stop("All elements in second last column must be character.")
                }
            }
            
            else if(frequencyDependentDeath || (frequencyDependentBirth && !deathSpec)) {
                if(!all(unlist(lapply(x[,-ncol(x)], is.numeric)))){
                    stop("All columns except from the last one must be numeric.")
                }
                if(is.factor(x[, ncol(x)])) {
                    warning("Last column of genotype birth-death is a factor. ",
                            "Converting to character.")
                    x[, ncol(x)] <- as.character(x[, ncol(x)])
                }
                if(!all(unlist(lapply(x[, ncol(x)], is.character)))){
                    stop("All elements in last column must be character.")
                }
            }
        }

        ## We are expecting here a matrix of 0/1 where columns are genes
        ## except for the last column, that is Fitness
        ## Of course, can ONLY work with epistastis, NOT order
        ## return(genot_fitness_to_epistasis(x))
        if(any(duplicated(colnames(x))))
            stop("duplicated column names")
        
        if(deathSpec) {
            cnfl <- which(colnames(x)[-c((ncol(x)-1):ncol(x))] == "")
        }
        else {
            cnfl <- which(colnames(x)[-ncol(x)] == "")
        }
        
        if(length(cnfl)) {
            freeletter <- setdiff(LETTERS, colnames(x))[1]
            if(length(freeletter) == 0) stop("Renaiming failed")
            warning("One column named ''. Renaming to ", freeletter)
            colnames(x)[cnfl] <- freeletter
        }
        if(!is.null(colnames(x)) && sort_gene_names) {
            ncx <- ncol(x)
            if(deathSpec) {
                cnx <- colnames(x)[-c((ncx-1):ncx)]
            }
            else {
                cnx <- colnames(x)[-ncx]
            }
            
            ocnx <- gtools::mixedorder(cnx)
            if(!(identical(cnx[ocnx], cnx))) {
                message("Sorting gene column names alphabetically")
                
                if(!is.null(frequencyDependentFitness)) {
                    x <- cbind(x[, ocnx, drop = FALSE], Fitness = x[, (ncx)])
                } else {
                    x <- cbind(x[, ocnx, drop = FALSE], Birth = x[, (ncx)])
                }
                
            }
        }

        if(is.null(colnames(x))) {
            if(deathSpec) {
                ncx <- (ncol(x) - 2)
            }
            else {
                ncx <- (ncol(x) - 1)
            }
            
            message("No column names: assigning gene names from LETTERS")
            if(ncx > length(LETTERS))
                stop("More genes than LETTERS; please give gene names",
                     " as you see fit.")
            if(deathSpec) {
                colnames(x) <- c(LETTERS[1:ncx], "Birth", "Death")
            }
            
            else {
                
                if (!is.null(frequencyDependentFitness)) {
                    colnames(x) <- c(LETTERS[1:ncx], "Fitness")
                } else {
                    colnames(x) <- c(LETTERS[1:ncx], "Birth")
                }
                    
            }
            
        }
        
        if(deathSpec) {
            if(!all(as.matrix(x[, -c((ncol(x)-1):ncol(x))]) %in% c(0, 1) ))
                stop("First ncol - 2 entries not in {0, 1}.")
        }
        else {
            if(!all(as.matrix(x[, -ncol(x)]) %in% c(0, 1) ))
                stop("First ncol - 1 entries not in {0, 1}.")
        }

    } else {

        if(!inherits(x, "data.frame"))
            stop("genotFitness: if genotype is specified, it must be data frame")
        if(ncol(x) == 0){
            stop("You have an empty data.frame")
        }
        ## Make sure no factors
        if(is.factor(x[, 1])) {
            warning("First column of genotype birth-death is a factor. ",
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
        if(nogoodepi && emarker) stop("Specify the genotypes separated by a ',', not ':'.")
        if(nogoodepi && !emarker) stop("Specify the genotypes separated by a ',', not ':'.")
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
            if(deathSpec) {
                x <- x[, c(1, 2, 3), drop = FALSE]
                if(!all(colnames(x) == c("Genotype", "Birth", "Death"))) {
                    message("Column names of object not Genotype, Birth and Death.",
                            " Renaming them assuming that is what you wanted")
                    colnames(x) <- c("Genotype", "Birth", "Death")
                }
            }
            else {
                x <- x[, c(1, 2), drop = FALSE]
                if (!is.null(frequencyDependentFitness)) {
                    if(!all(colnames(x) == c("Genotype", "Fitness"))) {
                        message("Column names of object not Genotype and Birth",
                                " Renaming them assuming that is what you wanted")
                        colnames(x) <- c("Genotype", "Fitness")
                    }
                    
                } else {
                    if(!all(colnames(x) == c("Genotype", "Birth"))) {
                        message("Column names of object not Genotype and Birth",
                                " Renaming them assuming that is what you wanted")
                        colnames(x) <- c("Genotype", "Birth")
                    }
                }
                
            }
            
            
            if((!omarker) && (!emarker) && (!nogoodepi)) {
                message("All single-gene genotypes as input to to_genotFitness_std")
            }
            ## Yes, we need to do this to  scale the fitness and put the "-"
            if(frequencyDependentBirth || frequencyDependentDeath){
                anywt <- which(x[, 1] == "WT")
                if (length(anywt) > 1){
                    stop("WT should not appear more than once in birth-death specification")
                }
                
                if (frequencyDependentBirth) {
                    if(is.factor(x[, 2])) {
                        warning("Second column of genotype birth-death is a factor. ",
                                "Converting to character.")
                        x[, 2] <- as.character(x[, 2])
                    }
                }
                if (frequencyDependentDeath) {
                    if(is.factor(x[, 3])) {
                        warning("Third column of genotype birth-death is a factor. ",
                                "Converting to character.")
                        x[, 3] <- as.character(x[, 3])
                    }
                }
                
            }

            x <- allGenotypes_to_matrix(x, frequencyDependentBirth,
                                        frequencyDependentDeath, 
                                        frequencyDependentFitness, deathSpec)
        }
    }
    ## And, yes, scale all births and deaths by that of the WT

    if (!frequencyDependentBirth && !frequencyDependentDeath){
        if (deathSpec) {
            whichroot <- which(rowSums(x[, -c((ncol(x)-1):ncol(x)), drop = FALSE]) == 0)
        }
        else {
            whichroot <- which(rowSums(x[, -ncol(x), drop = FALSE]) == 0)
        }
        
        if(length(whichroot) == 0) {
            if(deathSpec) {
                warning("No wildtype in the fitness landscape!!! Adding it with birth and death 1.")
                x <- rbind(c(rep(0, ncol(x) - 2), 1, 1), x)
            }
            else {
                warning("No wildtype in the fitness landscape!!! Adding it with birth 1.")
                x <- rbind(c(rep(0, ncol(x) - 1), 1), x)
            }
            
        } else {
            if(x[whichroot, ncol(x)] != 1) {
                if(deathSpec) {
                    warning("Death of wildtype != 1.",
                            " Dividing all deaths by death of wildtype.")
                }
                else {
                    warning("Birth of wildtype != 1.",
                            " Dividing all births by birth of wildtype.")
                }
                
                vwt <- x[whichroot, ncol(x)]
                x[, ncol(x)] <- x[, ncol(x)]/vwt
            }
            
            if (deathSpec && x[whichroot, ncol(x)-1] != 1) {
                warning("Birth of wildtype != 1.",
                        " Dividing all births by birth of wildtype.")
                
                vwt <- x[whichroot, ncol(x)-1]
                x[, ncol(x)-1] <- x[, ncol(x)-1]/vwt
            }
        }
        
        if(!is.null(frequencyDependentFitness))
            colnames(x)[which(colnames(x) == "Birth")] <- "Fitness"
    }

    if(any(is.na(x)))
        stop("NAs in fitness matrix")
    
    if(!frequencyDependentBirth && !frequencyDependentDeath) {
        if(is.data.frame(x)) 
            x <- as.matrix(x)
        stopifnot(inherits(x, "matrix"))
    }
    
    if(simplify) {
        if ((!frequencyDependentBirth && !deathSpec) || (!frequencyDependentDeath && deathSpec)) {
            x <- x[x[, ncol(x)] > min_filter_birth_death, , drop = FALSE]
        }
        
        if (!frequencyDependentBirth && deathSpec) {
            x <- x[x[, ncol(x)-1] > min_filter_birth_death, , drop = FALSE]
        }
        
    }
    
    if (!frequencyDependentBirth && !frequencyDependentDeath)
        class(x) <- c("matrix", "genotype_birth_death_matrix")
    
    if(frequencyDependentBirth) { ## frequency-dependent fitness
        
        if(frequencyType == "auto"){
            if (deathSpec) {
                ch <- paste(as.character(x[, ncol(x)-1]), collapse = "")
            }
            else {
                ch <- paste(as.character(x[, ncol(x)]), collapse = "")
            }
            
            if( grepl("f_", ch, fixed = TRUE) ){
                frequencyType = "rel"
                pattern <- stringr::regex("f_(\\d*_*)*")
                
            } else if ( grepl("n_", ch, fixed = TRUE) ){
                frequencyType = "abs"
                pattern <- stringr::regex("n_(\\d*_*)*")
                
            } else { stop("No pattern found when frequencyType set to 'auto'") }
            
        } else if(frequencyType == "abs"){
            pattern <- stringr::regex("n_(\\d*_*)*")
        } else {
            pattern <- stringr::regex("f_(\\d*_*)*")
        }
        
        if(deathSpec) {
            regularExpressionVectorBirth <-
                unique(unlist(lapply(x[, ncol(x)-1],
                                     function(z) {stringr::str_extract_all(string = z,
                                                                           pattern = pattern)})))
            
            if((!all(regularExpressionVectorBirth %in% fVariablesN(ncol(x) - 2, frequencyType))) |
               !(length(intersect(regularExpressionVectorBirth,
                                  fVariablesN(ncol(x) - 2, frequencyType)) >= 1) )){
                stop("There are some errors in birth column")
            }
        }
        else {
            regularExpressionVectorBirth <-
                unique(unlist(lapply(x[, ncol(x)],
                                     function(z) {stringr::str_extract_all(string = z,
                                                                           pattern = pattern)})))
            if((!all(regularExpressionVectorBirth %in% fVariablesN(ncol(x) - 1, frequencyType))) |
               !(length(intersect(regularExpressionVectorBirth,
                                  fVariablesN(ncol(x) - 1, frequencyType)) >= 1) )){
                stop("There are some errors in birth column")
            }
        }
    }
    
    if(frequencyDependentDeath) { ## frequency-dependent fitness
        
        if(frequencyType == "auto"){
            ch <- paste(as.character(x[, ncol(x)]), collapse = "")
            
            if( grepl("f_", ch, fixed = TRUE) ){
                frequencyType = "rel"
                pattern <- stringr::regex("f_(\\d*_*)*")
                
            } else if ( grepl("n_", ch, fixed = TRUE) ){
                frequencyType = "abs"
                pattern <- stringr::regex("n_(\\d*_*)*")
                
            } else { stop("No pattern found when frequencyTypeDeath set to 'auto'") }
            
        } else if(frequencyType == "abs"){
            pattern <- stringr::regex("n_(\\d*_*)*")
        } else {
            pattern <- stringr::regex("f_(\\d*_*)*")
        }
        
        regularExpressionVectorDeath <-
            unique(unlist(lapply(x[, ncol(x)],
                                 function(z) {stringr::str_extract_all(string = z,
                                                                       pattern = pattern)})))
        if((!all(regularExpressionVectorDeath %in% fVariablesN(ncol(x) - 1, frequencyType))) |
           !(length(intersect(regularExpressionVectorDeath,
                              fVariablesN(ncol(x) - 1, frequencyType)) >= 1) )){
            stop("There are some errors in death column")
        }
    }
    return(x)
}





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




allGenotypes_to_matrix <- function(x,
                                   frequencyDependentBirth = FALSE,
                                   frequencyDependentDeath = FALSE,
                                   frequencyDependentFitness = NULL,
                                   deathSpec = FALSE) {
    ## Makes no sense to allow passing order: the matrix would have
    ## repeated rows. A > B and B > A both have exactly A and B

    ## Take output of evalAllGenotypes or identical data frame and return
    ## a matrix with 0/1 in a column for each gene and a final column of
    ## Fitness
    
    if(!is.null(frequencyDependentFitness))
        frequencyDependentBirth <- frequencyDependentFitness
    
    if (is.factor(x[, 1])) {
        warning(
            "First column of genotype birth-death is a factor. ",
            "Converting to character."
        )
        x[, 1] <- as.character(x[, 1])
    }
    if (frequencyDependentBirth) {
        if (is.factor(x[, 2])) {
            warning(
                "Second column of genotype birth-death is a factor. ",
                "Converting to character."
            )
            x[, 2] <- as.character(x[, 2])
        }
    }
    
    if (frequencyDependentDeath) {
        if (is.factor(x[, 3])) {
            warning(
                "Third column of genotype birth-death is a factor. ",
                "Converting to character."
            )
            x[, 3] <- as.character(x[, 3])
        }
    }

    ## A WT can be specified with string "WT"
    anywt <- which(x[, 1] == "WT")
    if (length(anywt) > 1) stop("More than 1 WT")
    if (length(anywt) == 1) {
        bwt <- x[anywt, 2]
        if (deathSpec) {
            dwt <- x[anywt, 3]
        }
        
        x <- x[-anywt, ]
        ## Trivial case of passing just a WT?
    } else {
        if (!frequencyDependentBirth && !frequencyDependentDeath) {
            bwt <- 1
            if(deathSpec) {
                dwt <- 1
                warning("No WT genotype. Setting its birth and death to 1.")
            }
            else {
                warning("No WT genotype. Setting its birth to 1.")
            }
            
            
        } else {
            bwt <- NA
            
            if(deathSpec) 
                dwt <- NA
            ##   message("No WT genotype in FDF: setting it to 0.")
        }
    }
    splitted_genots <- lapply(
        x$Genotype,
        function(z) nice.vector.eo(z, ",")
    )

    all_genes <- sort(unique(unlist(splitted_genots)))
    if(length(all_genes) < 2) stop(paste("There must be at least two genes (loci)",
                                         "in the fitness effects.",
                                         "If you only care about a case with",
                                         "a single one (or none) enter gene(s)",
                                         "with a fitness effect of zero.",
                                         "For freq.dep.fitness, create another ",
                                         "genotype that always has fitness zero."))
    m <- matrix(0, nrow = length(splitted_genots), ncol = length(all_genes))
    the_match <- lapply(
        splitted_genots,
        function(z) match(z, all_genes)
    )
    ## A lot simpler with a loop
    for (i in 1:nrow(m)) {
        m[i, the_match[[i]]] <- 1
    }
    if(deathSpec) {
        m <- cbind(m, x[, 2], x[, 3])
        colnames(m) <- c(all_genes, "Birth", "Death")
    }
    else {
        m <- cbind(m, x[, 2])
        
        if (!is.null(frequencyDependentFitness)) {
            colnames(m) <- c(all_genes, "Fitness")
        } else {
            colnames(m) <- c(all_genes, "Birth")
        }
    }
    
    
    if(!is.na(bwt))
        if(deathSpec) {
            m <- rbind(c(rep(0, length(all_genes)), bwt, dwt), m)
        }
        else{
            m <- rbind(c(rep(0, length(all_genes)), bwt), m)
        }
        

    if (frequencyDependentBirth || frequencyDependentDeath) {
        m <- as.data.frame(m)
        m[, 1:length(all_genes)] <- apply(
            m[, 1:length(all_genes), drop = FALSE],
            2,
            as.numeric
        )
        if(frequencyDependentBirth) {
            m[, length(all_genes)+1] <- as.character(m[, length(all_genes)+1])
        }
        else {
            m[, length(all_genes)+1] <- as.numeric(m[, length(all_genes)+1])
        }
        
        if(frequencyDependentDeath) {
            m[, length(all_genes)+2] <- as.character(m[, length(all_genes)+2])
        }
        else if(deathSpec){
            m[, length(all_genes)+2] <- as.numeric(m[, length(all_genes)+2])
        }
        
    }
    ## Ensure sorted
    ## m <- data.frame(m)
    if(deathSpec) {
        rs <- rowSums(m[, -c((ncol(m)-1):ncol(m)), drop = FALSE])
    }
    else {
        rs <- rowSums(m[, -ncol(m), drop = FALSE])
    }
    m <- m[order(rs), , drop = FALSE]
    ## m <- m[do.call(order, as.list(cbind(rs, m[, -ncol(m)]))), ]
    return(m)
}




Magellan_stats <- function(x, max_num_genotypes = 2000,
                           verbose = FALSE,
                           use_log = FALSE,
                           short = TRUE,
                           replace_missing = FALSE) {
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
    call_M <- system2(fl_statistics_binary(),
                      args = paste(shortarg, logarg, zarg, "-o", fnret, fn),
                      stdout = NULL)
    if(short) {
        tmp <- unlist(read.table(fnret, skip = 1, header = TRUE)[c(-1)])
        ## ## Make names more explicit, but check we have what we think we have
        ## ## New versions of Magellan produce different output apparently of variable length
        ## stopifnot(length(tmp) >= 23) ## 23) ## variable length
        ## stopifnot(identical(names(tmp)[1:13], ## only some
        ##                     c("ngeno", "npeaks", "nsinks", "gamma", "gamma.", "r.s",
        ##                       "nchains", "nsteps", "nori", "depth", "magn", "sign",
        ##                       "rsign"))) ## , "w.1.", "w.2.", "w.3..", "mode_w", "s.1.",
        ## ## "s.2.", "s.3..", "mode_s", "pairs_s", "outD_v")))
        ## if(length(tmp) >= 24) ## the new version
        ##     stopifnot(identical(names(tmp)[c(20, 24)],
        ##                         c("steps_m", "mProbOpt_0")))
        ## ## steps_m: the mean number of steps over the entire landscape to
        ## ## reach the global optimum
        ## ## mProbOpt_0: The mean probability to
        ## ## reach that optimum (again averaged over the entire landscape).
        ## ## but if there are multiple optima, there are many of these
        ## names(tmp)[1:13] <- c("n_genotypes", "n_peaks", "n_sinks", "gamma", "gamma_star",
        ##                 "r/s","n_chains", "n_steps", "n_origins", "max_depth",
        ##                 "epist_magn", "epist_sign", "epist_rsign")## ,
        ##                 ## "walsh_coef_1", "walsh_coef_2", "walsh_coef_3", "mode_walsh",
        ##                 ## "synerg_coef_1", "synerg_coef_2", "synerg_coef_3", "mode_synerg",
        ## ## "std_dev_pairs", "var_outdegree")
    } else {
        message("Output file: ", fnret)
        tmp <- readLines(fnret)
    }
    return(tmp)
}



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
