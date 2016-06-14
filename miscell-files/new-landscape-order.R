generateAllGenotypes_minimal <- function(fitnessEffects, max = 256) {
    ## This generates the strictly needed set of genotypes: it uses order
    ## when it should, and does not when it shouldn't. And adapts usage to
    ## the different parts of the fitnessEffects specification.

    ## We use a simplifying assumption that might give somewhat larger
    ## tables: if order is ever used, then we obtain all ordered genotypes
    ## for everything not in geneNoInt. And geneNoInt as unordered. And
    ## then we do the Cartesian product.
   
    tot_o <- function(n) {sum(sapply(seq.int(n),
                                     function(x) choose(n, x) * factorial(x)))}
    tot_no <- function(n) {2^n}

    order <- order_no_int <- no_order <- FALSE
    
    if(!is.null(fitnessEffects$orderEffects)) {
        order <- TRUE
        n_o <- nrow(fitnessEffects$geneModule) -1

        if(!is.null(fitnessEffects$noIntGenes)) {
            order_no_int <- TRUE
            n_no <- nrow(fitnessEffects$long.geneNoInt)
        }
        tnn <- tot_o(n_o) + tot_no(n_no)
    } else {
        no_order <- TRUE
        n_no <- nrow(fitnessEffects$geneModule) -1  + nrow(fitnessEffects$long.geneNoInt)
        tnn <- tot_no(n_no)
    }

    if(tnn > max) {
        m <- paste("There are", tnn, "genotypes.")
        m <- paste(m, "This is larger than max.")
        m <- paste(m, "Adjust max and rerun if you want")
        stop(m)
    }

    vido <- vidno <- vector("integer")
    if(order) {
        vido <- fitnessEffects$geneModule$GenNumID[-1] ## 0 not in there
        if(order_no_int)
            vidno <- fitnessEffects$long.geneNoInt$GeneNumID
    } else if(no_order) {
        vidno <- allNamedGenes(fitnessEffects)$GeneNumID
    }

    ## Yes, I keep thinking there's got to be a simpler way :-)
    ## Anyway, we obtain a list of the vectors
    list.of.vectors <- function(y) {
        ## there's got to be a simpler way
        lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}), recursive = FALSE),
               function(m) m[[1]])
    }
    
    ## Function to call with ordered part
    f1_o <- function(n) {
        lapply(seq.int(n), function(x) permutations(n = n, r = x, v = vido))}
    ## Ditto in non-ordered part
    f1_no <- function(n) {
        lapply(seq.int(n), function(x) combinations(n = n, r = x, v = vidno))}

    names <- c(fitnessEffects$geneModule$Gene[-1],
               fitnessEffects$long.geneNoInt$Gene)

    f_names <- function(ggnums) {
        unlist(lapply(lapply(ggnums, function(x) names[x]),
                      function(z)
                          paste(z,
                                collapse = if(order){" > "} else {", "} )))
    }
    
    if(order) {
        genotNums_o <- f1_o(n_o)
        genotNums_o <- list.of.vectors(genotNums_o)
        genotNames_o <- f_names(genotNums_o)
        if(order_no_int) {
            genotNums_no <- f1_no(n_no)
            genotNums_no <- list.of.vectors(genotNums_no)
            genotNames_no <- f_names(genotNums_no)
            }
    } else if(no_order) {
        genotNums_no <- f1_no(n_no)
        genotNums_no <- list.of.vectors(genotNums_no)
        genotNames_no <- f_names(genotNums_no)
    }

    
    

    

    
    ## ## With mutator, the ids of genes need not go from 1:n
    ## vid <- allNamedGenes(fitnessEffects)$GeneNumID
    ## if(order) {
    ##     f1 <- function(n) {
    ##         lapply(seq.int(n), function(x) permutations(n = n, r = x, v = vid))
    ##     }
    ## } else {
    ##     f1 <- function(n) {
    ##         lapply(seq.int(n), function(x) combinations(n = n, r = x, v = vid))}
        
    ## }
   
    ## genotNums <- f1(nn)
    ## list.of.vectors <- function(y) {
    ##     ## there's got to be a simpler way
    ##     lapply(unlist(lapply(y, function(x) {apply(x, 1, list)}), recursive = FALSE),
    ##            function(m) m[[1]])
    ## }
    ## genotNums <- list.of.vectors(genotNums)
    ## names <- c(fitnessEffects$geneModule$Gene[-1],
    ##            fitnessEffects$long.geneNoInt$Gene)
    
    ## genotNames <- unlist(lapply(lapply(genotNums, function(x) names[x]),
    ##                             function(z)
    ##                                 paste(z,
    ##                                       collapse = if(order){" > "} else {", "} )))
    ## This ain't the best, as we repeatedly read and convert
    ## fitnessEffects.  If this were slow, prepare C++ function that takes
    ## vectors of genotypes and uses same fitnessEffects.
    return(list(genotNums = genotNums,
                genotNames = genotNames
                ))
}
