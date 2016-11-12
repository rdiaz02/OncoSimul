## allready added properly

library(slam)

m1 <- simple_triplet_zero_matrix(100, mode = "double")

m1[1:20, 20:40] <- 1




wrm <- which(apply(m1[, -1, drop = FALSE], 2,
                   function(y) {all(na.omit(y) < 0.5)})) + 1

m2 <- m1[-wrm, wrm]



## If I only need to count accessible, I do not need this function in
## full. I can only examine those genotypes that are already not excluded,
## so set of "rs" is smaller each time

## I can also pre-compute the rs[i] + 1

## slam offers a direct replacement here?
genot_to_adj_mat <- function(x) {
    ## FIXME this can take about 23% of the time of the ggplot call.
    ## But them, we are quickly constructing a 2000*2000 matrix
    ## Given a genotype matrix, as given by allGenotypes_to_matrix, produce a
    ## directed adjacency matrix between genotypes connected (i.e., those
    ## that differ by gaining one mutation) with value being the
    ## difference in fitness between destination and origin

    ## Make sure sorted, so ancestors always before descendants
    rs0 <- rowSums(x[, -ncol(x)])
    x <- x[order(rs0), ]
    rm(rs0)
    
    y <- x[, -ncol(x)]
    f <- x[, ncol(x)]
    rs <- rowSums(y)

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
    return(adm)
}
