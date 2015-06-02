## keepTheseMany

##


pruneK <- function(x, keepTheseMany = 200) {
    ## x is the object
    tp <- x$pops.by.time

    kept <- nrow(tp)
    if(kept <= keepTheseMany) {
        return(x)
    }

    to.keep <- rev(round(seq(from = kept, to = 1, length.out = keepTheseMany)))

    tp <- tp[to.keep, ]
    remove.genotypes <- which(colSums(tp) == 0)

    x$pops.by.time <- tp[, remove.genotypes]
    x$Genotypes <- x$Genotypes[, remove.genotypes]
}
