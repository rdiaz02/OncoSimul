
checkRegimeClonalInterference <- function(N, mu, numDrivers, s) {
    ## From Darayian and Shraiman, 2014, p. 914, where they cite other
    ## authors, this is when we are in a regime of clonal interference
    return( (N * mu * numDrivers) > (1/log(N * s)) )
}

## FIXME: prepare a plot from that function? range of N?
