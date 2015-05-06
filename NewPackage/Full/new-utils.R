
regimeClonalInterference <- function(N, mu, numDrivers, s) {
    ## From Darayian and Shraiman, 2014, p. 914, where they cite other
    ## authors, this is when we are in a regime of clonal interference
    return( (N * mu * numDrivers) > (1/log(N * s)) )
}

## FIXME: prepare a plot from that function? range of N?


plotRegimeClonalInterference <- function(N, numDrivers, mu, s, npoints = 100) {
    ## N is given as a two-element vector
    ## numDrivers the number of drivers to examine
    ns <- seq(from = N[1], to = N[2], length.out = npoints)
    data <- sapply(numDrivers, function(x)
                   regimeClonalInterference(ns, mu, x, s))

    matplot(x = ns, y = data, type = "l", ylab = c("No", "Yes"),
            lty = seq_along(numDrivers), col = seq_along(numDrivers),
            axes = FALSE)
    box()
    axis(1)
    axis(2, at = c(0, 1), labels = c("No", "Yes"))
    legend(N[1], 1, numDrivers,
           lty = seq_along(numDrivers), col = seq_along(numDrivers))
}

## do it for McFL model. Remember that a relationship between N and number
## of drivers.
