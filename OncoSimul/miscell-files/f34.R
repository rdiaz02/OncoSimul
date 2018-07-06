f3 <- function(x, k, baseline = 2000) {
    if(x <= baseline) return(0)
    else return( 1 - exp(-k * (x - baseline)))
}

f4 <- function(x, k, baseline = 2000) {
    if(x <= baseline) return(0)
    else return( 1 - exp(-k * ((x - baseline)/baseline)))
}


## For R, vetorized



f4 <- function(x, k, baseline = 1000) {pmax(0, 1 - exp(-k*((x - baseline) / baseline)))}

f3 <- function(x, k, baseline = 1000) {pmax(0, 1 - exp(-k*(x - baseline)))}


##

f3auto <- function(x, k, n2 = NULL, p2 = NULL, baseline = 1000) {
    if(!is.null(n2) && !is.null(p2))
        k <- -log(1 - p2)/(n2 - baseline)
    pmax(0, 1 - exp(-k*(x - baseline)))
}


f4auto <- function(x, k, n2 = NULL, p2 = NULL, baseline = 1000) {
    if(!is.null(n2) && !is.null(p2))
        k <- -log(1 - p2)/((n2 - baseline)/baseline)
    pmax(0, 1 - exp(-k*((x - baseline) / baseline)))
}
