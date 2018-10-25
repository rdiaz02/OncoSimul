## rfitness is a lot nicer. This is a simple implementation of the models
## as explained in Ferretti et al., "Measuring epistasis in fitness
## landscapes: The correlation of fitness effects of mutations", JTB, 2016


HoC <- function(g, v = 1) {
    ## v is variance of the log of fitness
    m <- generate_matrix_genotypes(g)
    return(cbind(m,
                 Fitness = rlnorm(nrow(m), meanlog = 0, sdlog = sqrt(v))))
}

mult <- function(g, mu = 1, v = 1) {
    m <- generate_matrix_genotypes(g)
    f_hoc <- rlnorm(nrow(m), meanlog = 0, sdlog = sqrt(v_hoc))
    s_mult <- rlnorm(g, meanlog = mu_mult, sdlog = sqrt(v_mult))
    f_mult <- m %*% s_mult
    m <- cbind(m, Fitness = f_hoc + f_mult)
    return(m)
}

RMF <- function(g, v_hoc = 1, mu_mult = .5, v_mult = .1) {
    m <- generate_matrix_genotypes(g)
    f_hoc <- rlnorm(nrow(m), meanlog = 0, sdlog = sqrt(v_hoc))
    s_mult <- rlnorm(g, meanlog = mu_mult, sdlog = sqrt(v_mult))
    f_mult <- m %*% s_mult
    m <- cbind(m, Fitness = f_hoc + f_mult)
    return(m)
}


## rfitness <- function(g, model = c("HoC", "lHoC", "additive", "RMF"),
##                      v_hoc = .1, mu_mult = .1, v_mult = .01) {
##     browser()
##     ## v_hoc is variance of the log of fitness
##     model <- match.arg(model)
##     m <- generate_matrix_genotypes(g)
##     if(model == "lHoC")
##         f <- rlnorm(nrow(m), meanlog = 0, sdlog = sqrt(v_hoc))
##     else if(model == "HoC")
##         f <- rnorm(nrow(m), mean = 0, sd = sqrt(v_hoc))
##     else if(model == "additive") {
##         s_mult <- rlnorm(g, meanlog = mu_mult, sdlog = sqrt(v_mult))
##         f <- m %*% s_mult
##     } else if(model == "RMF") {
##         s_mult <- rlnorm(g, meanlog = mu_mult, sdlog = sqrt(v_mult))
##         f_mult <- m %*% s_mult
##         f_hoc <- rlnorm(nrow(m), meanlog = 0, sdlog = sqrt(v_hoc))
##         f <- f_hoc + f_mult
##     }
##     m <- cbind(m, Fitness = f)
##     return(m)
## }

