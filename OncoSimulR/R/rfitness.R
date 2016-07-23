rfitness <- function(g, c= 0.5,
                     sd = 1,
                     reference = "random", ## "random", "max", or the vector,
                                     ## e.g., rep(1, g). If random, a
                                     ## random genotype is chosen as
                                     ## reference. If "max" this is rep(1, g)
                     scale = NULL, ## a two-element vector: min and max
                     wt_is_1 = TRUE, ## wt has fitness 1
                     log = FALSE, ## log only makes sense if all values >
                                 ## 0. scale with min > 0, and/or set
                                 ## wt_is_1 = TRUE
                     min_accessible_genotypes = 0,
                     accessible_th = 0) {
    ## Like Franke et al., 2011 and others of Krug. Very similar to Greene
    ## and Crona, 2014. And this allows moving from HoC to purely additive
    ## changing c and sd.

    ## FIXME future: do this with order too?
    ##    - do not generate a matrix of genotypes but a matrix of ordered genot.
    m <- generate_matrix_genotypes(g)
    done <- FALSE
    while(!done) {
        f_r <- rnorm(nrow(m), mean = 0, sd = sd)
        if(inherits(reference, "character") && length(reference) == 1) {
            if(reference == "random") {
                referenceI <- m[sample(nrow(m), 1), ]
            } else if(reference == "max") {
                referenceI <- rep(1, g)
            } 
        } else {
            referenceI <- reference
            }
        d_reference <- apply(m, 1, function(x) sum(abs(x - referenceI)))
        f_det <- -c * d_reference
        ## f_det <- rowSums(m) * slope/nrow(m) ## this is Greene and Krona
        fi <- f_r + f_det
        
        if(!is.null(scale)) {
            fi <- (fi - min(fi))/(max(fi) - min(fi))
            fi <- scale[1] + fi * (scale[2] - scale[1])
        }
        if(wt_is_1) {
            ## We need to shift to avoid ratios involving negative numbers and
            ## we need to avoid having any fitness at 0, specially the wt.  If
            ## any negative value, add the min, and shift by the magnitude of
            ## the min to avoid any at 0.

            ## If you use scale and wt_is_1, this will move the scale. It is
            ## not possible to obtain a linear transformation that will keep
            ## the min and max of the scale, and wt at 1.
            min_fi <- min(fi)
            if(min_fi < 0)
                fi <- fi + 2 * abs(min(fi))
            fi <- fi/fi[1]
        }
        
        if(log) {
            fi <- log(fi/fi[1]) + 1
        }
        m <- cbind(m, Fitness = fi)
        if(min_accessible_genotypes > 0) {
            num_accessible_genotypes <- count_accessible_g(m, accessible_th)
            if(num_accessible_genotypes >= min_accessible_genotypes) {
                done <- TRUE
                attributes(m) <- c(attributes(m),
                                   accessible_genotypes = num_accessible_genotypes,
                                   accessible_th = accessible_th)
            } else {
                ## Cannot start again with a fitness column
                m <- m[, -ncol(m), drop = FALSE]
            }
        } else {
            done <- TRUE
        }
    }
    class(m) <- c(class(m), "genotype_fitness_matrix")
    return(m)
}
