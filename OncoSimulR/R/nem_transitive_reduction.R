###############################################
## <code from nem package>

## Title: (Dynamic) Nested Effects Models and Deterministic Effects
##         Propagation Networks to reconstruct phenotypic hierarchies
## Version: 2.60.0
## Author: Holger Froehlich, Florian Markowetz, Achim Tresch, Theresa
## Niederberger, Christian Bender, Matthias Maneck, Claudio Lottaz, Tim
## Beissbarth
## Maintainer: Holger Froehlich <frohlich at bit.uni-bonn.de>
## https://bioconductor.org/packages/release/bioc/html/nem.html
## License: GPL (>= 2)

nem_transitive.reduction <- function(g){
    if (!any(class(g)%in%c("matrix","graphNEL"))) stop("Input must be an adjacency matrix or graphNEL object")
    if(any(class(g) == "graphNEL")){
        g = as(g, "matrix")		
    }

    g = nem_transitive.closure(g, mat=TRUE) # bug fix: algorithm only works for transitively closed graphs!
    g = g - diag(diag(g))
    type = (g > 1)*1 - (g < 0)*1	
    for(y in 1:nrow(g)){
        for(x in 1:nrow(g)){
            if(g[x,y] != 0){
                for(j in 1:nrow(g)){
                    if((g[y,j] != 0) & sign(type[x,j])*sign(type[x,y])*sign(type[y,j]) != -1){ 
                        g[x,j] = 0
                    }
                }
            }
        }
    }
    g		
}



nem_transitive.closure <- function(g,mat=FALSE,loops=TRUE){

    if (!any(class(g)%in%c("graphNEL","matrix"))) stop("Input must be either graphNEL object or adjacency matrix")
    g <- as(g, "matrix")
    
                                        #-- adjacency matrix
                                        #     if (class(g)=="matrix"){
    n <- ncol(g)
    matExpIterativ <- function(x,pow,y=x,z=x,i=1) {
        while(i < pow) {
            z <- z %*% x
            y <- y+z
            i <- i+1
        }
        return(y)
    }
    
    h <- matExpIterativ(g,n)
    h <- (h>0)*1   
    dimnames(h) <- dimnames(g)
    if (!loops) diag(h) <- rep(0,n) else diag(h) <- rep(1,n)
    if (!mat) h <- as(h,"graphNEL")	
                                        #     }

    return(h)
}
