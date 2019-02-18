## Examples of cPDetect, PDBaseline, etc.

cPDetect <- function(n2, p2, PDBaseline) {
    return (-log(1.0 - p2)/(n2 - PDBaseline));
}



## comparing with different initSizes fixing n2.

initSize <- 2000
cPDetect(n2 = 2  * initSize, p2 = 0.1, PDBaseline = 1.3 * initSize)
## 7.526e-05

initSize <- 5e4
cPDetect(n2 = 2  * initSize, p2 = 0.1, PDBaseline = 1.3 * initSize)
## 3.01e-06

initSize <- 1e6
cPDetect(n2 = 2  * initSize, p2 = 0.1, PDBaseline = 1.3 * initSize)
## 1.505e-07


## Now, the other way around. If we set the value at
## cPDetect = 7.526e-05, what n2 is that?

n_for_p2 <- function(p2, cPDetect, PDBaseline) {
    return( PDBaseline -  (log(1 - p2)/cPDetect)  )
}

