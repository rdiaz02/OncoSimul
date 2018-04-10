## Examples of cPDetect, PDBaseline, etc.

cPDetect <- function(n2, p2, PDBaseline) {
    return (-log(1.0 - p2)/(n2 - PDBaseline));
}



## comparing with different initSizes fixing n2.
## Note how the cPDetect changes a lot. As we increase initSize,
## we delay the speed of detection.

initSize <- 2000
cPDetect(n2 = 2  * initSize, p2 = 0.1, PDBaseline = 1.3 * initSize)
## 7.526e-05

initSize <- 5e4
cPDetect(n2 = 2  * initSize, p2 = 0.1, PDBaseline = 1.3 * initSize)
## 3.01e-06

initSize <- 1e6
cPDetect(n2 = 2  * initSize, p2 = 0.1, PDBaseline = 1.3 * initSize)
## 1.505e-07


## original expression
dp1 <- function(n, B, C) {
    if(n < B) return(0)
    1 - exp(-C * ( (n - B) ))
}


## Now, the other way around. If we set the value at
## cPDetect, what n2 is that?

n_for_p2_dp1 <- function(p2, cPDetect, PDBaseline) {
    return( PDBaseline -  (log(1 - p2)/cPDetect)  )
}


## check
n_for_p2_dp1(p2 = 0.1, cPDetect = cPDetect(2000, 0.1, 2000 * 1.3), 2000 * 1.3)

n_for_p2_dp1(p2 = 0.1, cPDetect = cPDetect(2000, 0.1, 2000 * 1.3), 5e4 * 1.3)
## 64400

n_for_p2_dp1(p2 = 0.1, cPDetect = cPDetect(2000, 0.1, 2000 * 1.3), 1e6 * 1.3)
## 1299400


## So when using this mechanism, you probably want to fix cPDetect, not n2
## and p2, if you vary over a large range of initSizes???  That is tricky:
## the formula uses number of cells of increase over PDBasline, not
## proportion of increase over PDBaseline.



######
## Different process:

## pop size, Baseline, C (cpdetect)
dp2 <- function(n, B, C) {
    if(n < B) return(0)
    1 - exp(-C * ( (n - B)/B ))
}

## the C for a given prob detection
C_dp2 <- function(n, p, B) {
    (-log(1 - p)) * (B/(n-B)) 
}

C_dp2(4000, 0.1, 2000 * 1.3)
## 0.1957


## The N for a given C and p2

n_for_p2_dp2 <- function(p2, C, B) {
    B * ( 1 - ( log(1-p2)/C ) )
}


## check
dp2(4000, 2000 * 1.3, 0.2)
n_for_p2_dp2(0.1021, 0.2, 2000 * 1.3)




## compare the p2 at same n2 factor
is <- 2000
dp2(2 * is, 1.3 * is, 0.1957) ## 0.1

is <- 5e4
dp2(2 * is, 1.3 * is, 0.1957) ## 0.1

is <- 1e6
dp2(2 * is, 1.3 * is, 0.1957) ## 0.1

## and of course in all cases we have 0.1957
initSize <- 2000
C_dp2(2 * initSize, p = 0.1, B = 1.3 * initSize)

initSize <- 5e4
C_dp2(2 * initSize, p = 0.1, B = 1.3 * initSize)

initSize <- 1e6
C_dp2(2 * initSize, p = 0.1, B = 1.3 * initSize)






dp2(3000, 1600, .4)
dp1(3000, 1600, 1e-4)

op <- par(mfrow = c(2, 1))
initSize <- 2000
curve(dp2(x, B = 1.3 * initSize, C = 0.1957),
      from = 1.3 * initSize, to = 20 * initSize)
curve(dp1(x, B = 1.3 * initSize, C = 7.526e-05),
      from = 1.3 * initSize, to = 20 * initSize)

## basically identical
x <- seq(from = 2000, to = 2000 * 5, length.out = 200)
y_dp1 <- sapply(x, function(z) dp1(z, 2000 * 1.3, 7.526e-05))
y_dp2 <- sapply(x, function(z) dp2(z, 2000 * 1.3, 0.1957))
plot(y_dp2 ~ y_dp1)                



dp1(4000, 2000 * 1.3, 7.526e-05) ## 0.1










