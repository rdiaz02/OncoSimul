## Miscell attempts. Most of this is not working.









## Stream does not always give you what you'd like conveying newer mutants.
## For example, using b11

## Testcode for log = y and log = x


### Why not stacked or stream plots with log transformed data:
## log(a + b) != log(a) + log(b)

## but stacked plots depend on additivity for the output to make
## sense. You can have an increase in total population but the stacked
## plot with log data would suggest you get a decrease.

## plot.stacked(1:2, log10(cbind(c(5, 1), c(5, 11))))
## plot.stacked(1:2, log10(cbind(c(6, 2), c(8, 14))))



###############

sa <- 0.1
sb <- -0.2
sab <- 0.25
sac <- -0.1
sbc <- 0.25

sv2 <- allFitnessEffects(epistasis = c("-A : B" = sb,
                                       "A : -B" = sa,
                                       "A : C" = sac,
                                       "A:B" = sab,
                                       "-A:B:C" = sbc),
                         geneToModule = c(
                             "Root" = "Root",
                             "A" = "a1, a2",
                             "B" = "b",
                             "C" = "c"))
evalAllGenotypes(sv2, order = FALSE, addwt = TRUE)


e1 <- oncoSimulIndiv(sv2, model = "McFL",
                     mu = 5e-6,
                     sampleEvery = 0.02,
                     keepEvery = 1,
                     initSize = 2000,
                     finalTime = 3000,
                     onlyCancer = FALSE)


########## colors, etc

### Notes for me




## for distinct colors, from gplots
rich.colors(5)
?rich.colors ## from gplots
plot(1:8, col = rich.colors(8), pch = 16, cex = 3)
## From the example on the colorspace vignette
plot(1:6, col = rainbow_hcl(6, c = 60, l = 75), pch = 16, cex = 3)

plot(1:6, col = hsv(seq(from = 0, to = 1, length.out = 10), 1, 1), pch = 16, cex = 3)



############### with streamgraph package







devtools::install_github("hrbrmstr/streamgraph")
library(streamgraph)

b4 <-  oncoSimulIndiv(mcf1,
                              model = "McFL", 
                              mu = 1e-7,
                              detectionSize = 1e8, 
                              detectionDrivers = 100,
                              sampleEvery = 0.02,
                              keepEvery = 10,
                              initSize = 2000,
                              finalTime = 2000,
                              onlyCancer = FALSE)
dd1 <- dwla(b4)

## dd4 <- dd1
## dd4 <- dd1[dd1$Drivers %in% c(1),  ]
## summary(dd4); nrow(dd4)
## length(unique(dd4$Genotype))

streamgraph(dd1, Genotype, Y, Time, scale = "continuous") %>% sg_fill_brewer("PuOr")

streamgraph(dd1, Genotype, Y, Time, scale = "continuous", offset = "zero")

## But probably hard to use a range of colors by drivers. And this ain't in CRAN yet.


dwl <- function(time, y, ndr, genotLabels) {
    nc <- ncol(y)
    nr <- nrow(y)
    y <- as.vector(y)
    df <- data.frame(Time = rep(time, nc),
                     Y = y,
                     Drivers = rep(ndr, rep(nr, nc)),
                     Genotype = rep(genotLabels, rep(nr, nc)))
    return(df)
}


dwlb <- function(x, min.col = .2, max.col = .7) {

    ## What if I do not care about drivers? All have same drivers
    ## What if I only want genotypes?
    if(!inherits(x, "oncosimul2")) {
        ndr <- colSums(x$Genotypes[1:x$NumDrivers, , drop = FALSE])
    } else {
        ndr <- colSums(x$Genotypes[x$Drivers, , drop = FALSE])
    }

    y <- x$pops.by.time[, 2:ncol(x$pops.by.time), drop = FALSE]

    ## oo <- order(ndr)
    ## y <- y[, oo, drop = FALSE]
    ## ndr <- ndr[oo]

    nc <- ncol(y)
    nr <- nrow(y)
    
    df <- data.frame(Time = rep(x$pops.by.time[, 1], nc),
                     Y = as.vector(y),
                     Drivers = rep(ndr, rep(nr, nc)), ## maybe not needed?
                     Genotype = rep(x$GenotypesLabels, rep(nr, nc)),
                     Cl = NA)

    ldr <- unique(ndr)

    for( i in ldr ) {
        ngenot <- length(unique(x$GenotypesLabels[ndr == i])) ## what if former class objects? deprecate
        vl <- seq(from = min.col, to = max.col, length.out = ngenot)
        colour <- rep(i + vl, rep(nr, ngenot))
        df$Cl[df$Drivers == i] <- colour
    }

    return(df)
}


plot.oncosimul2(b4)

dd2 <- dwlb(b4)


ggplot(dd2, aes(x = Time, y = Y, fill = Genotype)) + geom_area() +
    guides(fill = FALSE) 


ggplot(dd2, aes(x = Time, y = Y, fill = as.character(as.integer(Cl * 100000)))) + geom_area() +
    guides(fill = FALSE) 

## http://novyden.blogspot.com.es/2013/09/how-to-expand-color-palette-with-ggplot.html
colourCount <-  length(unique(dd2$Cl))
library(RColorBrewer)
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

## works, but not the legend I want

ggplot(dd2, aes(x = Time, y = Y, fill = as.character(as.integer(Cl * 100000)))) + geom_area() +
    scale_fill_manual(values = getPalette(colourCount))

## no legend, but I am not using breaks correctly
unix.time( {
    
    p <- ggplot(dd2, aes(x = Time, y = Y, fill = as.character(as.integer(Cl * 100000)))) + geom_area() +
        scale_fill_manual(values = getPalette(colourCount), breaks = levels(dd2$Drivers))
    p
}
)

## But I want the fill to move smoothly, like this!
## http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/

## The counterpart of this
ggplot(dd2, aes(x=Time, y = Y, colour=Cl)) + geom_point() + 
    scale_colour_gradientn(colours=rainbow(length(unique(dd2$Drivers))),
                           breaks = 0:15)

## Nope
ggplot(dd2, aes(x=Time, y = Y, colour=Cl)) + geom_area() + 
    scale_colour_gradientn(colours=rainbow(length(unique(dd2$Drivers))),
                           breaks = 0:15)

## Nope
ggplot(dd2, aes(x=Time, y = Y, fill = Cl)) + geom_area() + 
    scale_colour_gradientn(colours=rainbow(length(unique(dd2$Drivers))),
                           breaks = 0:15)

### I think this is hopeless:  https://groups.google.com/forum/#!topic/ggplot2/lihA-I8INpY

### or yes? http://stackoverflow.com/questions/27189453/shade-fill-or-color-area-under-density-curve-by-quantile

### we need a huge length.out to get cl3 to coincide with actual Cl
dd2$cl3 <- cut(dd2$Cl, seq(0, 15, length.out = 100000))

## these two take forever
ggplot(dd2, aes(x=Time, y = Y)) +
    geom_area(aes(x = Time, y = Y, group = cl3,
                  fill = cl3)) + guides(fill = FALSE)
ggplot(dd2, aes(x=Time, y = Y, fill = cl3)) +
    geom_area() + guides(fill = FALSE)


## this is faster (not fast)
ggplot(dd2, aes(x=Time, y = Y, fill = as.character(Cl))) + geom_area() + guides(fill = FALSE)




## nope, won't work
ggplot(dd2, aes(x=Time, y = Y, fill = Genotype)) + geom_area() + guides(fill = FALSE) +
    scale_fill_gradient2(midpoint = Cl)


 




## no legend
ggplot(dd2, aes(x = Time, y = Y, fill = factor(Cl))) + geom_area() +
    scale_fill_manual(values = getPalette(colourCount),
                      breaks = levels(as.factor(dd2$Drivers)),
                      )

tt <- seq(from = 0, to = 2000, by = 50)
dd3 <- dd2[dd2$Time %in% tt, ]


ggplot(dd3, aes(x = Time, y = Y, colour = as.character(Colour))) + geom_area() +
    scale_colour_gradient(colours = rainbow(length(unique(dd3$Drivers))))



ggplot(dd3, aes(x = Time, y = Y, fill = Genotype)) + geom_area() + guides(fill = FALSE) +
    scale_fill_gradientn(colours = rainbow(length(unique(dd3$Drivers))))





## ok
ggplot(dd3, aes(x = Time, y = Y, fill = Genotype)) + geom_area() + guides(fill = FALSE)

## this does not, as we hit the maximum for the palette
ggplot(dd3, aes(x = Time, y = Y, fill = as.character(as.integer(Cl * 100000)))) + geom_area() +
    scale_fill_brewer(breaks = levels(dd3$Drivers))


## ok
ggplot(dd3, aes(x = Time, y = Y, fill = Genotype)) + geom_area() + guides(fill = FALSE) +
    scale_fill_manual(values = rainbow(727))


##
ggplot(dd3, aes(x = Time, y = Y, group = Cl, fill = Cl)) + geom_area() + 
    scale_fill_manual(values = rainbow(727))




## ok, does something, but wrong.
ggplot(dd3, aes(x = Time, y = Y, fill = as.character(as.integer(Cl * 100000)),
                colour = Drivers)) + geom_area() + guides(fill = FALSE) +
    


## ok
ggplot(dd3, aes(x = Time, y = Y, colour = Cl)) + geom_point() +
    scale_colour_gradientn(colours = rainbow(length(unique(dd3$Drivers))))
ggplot(dd3, aes(x = Time, y = Y, colour = as.character(Cl))) + geom_area() +
    scale_colour_gradientn(colours = rainbow(length(unique(dd3$Drivers))))
ggplot(dd3, aes(x = Time, y = Y)) + geom_area(fill = dd3$Genotype)



## Simple use of ggplot, but slow 
ggplot(dd1, aes(x = Time, y = Y, fill = Genotype)) + geom_area() +
    guides(fill = FALSE)


myhsv <- function(major = 5, minor = 6) {
    h <- seq(from = 0, to = 1, length.out = major + 1)[-1]
    ## do not keep close similar hues
    oe <- seq_along(h) %% 2
    h <- h[order(oe, h)]
    
    hh <- rep(h, rep(minor, major))
    
    ## sr <- c(seq(from = 0.4, to = 0.5, length.out = ((minor - 1)/2) - 1),
    ##         1,
    ##         seq(from = 0.6, to = 0.8, length.out = (minor - 1)/2))

    ## sv <- c(seq(from = 0.4, to = 0.5, length.out = ((minor - 1)/2) - 1),
    ##         1,
    ##         seq(from = 0.6, to = 0.8, length.out = (minor - 1)/2))

    sr <- seq(from = 0.4, to = 1, length.out = minor)
    sv <- seq(from = 0.8, to = 1, length.out = minor)

    ## Nope, too many jumps
    ## oe <- seq_along(sr) %% 2
    ## sr <- sr[order(oe, sr)]
    ## oe <- seq_along(sv) %% 2
    ## sv <- sv[order(oe, sv)]
    ## hummm...
    sv <- rev(sv)
    
    srr <- rep(sr, major)
    svv <- rep(sv, major)

    col <- hsv(hh, srr, svv)

    plot(seq_along(hh), col = col, pch = 16, cex = 3)
    print(sr)
    print(sv)
}


myhsv2 <- function(major = 5, minor = 6,
                   srange = c(0.4, 1),
                   vrange = c(0.8, 1)) {
    ## major is a number, minor can be a vector
    if(length(minor) == 1) minor <- rep(minor, major)

    h <- seq(from = 0, to = 1, length.out = major + 1)[-1]
    ## do not keep close similar hues
    oe <- seq_along(h) %% 2
    h <- h[order(oe, h)]

    hh <- rep(h, minor)
    
    sr <- unlist(lapply(minor, function(x) 
        seq(from = srange[1], to = srange[2], length.out = x)))
    sv <- unlist(lapply(minor, function(x) 
        rev(seq(from = vrange[1], to = vrange[2], length.out = x)))
        )

    col <- hsv(hh, srr, svv)
    
    plot(1:(sum(minor)), col = col, pch = 16, cex = 3)
    print(sr)
    print(sv)
}



clcl <- function(ramp) {
    apply(ramp, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
}

dwl <- function(time, y, ndr, genotLabels) {
    ## Put data in long format, for ggplot et al
    nc <- ncol(y)
    nr <- nrow(y)
    y <- as.vector(y)
    df <- data.frame(Time = rep(time, nc),
                     Y = y,
                     Drivers = rep(ndr, rep(nr, nc)),
                     Genotype = rep(genotLabels, rep(nr, nc)))
    return(df)
}

dwla <- function(x) {
    if(!inherits(x, "oncosimul2"))
        ndr <- colSums(x$Genotypes[1:x$NumDrivers, , drop = FALSE])
    else {
        ndr <- colSums(x$Genotypes[x$Drivers, , drop = FALSE])
    }

    y <- x$pops.by.time[, 2:ncol(x$pops.by.time), drop = FALSE]

    oo <- order(ndr)
    y <- y[, oo, drop = FALSE]
    ndr <- ndr[oo]
    
    dwl(x$pops.by.time[, 1],
        y,
        ndr,
        x$GenotypesLabels)

}

