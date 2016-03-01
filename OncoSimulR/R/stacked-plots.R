## Functions for stacked and stream plots. Originals from Mark Taylor; see
## https://github.com/marchtaylor/sinkr and
## http://menugget.blogspot.com.es/2013/12/data-mountains-and-streams-stacked-area.html


#plot.stream makes a "stream plot" where each y series is plotted 
#as stacked filled polygons on alternating sides of a baseline.
#
#Arguments include:
#'x' - a vector of values
#'y' - a matrix of data series (columns) corresponding to x
#'order.method' = c("as.is", "max", "first") 
#  "as.is" - plot in order of y column
#  "max" - plot in order of when each y series reaches maximum value
#  "first" - plot in order of when each y series first value > 0
#'center' - if TRUE, the stacked polygons will be centered so that the middle,
#i.e. baseline ("g0"), of the stream is approximately equal to zero. 
#Centering is done before the addition of random wiggle to the baseline. 
#'frac.rand' - fraction of the overall data "stream" range used to define the range of
#random wiggle (uniform distrubution) to be added to the baseline 'g0'
#'spar' - setting for smooth.spline function to make a smoothed version of baseline "g0"
#'col' - fill colors for polygons corresponding to y columns (will recycle)
#'border' - border colors for polygons corresponding to y columns (will recycle) (see ?polygon for details)
#'lwd' - border line width for polygons corresponding to y columns (will recycle)
#'...' - other plot arguments

plot.stream <- function(
	x, y, 
	order.method = "as.is", frac.rand=0.1, spar=0.2,
	center=TRUE,
	ylab="", xlab="",  
	border = NULL, lwd=1, 
	col=rainbow(length(y[1,])),
	ylim=NULL, 
	...
){

    if(sum(y < 0, na.rm = TRUE) > 0) error("y cannot contain negative numbers")

    ## if(log %in% c("y", "xy", "yx") ) {
        

    ## }
    
    if(is.null(border)) border <- par("fg")
    border <- as.vector(matrix(border, nrow=ncol(y), ncol=1))
    col <- as.vector(matrix(col, nrow=ncol(y), ncol=1))
    lwd <- as.vector(matrix(lwd, nrow=ncol(y), ncol=1))
    
	if(order.method == "max") {
		ord <- order(apply(y, 2, which.max, na.rm = TRUE))
		y <- y[, ord]
		col <- col[ord]
		border <- border[ord]
	}

	if(order.method == "first") {
		ord <- order(apply(y, 2, function(x) min(which(r>0), na.rm = TRUE)))
		y <- y[, ord]
		col <- col[ord]
		border <- border[ord]
	}

	bottom.old <- x*0
	top.old <- x*0
	polys <- vector(mode="list", ncol(y))
	for(i in seq(polys)){
		if(i %% 2 == 1){ #if odd
			top.new <- top.old + y[,i]
			polys[[i]] <- list(x=c(x, rev(x)), y=c(top.old, rev(top.new)))
			top.old <- top.new
		}
		if(i %% 2 == 0){ #if even
			bottom.new <- bottom.old - y[,i]
			polys[[i]] <- list(x=c(x, rev(x)), y=c(bottom.old, rev(bottom.new)))
			bottom.old <- bottom.new
		}
	}

	ylim.tmp <- range(sapply(polys, function(x) range(x$y, na.rm=TRUE)), na.rm=TRUE)
	outer.lims <- sapply(polys, function(r) rev(r$y[(length(r$y)/2+1):length(r$y)]))
	mid <- apply(outer.lims, 1, function(r) mean(c(max(r, na.rm=TRUE), min(r, na.rm=TRUE)), na.rm=TRUE))
	
	#center and wiggle
	if(center) {
		g0 <- -mid + runif(length(x), min=frac.rand*ylim.tmp[1], max=frac.rand*ylim.tmp[2])
	} else {
		g0 <- runif(length(x), min=frac.rand*ylim.tmp[1], max=frac.rand*ylim.tmp[2])
	}
	browser()
	fit <- smooth.spline(g0 ~ x, spar=spar)

	for(i in seq(polys)){
		polys[[i]]$y <- polys[[i]]$y + c(fit$y, rev(fit$y))
	}

	if(is.null(ylim)) ylim <- range(sapply(polys, function(x) range(x$y, na.rm=TRUE)), na.rm=TRUE)
	plot(x,y[,1], ylab=ylab, xlab=xlab, ylim=ylim, t="n", ...)
	for(i in seq(polys)){
		polygon(polys[[i]], border=border[i], col=col[i], lwd=lwd[i])
	}

}



#plot.stacked makes a stacked plot where each y series is plotted on top
#of the each other using filled polygons
#
#Arguments include:
#'x' - a vector of values
#'y' - a matrix of data series (columns) corresponding to x
#'order.method' = c("as.is", "max", "first") 
#  "as.is" - plot in order of y column
#  "max" - plot in order of when each y series reaches maximum value
#  "first" - plot in order of when each y series first value > 0
#'col' - fill colors for polygons corresponding to y columns (will recycle)
#'border' - border colors for polygons corresponding to y columns (will recycle) (see ?polygon for details)
#'lwd' - border line width for polygons corresponding to y columns (will recycle)
#'...' - other plot arguments

plot.stacked <- function(
	x, y, 
	order.method = "as.is",
	ylab="", xlab="", 
	border = NULL, lwd=1, 
	col=rainbow(length(y[1,])),
	ylim=NULL,
	...
){

    	if(sum(y < 0) > 0) error("y cannot contain negative numbers")

	if(is.null(border)) border <- par("fg")
	border <- as.vector(matrix(border, nrow=ncol(y), ncol=1))
	col <- as.vector(matrix(col, nrow=ncol(y), ncol=1))
	lwd <- as.vector(matrix(lwd, nrow=ncol(y), ncol=1))

	if(order.method == "max") {
		ord <- order(apply(y, 2, which.max))
		y <- y[, ord]
		col <- col[ord]
		border <- border[ord]
	}

	if(order.method == "first") {
		ord <- order(apply(y, 2, function(x) min(which(r>0))))
		y <- y[, ord]
		col <- col[ord]
		border <- border[ord]
	}

	top.old <- x*0
	polys <- vector(mode="list", ncol(y))
	for(i in seq(polys)){
		top.new <- top.old + y[,i]
		polys[[i]] <- list(x=c(x, rev(x)), y=c(top.old, rev(top.new)))
		top.old <- top.new
	}

	if(is.null(ylim)) ylim <- range(sapply(polys, function(x) range(x$y, na.rm=TRUE)), na.rm=TRUE)
	plot(x,y[,1], ylab=ylab, xlab=xlab, ylim=ylim, t="n", ...)
	for(i in seq(polys)){
		polygon(polys[[i]], border=border[i], col=col[i], lwd=lwd[i])
	}

}



## FIXME: do we ever use log of x?

plotClonesSt <- function(z, ndr = NULL, na.subs = TRUE,
                         log = "y",
                         lwd = lwd,
                         ## type = "l",
                         lty = 1:8, col = 1:9,
                         order.method = order.method,
                         stream.center = stream.center,
                         stream.frac.rand = stream.frac.rand,
                         stream.spar = stream.spar,
                         border = border,
                         srange = srange,
                         vrange = vrange,
                         type = type,
                         ...) {

    ## if given ndr, we order columns based on ndr, so clones with more
    ## drivers are plotted last

    y <- z$pops.by.time[, 2:ncol(z$pops.by.time), drop = FALSE]

    ## How best to deal with log with stacked and stream still not clear.
    ## So for now leave as this.
    if(type %in% c("stacked", "stream") )
        na.subs <- FALSE
    
    if(na.subs){
        y[y == 0] <- NA
    }

    browser()
    
    if(!is.null(ndr)) {
        ## could be done above, to avoid creating
        ## more copies
        oo <- order(ndr)
        y <- y[, oo, drop = FALSE]
        ndr <- ndr[oo]
        col <- col[ndr + 1]
    }
    ## FIXME: break if ndr is null?!
    if(type == "line") {
        matplot(x = z$pops.by.time[, 1],
                y = y,
                log = log, type = "l",
                col = col, lty = lty,
                lwd = lwd,
                ...)
        box()
    } else {
        cll <- myhsvcols(ndr, srange = srange, vrange = vrange)
        x <- z$pops.by.time[, 1]
        if(grepl("y", log)) {
            y <- log10(y + 1)
        }
        if(grepl("x", log)) {
            x <- log10(x + 1)
        }

        if (type == "stacked") {
            browser()
            plot.stacked(x = x,
                         y = y,
                         order.method = order.method,
                         border = border,
                         lwd = lwd,
                         col = cll$colors,
                         log = log) ## log thing
            ## legend thing
        } else if (type == "stream") {
            plot.stream(x = x,
                        y = y,
                        order.method = order.method,
                        border = border,
                        lwd = lwd,
                        col = cll$colors,
                        frac.rand = stream.frac.rand,
                        spar = stream.spar,
                        center = stream.center,
                        log = log)
            ## legend thing
        }
        browser()
        legend(x = "topleft",
               title = "Number of drivers",
               lty = 1,
               col = cll$colorsLegend$Color,
               legend = cll$colorsLegend$Drivers)
        
    }
}



mylogaxis <- function(side = 2) {
    po <- axis( side = side, labels = FALSE, tick = FALSE)
    labels <- 10^po
    

}




plot.oncosimul2 <- function(x,
                           type = "stacked",
                           col = "auto",
                           log = "y",
                           ltyClone = 2:6,
                           lwdClone = 0.9,
                           ltyDrivers = 1,
                           lwdDrivers = 3,
                           xlab = "Time units",
                           ylab = "Number of cells",
                           plotClones = TRUE,
                           plotDrivers = TRUE,
                           addtot = FALSE,
                           addtotlwd = 0.5,
                           yl = NULL,
                           thinData = FALSE,
                           thinData.keep = 0.1,
                           thinData.min = 2,
                           plotDiversity = FALSE,
                           order.method = "as.is",
                           stream.center = TRUE,
                           stream.frac.rand = 0.01,
                           stream.spar = 0.2,
                           border = NULL,
                           lwd = 1,
                           srange = c(0.4, 1),
                           vrange = c(0.8, 1),
                           ...
                           ) {

    ## FIXME: test this
    if(!(type %in% c("stacked", "stream", "line")))
        stop("Type of plot unknown: it must be one of",
             "stacked, stream or line")
    
    if(col == "auto" && (type == "line") )
        col = c(8, "orange", 6:1)
    
    
    if(thinData)
        x <- thin.pop.data(x, keep = thinData.keep, min.keep = thinData.min)

    ## uvx
    if(!inherits(x, "oncosimul2"))
        ndr <- colSums(x$Genotypes[1:x$NumDrivers, , drop = FALSE])
    else {
        ndr <- colSums(x$Genotypes[x$Drivers, , drop = FALSE])
    }
    
    if(is.null(yl)) {
        if(log %in% c("y", "xy", "yx") )
            yl <- c(1, max(apply(x$pops.by.time[, -1, drop = FALSE], 1, sum)))
        else
            yl <- c(0, max(apply(x$pops.by.time[, -1, drop = FALSE], 1, sum)))
    }
    if(plotDiversity) {
        par(fig = c(0, 1, 0.8, 1))
        m1 <- par()$mar
        m <- m1
        m[c(1, 3)] <- c(0, 0.7)
        op <- par(mar = m )
        plotShannon(x)
        par(op)
        m1[c(3)] <- 0.2
        op <- par(mar = m1)
        par(fig = c(0, 1, 0, 0.8), new = TRUE)  
    }
    if(plotClones) {
        plotClonesSt(x,
                     ndr = ndr, 
                     xlab = xlab,
                     ylab = ylab,
                     lty = ltyClone,
                     col = col, 
                     ylim = yl,
                     lwd = lwdClone,
                     axes = FALSE,
                     log = log,
                     order.method = order.method,
                     stream.center = stream.center,
                     stream.frac.rand = stream.frac.rand,
                     stream.spar = stream.spar,
                     border = border,
                     srange = srange,
                     vrange = vrange,
                     type = type,
                     ...)
    }

    if(plotClones && plotDrivers && (type == "line"))
        par(new = TRUE)
    
    if(plotDrivers && (type == "line")){
        plotDrivers0(x,
                     ndr,
                     timescale = 1,
                     trim.no.drivers = FALSE,
                     xlab = "", ylab = "",
                     lwd = lwdDrivers,
                     lty = ltyDrivers,
                     col = col, 
                     addtot = addtot,
                     addtotlwd = addtotlwd,
                     log = log, ylim = yl,
                     ...)
    }
    
}






########## colors, etc

### Notes for me

library(RColorBrewer)

display.brewer.all()

## For genotypes, use "categorical", like "Set1", "Paired", or "Dark2"



## for distinct colors, from gplots
rich.colors(5)
?rich.colors ## from gplots
plot(1:8, col = rich.colors(8), pch = 16, cex = 3)
## From the example on the colorspace vignette
plot(1:6, col = rainbow_hcl(6, c = 60, l = 75), pch = 16, cex = 3)


myhsvcols <- function(ndr, srange = c(0.4, 1),
                      vrange = c(0.8, 1)) {
    ## Generate a set of colors so that:
    ##  - easy to tell when we increase number of drivers
    ##  - reasonably easy to generate a legend
    ##  - different clones with same number of drivers have "similar" colors

    ## I use hsv color specification as this seems the easiest.
    
    minor <- table(ndr)
    major <- length(minor)
    
    h <- seq(from = 0, to = 1, length.out = major + 1)[-1]
    ## do not keep close similar hues
    oe <- seq_along(h) %% 2
    h <- h[order(oe, h)]

    hh <- rep(h, minor)
    
    sr <- unlist(lapply(minor, function(x) 
        seq(from = srange[1], to = srange[2], length.out = x)))
    sv <- unlist(lapply(minor, function(x) 
        seq(from = vrange[1], to = vrange[2], length.out = x))
        )

    colors <- hsv(hh, sr, sv)
 
    colorsLegend <- aggregate(list(Color = colors), list(Drivers = ndr),
                              function(x)
                                  as.character(x[((length(x) %/% 2) + 1 )]))
    ## To show what it would look like
    ## plot(1:(sum(minor)), col = colors, pch = 16, cex = 3)
    ## legend(1, length(ndr), col = colorsLegend$Color, legend = names(minor),
    ##        pch = 16)

    return(list(colors = colors,
                colorsLegend = colorsLegend))
}


plot(1:6, col = hsv(seq(from = 0, to = 1, length.out = 10), 1, 1), pch = 16, cex = 3)



## look at these examples for base
## http://stackoverflow.com/questions/27250542/how-to-make-gradient-color-filled-timeseries-plot-in-r


##   FIXME
#### Get things ready to use streamgraph, and add example to vignette


############### with streamgraph package



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




colors <- function(ndr) {
    tt <- table(ndr)
    
    nl <- length(unique(ndr))
    

}





n <- 20;  y <- -sin(3*pi*((1:n)-1/2)/n)
op <- par(mar = rep(1.5, 4))
plot(y, axes = FALSE, frame.plot = TRUE,
     xlab = "", ylab = "", pch = 21, cex = 30,
     bg = rainbow(n, start = .8, end = .5),
     main = "Red tones")


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




np <- 8
nr <- 2
pl <- "Dark2"
## plot(1:np, col = brewer.pal(np,pl)[-(seq(1:nr))])
plot(1:np, col = brewer.pal(np,pl))

pp <- brewer.pal(8,pl)
np2 <- 8

plot(1:np2, col = colorRampPalette(pp[c(1, 2)])( np2))

cl <- colorRamp(brewer.pal(8,pl), alpha = TRUE)( (0:20)/20 )


plot(1:20, col = clcl(cl), pch = 16, cex = 3)




clcl <- function(ramp) {
    apply(ramp, 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
}





    


    
    
plotClones3 <- function(z, ndr = NULL, na.subs = TRUE,
                       log = "y", type = "l",
                       lty = 1:8, col = 1:9, ...) {

    ## if given ndr, we order columns based on ndr, so clones with more
    ## drivers are plotted last

    y <- z$pops.by.time[, 2:ncol(z$pops.by.time), drop = FALSE]
    
    ## if(na.subs){
    ##     y[y == 0] <- NA
    ## }
    
    if(!is.null(ndr)) {
        ## could be done above, to avoid creating
        ## more copies
        oo <- order(ndr)
        y <- y[, oo, drop = FALSE]
        ndr <- ndr[oo]
        col <- col[ndr + 1]
    }
    plot.stream(x = z$pops.by.time[, 1],
                y = y,
                frac.rand = 0) ## work on the log axis
    ## matplot(x = z$pops.by.time[, 1],
    ##         y = y,
    ##         log = log, type = type,
    ##         col = col, lty = lty,
    ##         ...)
    ## box()
}

plotClones2 <- function(z, ndr = NULL, na.subs = TRUE,
                       log = "y", type = "l",
                       lty = 1:8, col = 1:9, ...) {

    ## if given ndr, we order columns based on ndr, so clones with more
    ## drivers are plotted last

    y <- z$pops.by.time[, 2:ncol(z$pops.by.time), drop = FALSE]
    
    ## if(na.subs){
    ##     y[y == 0] <- NA
    ## }
    
    if(!is.null(ndr)) {
        ## could be done above, to avoid creating
        ## more copies
        oo <- order(ndr)
        y <- y[, oo, drop = FALSE]
        ndr <- ndr[oo]
        col <- col[ndr + 1]
    }
    plot.stacked(x = z$pops.by.time[, 1],
                 y = y) ## work on the log axis
    ## matplot(x = z$pops.by.time[, 1],
    ##         y = y,
    ##         log = log, type = type,
    ##         col = col, lty = lty,
    ##         ...)
    ## box()
}


plot.oncosimul2 <- function(x, col = c(8, "orange", 6:1),
                            log = "y",
                            ltyClone = 2:6,
                            lwdClone = 0.9,
                           ltyDrivers = 1,
                           lwdDrivers = 3,
                           xlab = "Time units",
                           ylab = "Number of cells",
                           plotClones = TRUE,
                           plotDrivers = TRUE,
                           addtot = FALSE,
                           addtotlwd = 0.5,
                           yl = NULL,
                           thinData = FALSE,
                           thinData.keep = 0.1,
                           thinData.min = 2,
                           plotDiversity = FALSE,
                           ...
                           ) {

    if(thinData)
        x <- thin.pop.data(x, keep = thinData.keep, min.keep = thinData.min)

    ## uvx
    if(!inherits(x, "oncosimul2"))
        ndr <- colSums(x$Genotypes[1:x$NumDrivers, , drop = FALSE])
    else {
        ndr <- colSums(x$Genotypes[x$Drivers, , drop = FALSE])
    }
    
    if(is.null(yl)) {
        if(log %in% c("y", "xy", "yx") )
            yl <- c(1, max(apply(x$pops.by.time[, -1, drop = FALSE], 1, sum)))
        else
            yl <- c(0, max(apply(x$pops.by.time[, -1, drop = FALSE], 1, sum)))
    }
    ## if(plotDiversity) {
    ##     par(fig = c(0, 1, 0.8, 1))
    ##     m1 <- par()$mar
    ##     m <- m1
    ##     m[c(1, 3)] <- c(0, 0.7)
    ##     op <- par(mar = m )
    ##     plotShannon(x)
    ##     par(op)
    ##     m1[c(3)] <- 0.2
    ##     op <- par(mar = m1)
    ##     par(fig = c(0, 1, 0, 0.8), new = TRUE)  
    ## }
    
    ##    if(plotClones) {

    
    plotClones2(x,
               ndr = ndr, 
               xlab = xlab,
               ylab = ylab,
               lty = ltyClone,
               col = col, 
               ylim = yl,
               lwd = lwdClone,
               axes = FALSE,
               log = log,
               ...)
##    }

    ## if(plotClones && plotDrivers)
    ##     par(new = TRUE)
    
    ## if(plotDrivers){
    ##     plotDrivers0(x,
    ##                  ndr,
    ##                  timescale = 1,
    ##                  trim.no.drivers = FALSE,
    ##                  xlab = "", ylab = "",
    ##                  lwd = lwdDrivers,
    ##                  lty = ltyDrivers,
    ##                  col = col, 
    ##                  addtot = addtot,
    ##                  addtotlwd = addtotlwd,
    ##                  log = log, ylim = yl,
    ##                  ...)
    }
