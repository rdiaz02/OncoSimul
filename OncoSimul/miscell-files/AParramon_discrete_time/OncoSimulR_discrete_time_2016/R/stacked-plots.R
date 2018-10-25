## Copyright 2013, 2014, 2015, 2016 Marc Taylor

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.


## Functions for stacked and stream plots. Originals from Mark Taylor; see
## https://github.com/marchtaylor/sinkr and
## http://menugget.blogspot.com.es/2013/12/data-mountains-and-streams-stacked-area.html




## plot.stream makes a "stream plot" where each y series is plotted 
## as stacked filled polygons on alternating sides of a baseline.

## Arguments include:
## 'x' - a vector of values
## 'y' - a matrix of data series (columns) corresponding to x
## 'order.method' = c("as.is", "max", "first") 
##  "as.is" - plot in order of y column
##  "max" - plot in order of when each y series reaches maximum value
##  "first" - plot in order of when each y series first value > 0
## 'center' - if TRUE, the stacked polygons will be centered so that the middle,
## i.e. baseline ("g0"), of the stream is approximately equal to zero. 
## Centering is done before the addition of random wiggle to the baseline. 
## 'frac.rand' - fraction of the overall data "stream" range used to define the range of
## random wiggle (uniform distrubution) to be added to the baseline 'g0'
## 'spar' - setting for smooth.spline function to make a smoothed version of baseline "g0"
## 'col' - fill colors for polygons corresponding to y columns (will recycle)
## 'border' - border colors for polygons corresponding to y columns (will recycle) (see ?polygon for details)
## 'lwd' - border line width for polygons corresponding to y columns (will recycle)
## '...' - other plot arguments

plot.stream2 <- function(
                        x, y, 
                        order.method = "as.is", frac.rand=0.1, spar=0.2,
                        center=TRUE,
                        ylab="", xlab="",  
                        border = NULL, lwd=1, 
                        col=rainbow(length(y[1,])),
                        ylim=NULL,
                        log = "",
                        ...
                        ){

    if(sum(y < 0, na.rm = TRUE) > 0) stop("y cannot contain negative numbers")

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
        ord <- order(apply(y, 2, function(x) min(which(x>0), na.rm = TRUE)))
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

    ylim.tmp <- range(sapply(polys,
                             function(x) range(x$y, na.rm=TRUE)), na.rm=TRUE)
    outer.lims <- sapply(polys,
                         function(r) rev(r$y[(length(r$y)/2+1):length(r$y)]))
    mid <- apply(outer.lims, 1,
                 function(r) mean(c(max(r, na.rm=TRUE),
                                    min(r, na.rm=TRUE)), na.rm=TRUE))
    
    ## center and wiggle
    if(center) {
        g0 <- -mid + runif(length(x),
                           min=frac.rand*ylim.tmp[1],
                           max=frac.rand*ylim.tmp[2])
    } else {
        g0 <- runif(length(x),
                    min=frac.rand*ylim.tmp[1],
                    max=frac.rand*ylim.tmp[2])
    }
    fit <- smooth.spline(g0 ~ x, spar=spar)

    for(i in seq(polys)){
        polys[[i]]$y <- polys[[i]]$y + c(fit$y, rev(fit$y))
    }

    if(is.null(ylim)) ylim <- range(sapply(polys,
                                           function(x) range(x$y, na.rm=TRUE)),
                                    na.rm=TRUE)
    if(grepl("x", log))
        axes <- FALSE
    else
        axes <- TRUE
    plot(x,y[,1], ylab=ylab, xlab=xlab, ylim=ylim, t="n", axes = axes, ...)
    for(i in seq(polys)){
        polygon(polys[[i]], border=border[i], col=col[i], lwd=lwd[i])
    }
    if(!axes) {
        ## yes, we only allow transformation of x axis
        relabelLogaxis(1)
        axis(2)
    }
}



## plot.stacked makes a stacked plot where each y series is plotted on top
## of the each other using filled polygons

## Arguments include:
## 'x' - a vector of values
## 'y' - a matrix of data series (columns) corresponding to x
## 'order.method' = c("as.is", "max", "first") 
##  "as.is" - plot in order of y column
##  "max" - plot in order of when each y series reaches maximum value
##  "first" - plot in order of when each y series first value > 0
## 'col' - fill colors for polygons corresponding to y columns (will recycle)
## 'border' - border colors for polygons corresponding to y columns (will recycle) (see ?polygon for details)
## 'lwd' - border line width for polygons corresponding to y columns (will recycle)
## '...' - other plot arguments

plot.stacked2 <- function(
                         x, y, 
                         order.method = "as.is",
                         ylab="", xlab="", 
                         border = NULL, lwd=1, 
                         col=rainbow(length(y[1,])),
                         ylim=NULL,
                         log = "",
                         ...){
    if(sum(y < 0) > 0) stop("y cannot contain negative numbers")

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
        ord <- order(apply(y, 2, function(x) min(which(x>0))))
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

    if(is.null(ylim)) ylim <- range(sapply(polys,
                                           function(x) range(x$y, na.rm=TRUE)),
                                    na.rm=TRUE)
    if(grepl("x", log))
        axes <- FALSE
    else
        axes <- TRUE
    plot(x,y[,1], ylab=ylab, xlab=xlab, ylim=ylim, t="n",  axes = axes, ...)
    for(i in seq(polys)){
        polygon(polys[[i]], border=border[i], col=col[i], lwd=lwd[i])
    }
    if(!axes) {
        ## yes, we only allow transformation of x axis
        relabelLogaxis(1)
        axis(2)
    }
}
