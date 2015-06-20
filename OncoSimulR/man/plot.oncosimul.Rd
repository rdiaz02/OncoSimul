\name{plot.oncosimul}
\alias{plot.oncosimul}
\alias{plot.oncosimulpop}


\title{
  Plot simulated tumor progression data.
}
\description{
  Plots data generated from the simulations, either for a single
  individual or for a population of individuals, with time units in the
  x axis and nubmer of cells in the y axis. By default, all clones
  with the same number of drivers are plotted using the same colour (but
  different line types), and clones with different number of drivers are
  plotted in different colours.
}

\usage{

\method{plot}{oncosimul}(x, col = c(8, "orange", 6:1), log = "y",
                   ltyClone = 2:6, lwdClone = 0.9,
                   ltyDrivers = 1, lwdDrivers = 3,
                   xlab = "Time units",
                   ylab = "Number of cells", plotClones = TRUE,
                   plotDrivers = TRUE, addtot = FALSE,
                   addtotlwd = 0.5, yl = NULL, thinData = FALSE,
                   thinData.keep = 0.1, thinData.min = 2,
                   plotDiversity = FALSE, ...)

\method{plot}{oncosimulpop}(x, ask = TRUE, col = c(8, "orange", 6:1),
                      log = "y",
                      ltyClone = 2:6, lwdClone = 0.9,
                      ltyDrivers = 1, lwdDrivers = 3,
                      xlab = "Time units",
                      ylab = "Number of cells", plotClones = TRUE,
                      plotDrivers = TRUE, addtot = FALSE,
                      addtotlwd = 0.5, yl = NULL, thinData = FALSE,
                      thinData.keep = 0.1, thinData.min = 2,
                      plotDiversity = FALSE, ...)



}
\arguments{

  \item{x}{
    An object of class \code{oncosimul} (for \code{plot.oncosimul}) or
    \code{oncosimulpop} (for \code{plot.oncosimulpop}). 
  }
  
  \item{ask}{Same meaning as in \code{\link{par}}.    
  }
  
  \item{col}{
    Colour of the lines, where each type of clone (where type is defined
    by number of drivers) has a different color. If there are many
    drivers, col is recycled, so you might want to increase the number
    of possible colours.
  }
  \item{log}{
    See \code{log} in \code{\link{plot.default}}. The default, "y", will
    make the y axis logarithmic.
  }
  
  \item{ltyClone}{
    
    Line type for each clone. Recycled as needed. You probably do not
    want to use lty=1 for any clone, to differentiate from the clone
    type, unless you change the setting for \code{ltyDrivers}.
    
  }
  
  \item{lwdClone}{
    Line width for clones.
    
  }
  \item{ltyDrivers}{
    Line type for the driver type.
    
  }
  \item{lwdDrivers}{
    Line width for the driver type.
  }
  \item{xlab}{
    Same as \code{xlab} in \code{\link{plot.default}}.
    
  }
  \item{ylab}{
    Same as \code{ylab} in \code{\link{plot.default}}.
    
  }
  \item{plotClones}{
    Should clones be plotted?
    
  }
  \item{plotDrivers}{
    Should clone types (which are defined by number of drivers), be plotted?
    
  }
  \item{addtot}{
    If TRUE, add a line with the total populatino size.
    
  }
  \item{addtotlwd}{
    Line width for total population size.
    
  }
  \item{yl}{
    If non NULL, limits of the y axis. Same as in
    \code{\link{plot.default}}. If NULL, the limits are calculated
    automatically.
    
  }
  \item{thinData}{
    If TRUE, the data plotted is a subset of the original data. The
    original data are "thinned" in such a way that the origin of each
    clone is not among the non-shown data (i.e., so that we can see when
    each clone/driver originates).
    
    Thining is done to reduce the plot size and to speed up plotting..
    
  }
  \item{thinData.keep}{
    The fraction of the data to keep (actually, a lower bound on the
    fraction of data to keep).
    
    
  }
  \item{thinData.min}{
    Any time point for which a clone has a population size < thinData.min
    will be kept (i.e., will not be removed from) in the data.
    
    
  }

  \item{plotDiversity}{If TRUE, we also show, on top of the main figure,
  Shannon's diversity index (and we considers as distinct those
  genotypes with different order of mutations when order matters).}

  
  \item{\dots}{
    Other arguments passed to \code{plots}. For instance, \code{main}.
    
  }
}
\author{
  Ramon Diaz-Uriarte
  
}



\seealso{
  \code{\link{oncoSimulIndiv}}
}
\examples{
data(examplePosets)
p701 <- examplePosets[["p701"]]

## Simulate and plot a single individual, including showing
## Shannon's diversity index
b1 <- oncoSimulIndiv(p701)
plot(b1, addtot = TRUE, plotDiversity = TRUE)


## simulate and plot 2 individuals
## (I set mc.cores = 2 to comply with --as-cran checks, but you
##  should either use a reasonable number for your hardware or
##  leave it at its default value).

p1 <- oncoSimulPop(2, p701, mc.cores = 2)

par(mfrow = c(1, 2))
plot(p1, ask = FALSE)

}


\keyword{hplot}


