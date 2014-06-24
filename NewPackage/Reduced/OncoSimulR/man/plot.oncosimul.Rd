\name{plot.oncosimul}
\alias{plot.oncosimul}
\alias{plot.oncosimulpop}

%- Also NEED an '\alias' for EACH other topic documented here.
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
                   ltyClone = 2:6, lwdClone = 0.2,
                   ltyDrivers = 1, lwdDrivers = 3,
                   xlab = "Time units",
                   ylab = "Number of cells", plotClones = TRUE,
                   plotDrivers = TRUE, addtot = FALSE,
                   addtotlwd = 0.5, yl = NULL, thinData = FALSE,
                   thinData.keep = 0.1, thinData.min = 2, ...)

\method{plot}{oncosimulpop}(x, ask = TRUE, col = c(8, "orange", 6:1),
                      log = "y",
                      ltyClone = 2:6, lwdClone = 0.2,
                      ltyDrivers = 1, lwdDrivers = 3,
                      xlab = "Time units",
                      ylab = "Number of cells", plotClones = TRUE,
                      plotDrivers = TRUE, addtot = FALSE,
                      addtotlwd = 0.5, yl = NULL, thinData = FALSE,
                      thinData.keep = 0.1, thinData.min = 2, ...)



}
%- maybe also 'usage' for other objects documented here.
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
    %%     ~~Describe \code{ltyDrivers} here~~
}
\item{lwdDrivers}{
  Line width for the driver type.
%%     ~~Describe \code{lwdDrivers} here~~
}
  \item{xlab}{
    Same as \code{xlab} in \code{\link{plot.default}}.
    %%     ~~Describe \code{xlab} here~~
}
\item{ylab}{
      Same as \code{ylab} in \code{\link{plot.default}}.
%%     ~~Describe \code{ylab} here~~
}
\item{plotClones}{
  Should clones be plotted?
%%     ~~Describe \code{plotClones} here~~
}
\item{plotDrivers}{
  Should clone types (which are defined by number of drivers), be plotted?
%%     ~~Describe \code{plotDrivers} here~~
}
\item{addtot}{
  If TRUE, add a line with the total populatino size.
%%     ~~Describe \code{addtot} here~~
}
\item{addtotlwd}{
  Line width for total population size.
%%     ~~Describe \code{addtotlwd} here~~
}
\item{yl}{
  If non NULL, limits of the y axis. Same as in
  \code{\link{plot.default}}. If NULL, the limits are calculated automatically.
%%     ~~Describe \code{yl} here~~
}
\item{thinData}{
  If TRUE, the data plotted is a subset of the original data. The
  original data are "thinned" in such a way that the origin of each
  clone is not among the non-shown data (i.e., so that we can see when
  each clone/driver originates).

  Thining is done to reduce the plot size and to speed up plotting..
%%     ~~Describe \code{thinPops} here~~
}
  \item{thinData.keep}{
    The fraction of the data to keep (actually, a lower bound on the
    fraction of data to keep).

    %%     ~~Describe \code{thinPops.keep} here~~
}
\item{thinData.min}{
  Any time point for which a clone has a population size < thinData.min
  will be kept (i.e., will not be removed from) in the data.
  
%%     ~~Describe \code{thinPops.min} here~~
}
\item{\dots}{
  Other arguments passed to \code{plots}. For instance, \code{main}.
%%     ~~Describe \code{\dots} here~~
}
}
\author{
  Ramon Diaz-Uriarte
  %%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{oncoSimulIndiv}}
}
\examples{
data(examplePosets)
p701 <- examplePosets[["p701"]]

## simulate and plot a single individual
b1 <- oncoSimulIndiv(p701, keepEvery = 5,
                     sampleEvery = 5,
                     detectionSize = 1e8)
plot(b1, addtot = TRUE)


## simulate and plot 9 individuals 
p1 <- oncoSimulPop(9, p701, detectionSize = 1e7,
                   sampleEvery = 10, mc.cores = 4)

par(mfrow = c(3, 3))
plot(p1, ask = FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}


