\name{plot.oncosimul}
\alias{plot.oncosimul}
\alias{plot.oncosimulpop}


\title{
  Plot simulated tumor progression data.
}
\description{
  Plots data generated from the simulations, either for a single
  individual or for a population of individuals, with time units in the
  x axis and nubmer of cells in the y axis.

  In "drivers" plots, by default, all clones with the same number of
  drivers are plotted using the same colour (but different line types),
  and clones with different number of drivers are plotted in different
  colours. Plots can alternatively display genotypes instead of drivers.

  Plots available are line plots, stacked area, and stream plots.
  
}

\usage{

\method{plot}{oncosimul}(x,
                           show = "drivers", 
                           type = ifelse(show == "genotypes",
                                         "stacked", "line"),
                           col = "auto",
                           log = ifelse(type == "line", "y", ""),
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
                           ylim = NULL,
                           xlim = NULL,
                           thinData = FALSE,
                           thinData.keep = 0.1,
                           thinData.min = 2,
                           plotDiversity = FALSE,
                           order.method = "as.is",
                           stream.center = TRUE,
                           stream.frac.rand = 0.01,
                           stream.spar = 0.2,
                           border = NULL,
                           lwdStackedStream = 1,
                           srange = c(0.4, 1),
                           vrange = c(0.8, 1),
                           breakSortColors = "oe",
                           legend.ncols = "auto", ...)

\method{plot}{oncosimulpop}(x,
                              ask = TRUE,
                              show = "drivers", 
                              type = ifelse(show == "genotypes",
                                            "stacked", "line"),
                              col = "auto",
                              log = ifelse(type == "line", "y", ""),
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
                              ylim = NULL,
                              xlim = NULL,
                              thinData = FALSE,
                              thinData.keep = 0.1,
                              thinData.min = 2,
                              plotDiversity = FALSE,
                              order.method = "as.is",
                              stream.center = TRUE,
                              stream.frac.rand = 0.01,
                              stream.spar = 0.2,
                              border = NULL,
                              lwdStackedStream = 1,
                              srange = c(0.4, 1),
                              vrange = c(0.8, 1),
                              breakSortColors = "oe",
                              legend.ncols = "auto",
                              ...)

}
\arguments{

  \item{x}{
    An object of class \code{oncosimul} (for \code{plot.oncosimul}) or
    \code{oncosimulpop} (for \code{plot.oncosimulpop}). 
  }
  
  \item{ask}{Same meaning as in \code{\link{par}}.    
  }

  \item{show}{
    One of "drivers" or "genotypes". If "drivers" the legend
    will reflect the number of drivers. If "genotypes" you will be shown
    genotypes. You probably want to limit "genotypes" to those cases
    where only a relatively small number of genotypes exist (or the plot
    will be an unmanageable mess). The default is "drivers".
  }

  \item{type}{One of "line", "stacked", "stream".

    If "line", you are shown lines for each genotype or clone. This
    means that to get an idea of the total population size you need to
    use \code{plotDrivers = TRUE} with \code{addtot = TRUE}, or do the
    visual calculation in your head.

    If "stacked" a stacked area plot. If "stream" a stream plot. Since
  these stack areas, you immediately get the total population. But that
  also means you cannot use \code{log}.

    The default is to use "line" for \code{show = "drivers"} and
    "stacked" for \code{show = "genotypes"}.
    
  }
  
  
  \item{col}{ Colour of the lines/areas. For \code{show = "drivers"}
    each type of clone (where type is defined by number of drivers) has
    a different color. For \code{show = "genotypes"} color refers to
    genotypes. The vector is recycled as needed.

    The default is "auto". If you have \code{show == "genotypes"} we
    start from the "Dark2" palette from \code{\link{brewer.pal}} in the
    \code{RColorBrewer} package and extend the palette via
    \code{\link{colorRampPalette}}.  For \code{show == "drivers"} and
    \code{type == "line"} we use a vector of eight colors (that are,
    then recycled as needed).  If you use "stacked" or "stream",
    however, instead of "line", then we generate colors via a HSV
    specification that tries to: a) make it easy to differentiate
    between different drivers (by not having like colors for adjacent
    numbers of drivers); b) make it easy to have a "representative"
    driver color while using sligtly different colors for different
    clones of a driver. See the code by doing
    \code{OncoSimulR:::myhsvcols}.

    You can specify your own vector of colors, but it will be ignored
  with \code{show == "drivers"}.
  }


  \item{log}{ See

    \code{log} in \code{\link{plot.default}}. The default is to have "y"
    for \code{type == "line"}, and that will make the y axis
    logarithmic. Stacked and stream area plots do not allow for
    logarithmic y axis (since those depend on the additivity of areas
    but \code{log(a + b) != log(a) + log(b)}).
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
    Should clone types (which are defined by number of drivers), be
  plotted? (Only applies when using \code{show = "drivers"}).
    
  }
  \item{addtot}{
    If TRUE, add a line with the total populatino size.
    
  }
  \item{addtotlwd}{
    Line width for total population size.
    
  }
  \item{ylim}{
    If non NULL, limits of the y axis. Same as in
    \code{\link{plot.default}}. If NULL, the limits are calculated
    automatically.
    
  }
  \item{xlim}{
    If non NULL, limits of the x axis. Same as in
    \code{\link{plot.default}}. If NULL, the limits are calculated
    automatically. Using a non-NULL range smaller than the range of observed
  values of time can also lead to speed ups of large figures (since we
  trim the data).    
  }
  
  \item{thinData}{
    If TRUE, the data plotted is a subset of the original data. The
    original data are "thinned" in such a way that the origin of each
    clone is not among the non-shown data (i.e., so that we can see when
    each clone/driver originates).
    
    Thinning is done to reduce the plot size and to speed up plotting.

    Note that thinning is carried out before dealing with the plot axis,
    so the actual number of points to be plotted could be a lot less (if
    you reduce the x-axis considerably) than those returned from the
    thinning. (In extreme cases this could lead to crashes when trying
    to use stream plots if, say, you end up plotting only three values).
    
  }
  \item{thinData.keep}{
    The fraction of the data to keep (actually, a lower bound on the
    fraction of data to keep).
    
    
  }
  \item{thinData.min}{
    Any time point for which a clone has a population size > thinData.min
    will be kept (i.e., will not be removed from) in the data.
    
    
  }

  \item{plotDiversity}{If TRUE, we also show, on top of the main figure,
  Shannon's diversity index (and we consider as distinct those
  genotypes with different order of mutations when order matters).

  If you set this to true, using \code{par(mfrow = c(2, 2))} and similar
  will not work (since we use \code{par(fig = )}) to display the
  diversity as the top plot).
}


\item{order.method}{ For stacked and stream plots. \code{c("as.is",
    "max", "first")}. \code{"as.is"}: plot in order of y column;
    \code{"max"}: plot in order of when each y series reaches maximum
    value. \code{"first"}: plot in order of when each y series first value > 0.
}

\item{stream.center}{ For stream plots. If TRUE, the stacked polygons
  will be centered so that the middle, i.e. baseline ("g0"), of the
  stream is approximately equal to zero.  Centering is done before the
  addition of random wiggle to the baseline.  }

\item{stream.frac.rand}{ For stream plots. Fraction of the overall data
  "stream" range used to define the range of random wiggle (uniform
  distribution) to be added to the baseline 'g0'.  }

\item{stream.spar}{Setting for smooth.spline function to make a smoothed
  version of baseline "g0".  }

\item{border}{For stacked and stream plots. Border colors for polygons
  corresponding to y columns (will recycle) (see \code{\link{polygon}}
  for details).  }

\item{lwdStackedStream}{border line width for polygons corresponding to
  y columns (will recycle).
}


\item{srange}{ Range of values of s in the HSV specification of colors
  (see \code{col} for details. Only applies when using "stacked" or
  "stream" plots and \code{col == "auto"}.)
}


\item{vrange}{ Range of values of v in the HSV specification of colors
  (see \code{col} for details.Only applies when using "stacked" or
  "stream" plots and \code{col == "auto"}. )
}


\item{breakSortColors}{How to try to minimize that similar colors be
  used for contiguous or nearby driver categories. The default is "oe"
  which resorts them in alternating way. The other two options are
  "distave", where we alternate after folding from the mean and "random"
  where the colors are randomly sorted. Only applies when using
  "stacked" or "stream" plots and \code{col == "auto"}.  }

\item{legend.ncols}{The number of columns of the legend. If
  \code{"auto"} (the default), will have one column for six or less
  entries, and two for more than six.
}


 
  \item{\dots}{
    Other arguments passed to \code{plots}. For instance, \code{main}.
    
  }
}
\author{
  Ramon Diaz-Uriarte. Marc Taylor for stacked and stream plots.
  
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

## A stacked area plot
plot(b1, type = "stacked", plotDiversity = TRUE)

## And what if I show a stream plot?
plot(b1, type = "stream", plotDiversity = TRUE)

## Simulate and plot 2 individuals
## (I set mc.cores = 2 to comply with --as-cran checks, but you
##  should either use a reasonable number for your hardware or
##  leave it at its default value).

p1 <- oncoSimulPop(2, p701, mc.cores = 2)

par(mfrow = c(1, 2))
plot(p1, ask = FALSE)

## Stacked; we cannot log here, and harder to see patterns
plot(p1, ask = FALSE, type = "stacked")


## Show individual genotypes and drivers for an
## epistasis case with at most eight genotypes

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

## Drivers and clones
plot(e1, show = "drivers")

## Make genotypes explicit
plot(e1, show = "genotypes")

## Oh, but I want other colors
plot(e1, show = "genotypes", col = rainbow(8))

## and actually I want a line plot
plot(e1, show = "genotypes", type = "line")


}


\keyword{hplot}


