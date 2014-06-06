\name{plot.oncosimul}
\alias{plot.oncosimul}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.oncosimul(x, col = c(8, "orange", 6:1), log = "y", ltyPop = 2:6, lwdPop = 0.2, ltyDrivers = 1, lwdDrivers = 3, xlab = "Time units", ylab = "Number of cells", plotClones = TRUE, plotDrivers = TRUE, addtot = FALSE, addtotlwd = 0.5, yl = NULL, thinPops = FALSE, thinPops.keep = 0.1, thinPops.min = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{log}{
%%     ~~Describe \code{log} here~~
}
  \item{ltyPop}{
%%     ~~Describe \code{ltyPop} here~~
}
  \item{lwdPop}{
%%     ~~Describe \code{lwdPop} here~~
}
  \item{ltyDrivers}{
%%     ~~Describe \code{ltyDrivers} here~~
}
  \item{lwdDrivers}{
%%     ~~Describe \code{lwdDrivers} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{plotClones}{
%%     ~~Describe \code{plotClones} here~~
}
  \item{plotDrivers}{
%%     ~~Describe \code{plotDrivers} here~~
}
  \item{addtot}{
%%     ~~Describe \code{addtot} here~~
}
  \item{addtotlwd}{
%%     ~~Describe \code{addtotlwd} here~~
}
  \item{yl}{
%%     ~~Describe \code{yl} here~~
}
  \item{thinPops}{
%%     ~~Describe \code{thinPops} here~~
}
  \item{thinPops.keep}{
%%     ~~Describe \code{thinPops.keep} here~~
}
  \item{thinPops.min}{
%%     ~~Describe \code{thinPops.min} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, col = c(8, "orange", 6:1), log = "y", ltyPop = 2:6, 
    lwdPop = 0.2, ltyDrivers = 1, lwdDrivers = 3, xlab = "Time units", 
    ylab = "Number of cells", plotClones = TRUE, plotDrivers = TRUE, 
    addtot = FALSE, addtotlwd = 0.5, yl = NULL, thinPops = FALSE, 
    thinPops.keep = 0.1, thinPops.min = 2, ...) 
{
    if (thinPops) 
        x <- thin.pops(x, keep = thinPops.keep, min.keep = thinPops.min)
    ndr <- apply(x$Genotypes[1:x$NumDrivers, , drop = FALSE], 
        2, sum)
    if (is.null(yl)) 
        yl <- c(1, max(apply(x$pops.by.time[, -1, drop = FALSE], 
            1, sum)))
    if (plotClones) {
        plotClones(x, ndr = ndr, xlab = xlab, ylab = ylab, lty = ltyPop, 
            col = col, ylim = yl, lwd = lwdPop, axes = FALSE, 
            log = log, ...)
    }
    if (plotClones && plotDrivers) 
        par(new = TRUE)
    if (plotDrivers) {
        plotDrivers0(x, timescale = 1, trim.no.drivers = FALSE, 
            xlab = "", ylab = "", lwd = lwdDrivers, lty = ltyDrivers, 
            col = col, addtot = addtot, addtotlwd = addtotlwd, 
            log = log, ylim = yl, ...)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
