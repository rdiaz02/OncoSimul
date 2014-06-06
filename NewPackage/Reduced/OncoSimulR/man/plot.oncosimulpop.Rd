\name{plot.oncosimulpop}
\alias{plot.oncosimulpop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.oncosimulpop(x, ask = TRUE, col = c(8, "orange", 6:1), log = "y", ltyPop = 2:6, lwdPop = 0.2, ltyDrivers = 1, lwdDrivers = 3, xlab = "Time units", ylab = "Number of cells", plotClones = TRUE, plotDrivers = TRUE, addtot = FALSE, addtotlwd = 0.5, yl = NULL, thinPops = FALSE, thinPops.keep = 0.1, thinPops.min = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{ask}{
%%     ~~Describe \code{ask} here~~
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
function (x, ask = TRUE, col = c(8, "orange", 6:1), log = "y", 
    ltyPop = 2:6, lwdPop = 0.2, ltyDrivers = 1, lwdDrivers = 3, 
    xlab = "Time units", ylab = "Number of cells", plotClones = TRUE, 
    plotDrivers = TRUE, addtot = FALSE, addtotlwd = 0.5, yl = NULL, 
    thinPops = FALSE, thinPops.keep = 0.1, thinPops.min = 2, 
    ...) 
{
    op <- par(ask = ask)
    on.exit(par(op))
    null <- lapply(x, function(z) plot.oncosimul(z, col = col, 
        log = log, ltyPop = ltyPop, lwdPop = lwdPop, ltyDrivers = ltyDrivers, 
        lwdDrivers = lwdDrivers, xlab = xlab, ylab = ylab, plotClones = plotClones, 
        plotDrivers = plotDrivers, addtot = addtot, addtotlwd = addtotlwd, 
        yl = yl, thinPops = thinPops, thinPops.keep = thinPops.keep, 
        thinPops.min = thinPops.min, ...))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
