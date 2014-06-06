\name{oncoSimulIndiv}
\alias{oncoSimulIndiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
oncoSimulIndiv(adjm, model = "Bozic", numGenes = 30, mu = 5e-07, detectionSize = 1e+06, detectionDrivers = 4, sampleEvery = 15, initSize = 2000, s = 0.1, sh = -1, K = initSize/(exp(1) - 1), keepEvery = sampleEvery, finalTime = 0.25 * 25 * 365, onlyCancer = TRUE, max.memory = 2000, max.wall.time = 200, endTimeEvery = -9, silent = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjm}{
%%     ~~Describe \code{adjm} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{numGenes}{
%%     ~~Describe \code{numGenes} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{detectionSize}{
%%     ~~Describe \code{detectionSize} here~~
}
  \item{detectionDrivers}{
%%     ~~Describe \code{detectionDrivers} here~~
}
  \item{sampleEvery}{
%%     ~~Describe \code{sampleEvery} here~~
}
  \item{initSize}{
%%     ~~Describe \code{initSize} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
  \item{sh}{
%%     ~~Describe \code{sh} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{keepEvery}{
%%     ~~Describe \code{keepEvery} here~~
}
  \item{finalTime}{
%%     ~~Describe \code{finalTime} here~~
}
  \item{onlyCancer}{
%%     ~~Describe \code{onlyCancer} here~~
}
  \item{max.memory}{
%%     ~~Describe \code{max.memory} here~~
}
  \item{max.wall.time}{
%%     ~~Describe \code{max.wall.time} here~~
}
  \item{endTimeEvery}{
%%     ~~Describe \code{endTimeEvery} here~~
}
  \item{silent}{
%%     ~~Describe \code{silent} here~~
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
function (adjm, model = "Bozic", numGenes = 30, mu = 5e-07, detectionSize = 1e+06, 
    detectionDrivers = 4, sampleEvery = 15, initSize = 2000, 
    s = 0.1, sh = -1, K = initSize/(exp(1) - 1), keepEvery = sampleEvery, 
    finalTime = 0.25 * 25 * 365, onlyCancer = TRUE, max.memory = 2000, 
    max.wall.time = 200, endTimeEvery = -9, silent = FALSE) 
{
    rt <- adjmat.to.restrictTable(adjm)
    call <- match.call()
    typeFitness <- switch(model, Bozic = "bozic1", Exp = "exp", 
        McFarland = "mcfarlandlog", McF = "mcfarlandlog", stop("No valid value for model"))
    if (typeFitness == "exp") {
        death <- 1
        mutatorGenotype <- 1
    }
    else {
        death <- -99
        mutatorGenotype <- 0
    }
    birth <- -99
    if ((typeFitness == "mcfarlandlog") && (sampleEvery > 0.05)) {
        warning("With the McFarland model you generally want smaller sampleEvery")
    }
    doneSimuls <- FALSE
    while (!doneSimuls) {
        op <- try(oncoSimul.internal(restrict.table = rt, numGenes = numGenes, 
            typeFitness = typeFitness, typeCBN = "-99", birth = birth, 
            s = s, sh = sh, death = death, mu = mu, initSize = initSize, 
            sampleEvery = sampleEvery, detectionSize = detectionSize, 
            mutatorGenotype = mutatorGenotype, finalTime = finalTime, 
            initSize_species = 2000, initSize_iter = 500, seed_gsl = NULL, 
            verbosity = 1, initMutant = -1, speciesFS = 40000, 
            ratioForce = 2, max.memory = max.memory, max.wall.time = max.wall.time, 
            keepEvery = keepEvery, alpha = 0.0015, K = K, endTimeEvery = endTimeEvery, 
            finalDrivers = detectionDrivers, silent = silent), 
            silent = silent)
        if (!inherits(op, "try-error")) {
            if (!silent) {
                cat("\n ... finished this run:")
                cat("\n       Total Pop Size = ", op$TotalPopSize)
                cat("\n       Drivers Last = ", op$MaxDriversLast)
                cat("\n       Final Time = ", op$FinalTime)
            }
            if (onlyCancer) {
                doneSimuls <- reachCancer(op, ndr = detectionDrivers, 
                  detectionSize = detectionSize, maxPopSize = 1e+15)
            }
            else {
                doneSimuls <- TRUE
            }
            if (!silent) {
                if (doneSimuls) 
                  cat("\n ... Keeping this one\n")
                else cat("\n ... Cancer not reached\n")
            }
        }
        else {
            if (!silent) 
                cat("\n Simulation aborted because of numerical problems.", 
                  "Proceeding to next one.\n")
        }
    }
    class(op) <- "oncosimul"
    attributes(op)$call <- call
    return(op)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
