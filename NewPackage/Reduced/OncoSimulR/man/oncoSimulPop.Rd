\name{oncoSimulPop}
\alias{oncoSimulPop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
oncoSimulPop(Nindiv, adjm, model = "Bozic", numGenes = 30, mu = 5e-07, detectionSize = 1e+08, detectionDrivers = 4, sampleEvery = 10, initSize = 2000, s = 0.1, sh = -1, K = 2000, keepEvery = sampleEvery, finalTime = 0.25 * 25 * 365, onlyCancer = TRUE, max.memory = 2000, max.wall.time = 200, endTimeEvery = -9, silent = FALSE, mc.cores = detectCores())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Nindiv}{
%%     ~~Describe \code{Nindiv} here~~
}
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
  \item{mc.cores}{
%%     ~~Describe \code{mc.cores} here~~
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
function (Nindiv, adjm, model = "Bozic", numGenes = 30, mu = 5e-07, 
    detectionSize = 1e+08, detectionDrivers = 4, sampleEvery = 10, 
    initSize = 2000, s = 0.1, sh = -1, K = 2000, keepEvery = sampleEvery, 
    finalTime = 0.25 * 25 * 365, onlyCancer = TRUE, max.memory = 2000, 
    max.wall.time = 200, endTimeEvery = -9, silent = FALSE, mc.cores = detectCores()) 
{
    if (.Platform$OS.type == "windows") {
        if (mc.cores != 1) 
            message("You are running Windows. Setting mc.cores = 1")
        mc.cores <- 1
    }
    pop <- mclapply(seq.int(Nindiv), function(x) oncoSimulIndiv(adjm = adjm, 
        model = model, numGenes = numGenes, mu = mu, detectionSize = detectionSize, 
        detectionDrivers = detectionDrivers, sampleEvery = sampleEvery, 
        initSize = initSize, s = s, sh = sh, K = K, keepEvery = keepEvery, 
        finalTime = finalTime, onlyCancer = onlyCancer, max.memory = max.memory, 
        max.wall.time = max.wall.time, endTimeEvery = endTimeEvery, 
        silent = silent), mc.cores = mc.cores)
    class(pop) <- "oncosimulpop"
    attributes(pop)$call <- match.call()
    return(pop)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
