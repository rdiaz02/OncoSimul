\name{evalGenotype}
\alias{evalGenotype}
\title{
Evaluate the fitness of a genotype.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
evalGenotype(genotype, fitnessEffects, verbose = FALSE, echo = FALSE, model = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genotype}{
%%     ~~Describe \code{genotype} here~~
}
  \item{fitnessEffects}{
%%     ~~Describe \code{fitnessEffects} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{echo}{
%%     ~~Describe \code{echo} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
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
function (genotype, fitnessEffects, verbose = FALSE, echo = FALSE, 
    model = "") 
{
    if (echo) 
        cat(paste("Genotype: ", genotype))
    if (!is.integer(genotype)) {
        if (length(grep(">", genotype))) {
            genotype <- nice.vector.eo(genotype, ">")
        }
        else if (length(grep(",", genotype))) {
            genotype <- nice.vector.eo(genotype, ",")
        }
        all.g.nums <- c(fitnessEffects$geneModule$GeneNumID, 
            fitnessEffects$long.geneNoInt$GeneNumID)
        all.g.names <- c(fitnessEffects$geneModule$Gene, fitnessEffects$long.geneNoInt$Gene)
        genotype <- all.g.nums[match(genotype, all.g.names)]
    }
    if (any(is.na(genotype))) 
        stop("genotype contains NAs or genes not in fitnessEffects")
    if (!length(genotype)) 
        stop("genotypes must have at least one mutated gene")
    if (any(genotype < 0)) {
        stop(paste("genotypes cannot contain negative values.", 
            "If you see this message, you found a bug."))
    }
    if (model \%in\% c("Bozic", "bozic1", "bozic2")) 
        prodNeg <- TRUE
    else prodNeg <- FALSE
    ff <- evalRGenotype(genotype, fitnessEffects, verbose, prodNeg)
    if (echo) {
        if (!prodNeg) 
            cat(" Fitness: ", ff, "\n")
        else cat(" Death rate: ", ff, "\n")
    }
    else {
        return(ff)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
