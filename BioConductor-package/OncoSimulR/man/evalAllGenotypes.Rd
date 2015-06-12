\name{evalAllGenotypes}
\alias{evalAllGenotypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate fitness of all possible genotypes.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
evalAllGenotypes(fitnessEffects, order = TRUE, max = 256, addwt = FALSE, model = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitnessEffects}{
%%     ~~Describe \code{fitnessEffects} here~~
}
  \item{order}{
%%     ~~Describe \code{order} here~~
}
  \item{max}{
%%     ~~Describe \code{max} here~~
}
  \item{addwt}{
%%     ~~Describe \code{addwt} here~~
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
function (fitnessEffects, order = TRUE, max = 256, addwt = FALSE, 
    model = "") 
{
    if (order) 
        tot <- function(n) {
            sum(sapply(seq.int(n), function(x) choose(n, x) * 
                factorial(x)))
        }
    else tot <- function(n) {
        2^n
    }
    nn <- nrow(fitnessEffects$geneModule) - 1 + nrow(fitnessEffects$long.geneNoInt)
    tnn <- tot(nn)
    if (tnn > max) {
        m <- paste("There are ", tnn, "genotypes.")
        m <- paste(m, " This is larger than max.")
        m <- paste(m, "Adjust max and rerun if you want")
        stop(m)
    }
    if (order) {
        f1 <- function(n) {
            lapply(seq.int(n), function(x) permutations(n = n, 
                r = x))
        }
    }
    else {
        f1 <- function(n) {
            lapply(seq.int(n), function(x) combinations(n = n, 
                r = x))
        }
    }
    genotNums <- f1(nn)
    list.of.vectors <- function(y) {
        lapply(unlist(lapply(y, function(x) {
            apply(x, 1, list)
        }), recursive = FALSE), function(m) m[[1]])
    }
    genotNums <- list.of.vectors(genotNums)
    names <- c(fitnessEffects$geneModule$Gene[-1], fitnessEffects$long.geneNoInt$Gene)
    genotNames <- unlist(lapply(lapply(genotNums, function(x) names[x]), 
        function(z) paste(z, collapse = if (order) {
            " > "
        }
        else {
            ", "
        })))
    if (model \%in\% c("Bozic", "bozic1", "bozic2")) 
        prodNeg <- TRUE
    else prodNeg <- FALSE
    allf <- vapply(genotNums, function(x) evalRGenotype(x, fitnessEffects, 
        FALSE, prodNeg), 1.1)
    df <- data.frame(Genotype = genotNames, Fitness = allf, stringsAsFactors = FALSE)
    if (addwt) 
        df <- rbind(data.frame(Genotype = "wt", Fitness = 1, 
            stringsAsFactors = FALSE), df)
    if (prodNeg) 
        colnames(df)[match("Fitness", colnames(df))] <- "Death_rate"
    return(df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
