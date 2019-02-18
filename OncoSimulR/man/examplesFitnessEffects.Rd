\name{examplesFitnessEffects}
\alias{examplesFitnessEffects}
\docType{data}
\title{
Examples of fitness effects
}

\description{

  Some examples of fitnessEffects objects. This is a collection, in a list,
  of most of the fitnessEffects created (using
  \code{\link{allFitnessEffects}}) for the vignette. See the vignette
  for descriptions and references.
  
  

}
\usage{data("examplesFitnessEffects")}
\format{
The format is a list of fitnessEffects objects.
}
%% \source{
%%   Gerstung et al., 2009. Quantifying cancer progression with conjunctive
%%   Bayesian networks. \emph{Bioinformatics}, 21: 2809--2815.
  
%%   Gerstung et al., 2011. The Temporal Order of Genetic and Pathway
%%   Alterations in Tumorigenesis. \emph{PLoS ONE}, 6.
  
%% }

\seealso{
  \code{\link{allFitnessEffects}}
}


\examples{
data(examplesFitnessEffects)

plot(examplesFitnessEffects[["fea"]])

evalAllGenotypes(examplesFitnessEffects[["cbn1"]], order = FALSE)
}
\keyword{datasets}


