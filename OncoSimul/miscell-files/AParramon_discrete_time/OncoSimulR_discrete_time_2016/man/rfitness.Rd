\name{rfitness}
\alias{rfitness}


\title{Generate random fitness.}

\description{ Generate random fitness landscapes under a House of Cards,
  Rough Mount Fuji, or additive model.  }


\usage{

rfitness(g, c = 0.5, sd = 1, reference = "random", scale = NULL,
         wt_is_1 = TRUE, log = FALSE, min_accessible_genotypes = 0,
         accessible_th = 0)
}




\arguments{

  \item{g}{Number of genes.}

  \item{c}{The decrease in fitness of a genotype per each unit increase
    in Hamming distance from the reference genotype (see \code{reference}).}

  \item{sd}{The standard deviation of the random component (a normal
  distribution of mean 0 and standard deviation \code{sd}).}

\item{reference}{The reference genotype: for the deterministic, additive
  part, this is the genotype with maximal fitness, and all other
  genotypes decrease their fitness by \code{c} for every unit of Hamming
  distance from this reference. If "random" a genotype will be randomly
  chosen as the reference. If "max" the genotype with all positions
  mutated will be chosen as the reference. If you pass a vector (e.g.,
  \code{reference = c(1, 0, 1, 0)}) that will be the reference genotype.}

\item{scale}{Either NULL (nothing is done) or a two-element vector. If a
  two-element vector, fitness is re-scaled between \code{scale[1]} (the
  minimum) and \code{scale[2]} (the maximum).}

\item{wt_is_1}{If TRUE, fitness will be scaled so that the wildtype (the
  genotype with no mutations) has fitness of 1. This is applied after
  \code{scale}, so if you specify both it is most likely that the final
  fitness will not respect the limits in \code{scale}.}


\item{log}{If TRUE, log-transform fitness.}

\item{min_accessible_genotypes}{If larger than 0, the minimum number of
  accessible genotypes in the fitness landscape. A genotype is
  considered accessible if you can reach if from the wildtype by going
  through at least one path where all changes in fitness are larger or
  equal to \code{accessible_th}. The changes in fitness are considered
  at each mutational step, i.e., at each addition of one mutation we
  compute the difference between the genotype with \code{k + 1}
  mutations minus the ancestor genotype with \code{k} mutations. Thus, a
  genotype is considered accessible if there is at least one path where
  fitness increases at each mutational step by at least
  \code{accessible_th}.

  If the condition is not satisfied, we continue generating random
  fitness landscapes with the specified parameters until the condition
  is satisfied.
}

\item{accessible_th}{The threshold for the minimal change in fitness at
  each mutation step (i.e., between successive genotypes) that allows a
  genotype to be regarded as accessible. This only applies if
  \code{min_accessible_genotypes} is larger than 0.  So if you want to
  allow small decreases in fitness in successive steps, use a small
  negative value for \code{accessible_th}.  }

} 


\details{

  The model used here follows the Rough Mount Fuji model in Szendro et
  al., 2013 or Franke et al., 2011. Fitness is given as

  \deqn{f(i) = -c d(i, reference) + x_i}

  where \eqn{d(i, j)} is the Hamming distance between genotypes \eqn{i}
  and \eqn{j} (the number of positions that differ) and \eqn{x_i} is a
  random variable (in this case, a normal deviate of mean 0 and standard
  deviation \code{sd}).

  Setting \eqn{c = 0} we obtain a House of Cards model. Setting \eqn{sd
    = 0} fitness is given by the distance from the reference and if the
    reference is the genotype with all positions mutated, then we have a
    fully additive model (fitness increases linearly with the number of
    positions mutated).
} 

\value{
  
  An matrix with \code{g + 1} columns. Each column corresponds to a
  gene, except the last one that corresponds to fitness. 1/0 in a gene
  column denotes gene mutated/not-mutated. (For ease of use in other
  functions, this matrix has class  "genotype_fitness_matrix".) 

  If you have specified \code{min_accessible_genotypes > 0}, the return
  object has added attributes \code{accessible_genotypes} and
  \code{accessible_th} that show the number of accessible
  genotypes under the specified  threshold.
}
  
\references{

  Szendro I.~G. et al. (2013). Quantitative analyses of empirical
fitness landscapes. \emph{Journal of Statistical Mehcanics: Theory and
  Experiment\/}, \bold{01}, P01005.

Franke, J. et al. (2011). Evolutionary accessibility of mutational
pathways. \emph{PLoS Computational Biology\/}, \bold{7}(8), 1--9.

}

\author{ Ramon Diaz-Uriarte
}

\seealso{
  
  \code{\link{oncoSimulIndiv}},
  \code{\link{plot.genotype_fitness_matrix}},
  \code{\link{evalAllGenotypes}}
  \code{\link{allFitnessEffects}}
  \code{\link{plotFitnessLandscape}}

}
\examples{

## Random fitness for four genes-genotypes,
## plotting and simulating an oncogenetic trajectory

r1 <- rfitness(4)
plot(r1)
oncoSimulIndiv(allFitnessEffects(genotFitness = r1))

}

\keyword{ datagen }

