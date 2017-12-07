\name{samplePop}
\alias{samplePop}
\alias{sampledGenotypes}
\alias{print.sampledGenotypes}
\title{
  Obtain a sample from a population of simulations.
  
}
\description{
  
  Obtain a sample (a matrix of individuals/samples by genes or, equivalently, a
  vector of "genotypes") from an oncosimulpop object (i.e., a simulation
  of multiple individuals) or a single oncosimul object. Sampling
  schemes include whole tumor and single cell sampling, and sampling at
  the end of the tumor progression or during the progression of the
  disease.

  
  \code{sampledGenotypes} shows the genotype frequencies from that
  sample; Shannon's diversity ---entropy--- of the genotypes is also
  returned.  Order effects are ignored.
  
}

\usage{
samplePop(x, timeSample = "last", typeSample = "whole",
          thresholdWhole = 0.5, geneNames = NULL, popSizeSample = NULL,
          propError = 0)

sampledGenotypes(y, genes = NULL)

}

\arguments{
  \item{x}{An object of class \code{oncosimulpop} or class \code{oncosimul2} (a
  single simulation).}

  \item{y}{The output from a call to \code{samplePop}.}
  
  \item{timeSample}{
    "last" means to sample each individual in the very last time period
    of the simulation. "unif" (or "uniform") means sampling each
    individual at a time choosen uniformly from all the times recorded
    in the simulation between the time when the first driver appeared
    and the final time period. "unif" means that it is almost sure that
    different individuals will be sampled at different times. "last"
    does not guarantee that different individuals will be sampled at the
    same time unit, only that all will be sampled in the last time unit
    of their simulation.

    You can, alternatively, specify the population size at which you
    want the sample to be taken. See argument \code{popSizeSample}.
  } 
  
  \item{typeSample}{
    "singleCell" (or "single") for single cell sampling, where the
    probability of sampling a cell (a clone) is directly proportional to
    its population size.  "wholeTumor" (or "whole") for whole tumor
    sampling (i.e., this is similar to a biopsy being the entire tumor).
    "singleCell-noWT" or "single-nowt" is single cell sampling, but
    excluding the wild type.
  }
  
  \item{thresholdWhole}{
    In whole tumor sampling, whether a gene is detected as mutated depends
    on thresholdWhole: a gene is considered mutated if it is altered in at
    least thresholdWhole proportion of the cells in that individual.
  }

  \item{geneNames}{An optional vector of gene names so as to label the
    column names of the output.}

  \item{popSizeSample}{An optional vector of total population sizes at
    which you want the samples to be taken. If you pass this vector,
    \code{timeSample} has no effect. The samples will be taken at the
    first time at which the population size gets as large as (or larger
    than) the size specified in \code{popSizeSample}.

    This allows you to specify arbitrary sampling schemes with respect
    to total population size.  }

  \item{propError}{The proportion of observations with error (for
    instance, genotyping error). If larger than 0, this proportion of
    entries in the sampled matrix will be flipped (i.e., 0s turned
    to 1s and 1s turned to 0s).}

  
  \item{genes}{If non-NULL, use only the genes in \code{genes} to create
  the table of genotypes. This can be useful if you only care about the
  genotypes with respect to a subset of genes (say, X), and want to
  collapse with respect to another subset of genes (say, Y), for
  instance if Y is a large set of passenger genes. For example, suppose
  the complete set of genes is 'a', 'b', 'c', 'd', and you specify
  \code{genes = c('a', 'b')}; then, genotypes 'a, b, c' and genotypes
  'a, b, d' will not be shown as different rows in the table of
  frequencies. Likewise, genotypes 'a, c' and genotypes 'a, d' will not
  be shown as different rows.  Of course, if what are actually different
  genotypes are not regarded as different, this will affect the
  calculation of the diversity.}

}

\details{
  samplePop simply repeats the sampling process in each individual of
  the oncosimulpop object.

  Please see \code{\link{oncoSimulSample}} for a much more efficient way
  of sampling when you are sure what you want to sample.

  Note that if you have set \code{onlyCancer = FALSE} in the call to
  \code{\link{oncoSimulSample}}, you can end up trying to sample from
  simulations where the population size is 0. In this case, you will get
  a vector/matrix of NAs and a warning.

  Similarly, when using \code{timeSample = "last"} you might end up with
  a vector of 0 (not NAs) because you are sampling from a population
  that contains no clones with mutated genes. This event (sampling from
  a population that contains no clones with mutated genes), by
  construction, cannot happen when \code{timeSample = "unif"} as
  "uniform" sampling is taken here to mean sampling at a time choosen
  uniformly from all the times recorded in the simulation between the
  time when the first driver appeared and the final time
  period. However, you might still get a vector of 0, with uniform
  sampling, if you sample from a population that contains only a few
  cells with any mutated genes, and most cells with no mutated genes.}

\value{

  A matrix. Each row is a "sample genotype", where 0 denotes no alteration
  and 1 alteration. When using v.2, columns are named with the gene
  names.

  We quote "sample genotype" because when not using single cell, a row
  (a sample genotype) need not be, of course, any really existing
  genotype in a population as we are genotyping a whole tumor. Suppose
  there are really two genotypes present in the population, genotype A,
  which has gene A mutated and genotype B, which has gene B
  mutated. Genotype A has a frequency of 60\% (so B's frequency is
  40\%). If you use whole tumor sampling with \code{thresholdWhole =
  0.4} you will obtain a genotype with A and B mutated.
  %% To make it more clear that genes/laterations are in
  %% columns, columns are named starting by "G."  (for "gene").


  For \code{sampledGenotypes} a data frame with two columns: genotypes
  and frequencies. This data frame has an additional attribute,
  "ShannonI", where Shannon's index of diversity (entropy) is
  stored. This is an object of class  "sampledGenotypes"
  with an S3 print method.

}

\references{
  Diaz-Uriarte, R. (2015). Identifying restrictions in the order of
  accumulation of mutations during tumor progression: effects of
  passengers, evolutionary models, and sampling
  \url{http://www.biomedcentral.com/1471-2105/16/41/abstract}

}

\author{
  Ramon Diaz-Uriarte
}

\seealso{
  \code{\link{oncoSimulPop}}, \code{\link{oncoSimulSample}}
  
}

\examples{
data(examplePosets)
p705 <- examplePosets[["p705"]]

## (I set mc.cores = 2 to comply with --as-cran checks, but you
##  should either use a reasonable number for your hardware or
##  leave it at its default value).

p1 <- oncoSimulPop(4, p705, mc.cores = 2)
(sp1 <- samplePop(p1))
sampledGenotypes(sp1)


## Sample at fixed sizes. Notice the requested size
## for the last population is larger than the any population size
## so we get NAs

(sp2 <- samplePop(p1, popSizeSample = c(1e7, 1e6, 4e5, 1e13)))
sampledGenotypes(sp2)


## Now single cell sampling

r1 <- oncoSimulIndiv(p705)
samplePop(r1, typeSample = "single")

sampledGenotypes(samplePop(r1, typeSample = "single"))

}

\keyword{manip}









