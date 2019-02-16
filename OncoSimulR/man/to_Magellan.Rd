\name{to_Magellan}
\alias{to_Magellan}

\title{ Create output for MAGELLAN.  }

\description{ Create a fitness landscape in a format that is understood
  by MAGELLAN \url{http://wwwabi.snv.jussieu.fr/public/Magellan/}. 

}

\usage{
to_Magellan(x, file,
            max_num_genotypes = 2000)
}

\arguments{
  
  \item{x}{ One of the following:
    \itemize{

      \item A matrix (or data frame) with g + 1 columns. Each of the
      first g columns contains a 1 or a 0 indicating that the gene of
      that column is mutated or not. Column g+ 1 contains the fitness
      values. This is, for instance, the output you will get from
      \code{\link{rfitness}}.

      \item A two column data frame. The second column is fitness, and
      the first column are genotypes, given as a character vector. For
      instance, a row "A, B" would mean the genotype with both A and B
      mutated.
      
      \item The output from a call to
      \code{\link{evalAllGenotypes}}. Make sure you use \code{order =
	FALSE} in that call.
      
      \item The output from a call to
      \code{\link{evalAllGenotypesMut}}. Make sure you use \code{order =
	FALSE}.
      
      \item The output from a call to \code{\link{allFitnessEffects}}
      (with no order effects in the specification).
      
    }
    
    The first two are the same as the format for the \code{genotFitness}
    component in \code{\link{allFitnessEffects}}.  }
  
  \item{file}{The name of the output file. If NULL, a name will be
  created using \code{\link{tempfile}}.}

  \item{max_num_genotypes}{Maximum allowed number of genotypes. For some
  types of input, we make a call to \code{\link{evalAllGenotypes}}, and
  use this as the maximum.}

}

\value{

  A file is written to disk. You can then plot and/or show summary
  statistics using MAGELLAN.
}

\note{
  If you try to pass a fitness specification with order effects you will
  receive an error, since that cannot be plotted with MAGELLAN.
}

\author{
Ramon Diaz-Uriarte
}

\references{
  MAGELLAN web site: \url{http://wwwabi.snv.jussieu.fr/public/Magellan/}

  Brouillet, S. et al. (2015). MAGELLAN: a tool to explore small fitness landscapes. \emph{bioRxiv}, \bold{31583}. \url{http://doi.org/10.1101/031583}

}



\seealso{
  \code{\link{allFitnessEffects}},
  \code{\link{evalAllGenotypes}},
  \code{\link{allFitnessEffects}},
  \code{\link{rfitness}}
}
\examples{

## Generate random fitness for four-genes genotype
## and export landscape.

r1 <- rfitness(4)
to_Magellan(r1, NULL)

## Specify fitness using a DAG and export it
cs <-  data.frame(parent = c(rep("Root", 3), "a", "d", "c"),
                      child = c("a", "b", "d", "e", "c", "f"),
                      s = 0.1,
                      sh = -0.9,
                      typeDep = "MN")

to_Magellan(allFitnessEffects(cs), NULL)

}

\keyword{ manip }

