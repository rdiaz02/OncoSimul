\name{POM}
\alias{POM}
\alias{LOD}
\alias{diversityPOM}
\alias{diversityLOD}
\alias{POM.oncosimul2}
\alias{LOD.oncosimul2}
\alias{POM.oncosimulpop}
\alias{LOD.oncosimulpop}


\title{
  Obtain Lines of Descent and Paths of the Maximum and their diversity from simulations.
}

\description{
  
  Compute Lines of Descent (LOD) and Path of the Maximum (POM) for a
  single simulation or a set of simulations (from \code{oncoSimulPop}).

  \code{diversityPOM} and \code{diversityLOD} return the Shannon's
  diversity (entropy) of the POM and LOD, respectively, of a set of
  simulations (it makes no sense to compute those from a single simulation).
  
}

\usage{

POM(x)
LOD(x)
diversityPOM(lpom)
diversityLOD(llod)
}

\arguments{ \item{x}{An object of class \code{oncosimulpop} (version >=
  2, so simulations with the old poset specification will not work) or
  class \code{oncosimul2} (a single simulation). }

%% \item{strict}{If TRUE, a single LOD as in Szendro et al. See Details.
%%   If FALSE, simulations must have been run with \code{keepPhylog = TRUE}
%%   to compute all possible LODs (see Details).}

\item{lpom}{A list of POMs, as returned from \code{POM} on an object of
  class \code{oncosimulpop}.}

\item{llod}{A list of LODs, as returned from \code{LOD} on an object of
  class \code{oncosimulpop}.}

% \item{...}{Other arguments passed to methods (ignored now).}
}

\details{

  Lines of Descent (LOD) and Path of the Maximum (POM) were defined in
  Szendro et al. (2013) and I follow those definitions here, as applied
  to a process in continuous time with sampling at user-specified
  periods.

  For POM, the results can depend strongly on how often we sample (i.e.,
  the \code{sampleEvery} argument to \code{oncoSimulIndiv} and
  \code{oncoSimulPop}), since the POM is computed by finding the clone
  with largest population size whenever we sample.%% from the values
  %% stored in the \code{pops.by.time} matrix.
  This also explains why
  it is generally meaningless to use POM on \code{oncoSimulSample} runs:
  these only keep the very last sample.


  For LOD, %% when using \code{strict = TRUE}, 
  a single LOD per simulation
  is returned, with the same meaning as that in p. 572 of Szendro et
  al. (2013). "A given genotype may undergo several episodes of colonization and extinction that are stored by the algorithm, and the last episode before the colonization of the final state is used to construct the step.",
  and I check that this genotype (which is the one that will become the
  most populated at final time) does not become extinct before the final
  colonization.

  %% If \code{strict = FALSE}, and if you have run the simulations with
  %% \code{keepPhylog = TRUE}, then a I return both \code{all_paths} and
  %% \code{lod_single}, with meanings as follow.  First, in case this might
  %% be useful, for each simulation I keep all the paths that
  %% "(...) arrive at the most populated genotype at the final time" (first
  %% paragraph in p. 572 of Szendro et al.), and these are stored in
  %% \code{all_paths}.  When \code{strict = FALSE} I also provide another
  %% single LOD for each run, too. This is the first path to arrive at the
  %% genotype that eventually becomes the most populated genotype at the
  %% final time (and, in this sense, agrees with the LOD of Szendro et
  %% al.). However, in contrast to what is done in Szendro
  %% ("A given genotype may undergo several episodes of colonization and extinction that are stored by the algorithm, and the last episode before the colonization of the final state is used to construct the step.")
  %% and when \code{strict = TRUE}, I do not check that this genotype
  %% (which is the one that will become the most populated at final time)
  %% does not become extinct before the final colonization. So there could
  %% be other paths (all in \code{all_paths}) that are actually the one(s)
  %% that are colonizers of the most populated genotype (with no extinction
  %% before the final colonization).

  Note \emph{breaking changes}: for LOD we used to return all lines of
  descent in a given simulation. In v. 2.9.1 we also returned the LOD
  as explained above. Now we only return the LOD as defined above.
  
  
}

\value{

  For \code{POM} either a character vector (if \code{x} is a single
  simulation) or a list of character vectors. Each character vector is
  the ordered set of genotypes that contain the largest subpopulation at
  the times of sampling.

  For \code{LOD}, if \code{x} is a single simulation, the line of
  descent as defined above (either an object of class "igraph.vs" (an
  igraph vertex sequence: see \code{\link[igraph]{vertex_attr}}) or a
  character vector if there were no descendants). If \code{x} is a list
  (population) of simulations, then a list where each element is a list
  as just explained.

  %% a two-element
  %% list. If \code{strict = TRUE}, only \code{lod_single} is returned. If
  %% \code{strict = FALSE} (and simulations were run with \code{keepPhylog
  %% = TRUE}), \code{all_paths} contains all paths to the maximum, and
  %% \code{lod_single} contains the single LOD which first arrives at the
  %% maximum.

  %% If \code{x} is a list (population) of simulations, then a list
  %% where each element is a two-element list, as just explained.
  %% All the lists
  %% contain objects of class "igraph.vs" (an igraph vertex sequence: see
  %% \code{\link[igraph]{vertex_attr}}).
 
  For \code{diversityLOD} and \code{diversityPOM} a single element
  vector with the Shannon's diversity (entropy) of the LODs (for
  \code{diversityLOD}) or of the POMs (for \code{diversityPOM}).

}

\references{

  Szendro, I. G., Franke, J., Visser, J. A. G. M. de, & Krug,
  J. (2013). Predictability of evolution depends nonmonotonically on
  population size. \emph{Proceedings of the National Academy of Sciences},
  110(2), 571-576. \url{https://doi.org/10.1073/pnas.1213613110}

}

\author{
  Ramon Diaz-Uriarte
}

\seealso{
  \code{\link{oncoSimulPop}}, \code{\link{oncoSimulIndiv}}
  
}

\examples{

######## Using a poset for pancreatic cancer from Gerstung et al.
###      (s and sh are made up for the example; only the structure
###       and names come from Gerstung et al.)

pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                          "TP53", "TP53", "MLL3"),
                                      child = c("KRAS","SMAD4", "CDNK2A", 
                                          "TP53", "MLL3",
                                          rep("PXDN", 3), rep("TGFBR2", 2)),
                                      s = 0.05,
                                      sh = -0.3,
                                      typeDep = "MN"))


pancr1 <- oncoSimulIndiv(pancr, model = "Exp")
pancr8 <- oncoSimulPop(8, pancr, model = "Exp",
                       mc.cores = 2)

POM(pancr1)
LOD(pancr1)

POM(pancr8)
LOD(pancr8)

diversityPOM(POM(pancr8))
diversityLOD(LOD(pancr8))



}

\keyword{manip}
\keyword{univar}









