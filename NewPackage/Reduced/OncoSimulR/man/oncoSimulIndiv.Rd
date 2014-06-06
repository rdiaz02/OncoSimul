\name{oncoSimulIndiv}
\alias{oncoSimulIndiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate tumor progression.
}
\description{
  
  Simulate tumor progression including possible restrictions in the
  order of driver mutations. Optionally add passenger
  mutations. Simulation is done using the BNB algorithm of Mather et
  al., 2012.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{
 oncoSimulIndiv(adjm, model = "Bozic", numGenes = 30, mu = 5e-07, detectionSize = 1e+06, detectionDrivers = 4, sampleEvery = 15, initSize = 2000, s = 0.1, sh = -1, K = initSize/(exp(1) - 1), keepEvery = sampleEvery, finalTime = 0.25 * 25 * 365, onlyCancer = TRUE, max.memory = 2000, max.wall.time = 200, silent = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjm}{
    
    The adjacency matrix that specifies the order restrictions. See \code{\link{poset2AdjMat}}.
%%     ~~Describe \code{adjm} here~~
}
\item{model}{
  One of "Bozic", "Exp", "McFarlandLog" (the last one can be abbreviated
  to "McFL").
%%     ~~Describe \code{model} here~~
}
\item{numGenes}{
  The total number of genes tracked. Must be smaller than 64. This
  includes both driver and passenger genes.

  All driver genes should be included in the adjacency matrix of
  restrictions (even if they depend on no one and no one depends on
  them), and will be numbered from 1 to the total number of driver
  genes. Thus, passenger genes will be numbered from (number of driver
  genes + 1):numGenes.
%%     ~~Describe \code{numGenes} here~~
}
\item{mu}{
  Mutation rate.
%%     ~~Describe \code{mu} here~~
}
\item{detectionSize}{
  What is the minimal number of cells for cancer to be detected.
%%     ~~Describe \code{detectionSize} here~~
}
\item{detectionDrivers}{
  The minimal number of drivers present in any clone for cancer to be detected.
%%     ~~Describe \code{detectionDrivers} here~~
}
\item{sampleEvery}{
  How often the whole population is sampled.
%%     ~~Describe \code{sampleEvery} here~~
}
\item{initSize}{
  Initial population size.
%%     ~~Describe \code{initSize} here~~
}
\item{s}{
  Selection coefficient for drivers.
%%     ~~Describe \code{s} here~~
}
\item{sh}{
  Selection coefficient for drivers with restrictions not satisfied. A
  value of 0 means there are no penalties for a driver appearing in a
  clone when its restrictions are not satisfied.

  To specify "sh=Inf" (in Diaz-Uriarte, 2014) use sh = -1.
%%     ~~Describe \code{sh} here~~
}
\item{K}{
  Initial population equilibrium size in the McFarland models.
%%     ~~Describe \code{K} here~~
}
\item{keepEvery}{
  Time interval between successive whole population samples that are
  actually stored. This must be larger or equal to sampleEvery. If keepEvery is
  not a multiple integer of sampleEvery, the keepEvery in use will be the
  smallest multiple integer of keepEvery larger than the specified keepEvery.
%%     ~~Describe \code{keepEvery} here~~
}
\item{finalTime}{
  What is the maximum number of time units that the simulation can run.
%%     ~~Describe \code{finalTime} here~~
}
\item{onlyCancer}{
  Return only simulations that reach cancer?

  If set to TRUE, only simulations that satisfy the detectionDrivers or
  the detectionSize will be returned (i.e., the simulation will be
  repeated until one which meets that cancer requirements is
  met). Otherwise, the simulation is returned regardless of final
  population size or number of drivers in any clone and this includes
  simulations where the population goes extinct.
%%     ~~Describe \code{onlyCancer} here~~
}
\item{max.memory}{
  The largest size (in MB) of the matrix of Populations by Time. If it
  creating it would use more than this amount of memory, it is not
  created. This prevents you from accidentally passing parameters that
  will return an enormous object.
%%     ~~Describe \code{max.memory} here~~
}
\item{max.wall.time}{
  Maximum wall time for each individual simulation run. If the
  simulation is not done in this time, it is aborted.
%%     ~~Describe \code{max.wall.time} here~~
}
%%   \item{endTimeEvery}{
%% %%     ~~Describe \code{endTimeEvery} here~~
%% }
\item{silent}{
  If FALSE, intermediate information about what is being done is
  provided during the simulation, including possible warnings about
  numerical problems.
%%     ~~Describe \code{silent} here~~
}
}
\details{

  The basic simulation algorithm implemented is the BNB one of Mather et
  al., 2012, where I have added modifications to fitness based on the
  restrictions in the order of mutations.

  Full details about the algorithm are provided in Mather et al.,
  2012. The evolutionary models, including references, and the rest of
  the parameters are explained in Diaz-Uriarte, 2014, especially in the
  Supplementary Material. The model called "Bozic" is based on Bozic et
  al., 2010, and the model called "McFarland" in McFarland et al., 2013.

  
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A list, of class "oncosimul", with the following components:
  
  \item{pops.by.time }{A matrix of the population sizes of the clones,
    with clones in columns and time in row. Not all clones are shown here,
    only those that were present in at least on of the keepEvery samples.}
  
  \item{NumClones   }{Total number of clones in the above matrix.  }

  \item{TotalPopSize}{Total population size at the end.}

  \item{Genotypes}{A matrix of genotypes. For each of the clones in the
    pops.by.time matrix, its genotype, with a 0 if the gene is not
    mutated and a 1 if it is mutated.}
  
  \item{MaxNumDrivers}{The largest number of mutated driver genes ever seen in the
    simulation in any clone.}

  \item{MaxDriversLast}{The largest number of mutated drivers in any clone at the
    end of the simulation.}

  \item{NumDriversLargestPop}{The number of mutated driver genes in the
    clone with largest population size. }

  \item{LargestClone}{Population size of the clone with largest number
    of population size.}

  \item{PropLargestPopLast}{Ratio of LargestClone/TotalPopSize}

  \item{FinalTime}{The time (in time units) at the end of the
    simulation.}

  \item{NumIter}{The number of iterations of the BNB algorithm.}

  \item{HittedWallTime}{TRUE if we reached the limit of max.wall.time. FALSE
    otherwise.}

  \item{TotalPresentDrivers}{The total number of mutated driver genes,
    whether or not in the same clone. The number of elements in
    \code{OccurringDrivers}, below.}

  \item{CountByDriver}{A vector of length number of drivers, with the
    count of the number of clones that have that driver mutated.}

  \item{OccurringDrivers}{The actual number of drivers mutated.}

  \item{PerSampleStats}{A 5 column matrix with a row for each sampling
    period. The columns are: total population size, population size of the
    largest clone, the ratio of the two, the largest number of drivers in
    any clone, and the number of drivers in the clone with the largest
    population size.}
  
  \item{other}{A list that contains statistics for an estimate of the
    simulation error when using the McFarland model. The relevant value
    is errorMF, which is -99 unless in the McFarland model. For the
    McFarland model it is the largest difference of successive death
    rates.}
  

  
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{

  Bozic, I., et al.,  (2010). Accumulation of driver and passenger mutations during tumor
  progression. {\em Proceedings of the National Academy of Sciences of the United
  States of America\/}, {\bf 107}, 18545--18550.

  Diaz-Uriarte, R. (2014). Inferring restrictions in the temporal order
of mutations during tumor progression: effects of passenger mutations,
evolutionary models, and sampling. \url{http://dx.doi.org/10.1101/005587}

  McFarland, C.~D. et al. (2013). Impact of deleterious passenger
  mutations on cancer progression.  {\em Proceedings of the National
  Academy of Sciences of the United States of America\/}, {\bf 110}(8),
  2910--5.

  Mather, W.~H., Hasty, J., and Tsimring, L.~S. (2012). Fast stochastic
  algorithm for simulating evolutionary population dynamics.  {\em
    Bioinformatics (Oxford, England)\/}, {\bf 28}(9), 1230--1238.

  
%% ~put references to the literature/web site here ~
}
\author{
  Ramon Diaz-Uriarte
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{plot.oncosimul}}, \code{\link{print.oncosimul}},
  \code{\link{summary.oncosimul}},
  \code{\link{oncoSimulPop}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
\keyword{iteration}% __ONLY ONE__ keyword per line
