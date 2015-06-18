OncoSimul
=========


Code for forward population genetic simulation in asexual populations,
with special focus on cancer progression.  Fitness can be an arbitrary
function of genetic interactions between multiple genes, including
epistasis, order restrictions in mutation accumulation, and order effects.
Simulations so far use continuous-time models and can include driver and
passenger genes.  Also included are functions for simulating random DAGs
of the type found in Oncogenetic Tress, Conjunctive Bayesian Networks, and
other tumor progression models, and for plotting and sampling from single
or multiple realizations of the simulations.



The /OncoSimulR directory contains the code for the [BioConductor](http://www.bioconductor.org) package
[OncoSimulR](http://www.bioconductor.org/packages/devel/bioc/html/OncoSimulR.html). The
/miscell-files directory contains additional files so far only related to
the above.

