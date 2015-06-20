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
or multiple realizations of the simulations, including whole-tumor and
single-cell sampling.


The /OncoSimulR directory contains the code for the [BioConductor](http://www.bioconductor.org) package
[OncoSimulR](http://www.bioconductor.org/packages/devel/bioc/html/OncoSimulR.html). The
/miscell-files directory contains additional files so far only related to
the above.


A former version of this code has been used in the paper "Identifying
restrictions in the order of accumulation of mutations during tumor
progression: effects of passengers, evolutionary models, and sampling",
[BMC Bioinformatics, 2015, 16:41](http://www.biomedcentral.com/1471-2105/16/41/abstract).




Installation
============

To use the latest code you first need to [use a development version of
Bioconductor](http://www.bioconductor.org/developers/how-to/useDevel/). Most
of the time, from R you only need to do:


```r
    library(BiocInstaller) 
    useDevel()
```

Then the next will install the development version of the package and its
dependencies, if needed:


```r
    source("http://bioconductor.org/biocLite.R")
    biocLite("OncoSimulR")
```

To start using it:

```r
library(OncoSimulR)
```



Documentation
=============

As any R/BioConductor package, OncoSimulR comes with documentation for its
user-visible functions and data sets (using the help is just standard R
usage). From BioConductor you can obtain the
[PDF reference manual](http://www.bioconductor.org/packages/3.2/bioc/manuals/OncoSimulR/man/OncoSimulR.pdf). A
better place to start, though, is the long vignette, with commented
examples (and created from the `OncoSimulR/vignettes/OncoSimulR.Rnw` file
that includes both text and code). Here is
[the vignette as PDF](http://www.bioconductor.org/packages/3.2/bioc/vignettes/OncoSimulR/inst/doc/OncoSimulR.pdf),
from the BioConductor site.


You can view the vignette from R itself doing


```r
browseVignettes("OncoSimulR")
```

and this gives you access to the PDF, the Rnw file (LaTeX + R), and the R code.


Licenses and copyright
======================

The R/BioConductor OncoSimulR package is licensed under the GPLv3
license. Except for the file `randutils.h` (see below), all of the code
for the OncoSimulR BioConductor package is Copyright 2012-2015 by Ramon
Diaz-Uriarte.


The code in `OncoSimulR/src/randutils.h` is copyright Melissa E. O'Neill,
and is licensed under "The MIT License (MIT)" in the terms explained in
the file itself. This is a license that is
[compatible with the GPL](http://directory.fsf.org/wiki/License:Expat).
The file randutils.hpp was downloaded from
(https://gist.github.com/imneme/540829265469e673d045) on 2015-06-20 and is
also referenced from the main article [Ease of Use without Loss of Power]
(http://www.pcg-random.org/posts/ease-of-use-without-loss-of-power.html). I
renamed it to randutils.h to conform to R's requirements (and changed the
`auto exit_func = hash(&_Exit);` line to keep R from complaining about the
Exit function).



The file under gitinfo-hooks is Copyright 2011 Brent Longborough, is
part of the
[gitinfo package](https://www.ctan.org/tex-archive/macros/latex/contrib/gitinfo?lang=en),
and is under the LaTeX Project Public License 1.3, which is
[incompatible with the GPL](http://directory.fsf.org/wiki/License:LPPLv1.3a). Note
this file is not part of the OncoSimulR BioConductor package.
