<!-- [![Travis-CI Build Status](https://travis-ci.org/rdiaz02/OncoSimul.svg?branch=master)](https://travis-ci.org/rdiaz02/OncoSimul) -->
<!-- [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/rdiaz02/OncoSimul?branch=master&svg=true)](https://ci.appveyor.com/project/rdiaz02/OncoSimul) -->
<!-- [![codecov.io](https://codecov.io/github/rdiaz02/OncoSimul/coverage.svg?branch=master)](https://codecov.io/github/rdiaz02/OncoSimul?branch=master) -->



# OncoSimul



Code for forward population genetic simulation in asexual populations,
with special focus on cancer progression.  Fitness can be an arbitrary
function of genetic interactions between multiple genes or modules of
genes, including epistasis, order restrictions in mutation accumulation,
and order effects.  Mutation rates can differ between genes, and we can
include mutator/antimutator genes (to model mutator
phenotypes). Simulations so far use continuous-time models and can include
driver and passenger genes and modules. Also included are functions for:
simulating random DAGs of the type found in Oncogenetic Tress, Conjunctive
Bayesian Networks, and other tumor progression models; plotting and
sampling from single or multiple realizations of the simulations,
including single-cell sampling; plotting the parent-child relationships of
the clones; generating random fitness landscapes (Rough Mount Fuji, House
of Cards, and additive models) and plotting them.



The /OncoSimulR directory contains the code for the [BioConductor](http://www.bioconductor.org) package
[OncoSimulR](http://www.bioconductor.org/packages/devel/bioc/html/OncoSimulR.html). The
/miscell-files directory contains additional files so far only related to
the above.


A former version of this code has been used in the paper "Identifying
restrictions in the order of accumulation of mutations during tumor
progression: effects of passengers, evolutionary models, and
sampling",
[BMC Bioinformatics, 2015, 16:41](http://www.biomedcentral.com/1471-2105/16/41).
OncoSimulR has also been used extensively in the simulations reported in
the [bioRxiv](biorxiv.org/)
preprint
["Cancer Progression Models And Fitness Landscapes: A Many-To-Many Relationship"](https://doi.org/10.1101/141465).


You can also find
[OncoSimulR](https://popmodels.cancercontrol.cancer.gov/gsr/packages/oncosimulr/)
on the [Genetic Simulation Resources
catalogue](https://popmodels.cancercontrol.cancer.gov/gsr/).
<a href="http://popmodels.cancercontrol.cancer.gov/gsr/"><img src="http://popmodels.cancercontrol.cancer.gov/gsr/static/img/gsr_tile.jpg" alt="Catalogued on GSR" width="180" height="60" /></a>

# Installation


To use the most recent code that I regard as stable you first need to [use
a development version of Bioconductor](http://www.bioconductor.org/developers/how-to/useDevel/). Most
of the time, from R you only need to do:


```r
    library(BiocInstaller) 
    useDevel()
```

Then the next code will install the development version of the package and
its dependencies, if needed:


```r
    source("http://bioconductor.org/biocLite.R")
    biocLite("OncoSimulR")
```

To start using it:

```r
library(OncoSimulR)
```



You can directly install from github:

```r
install.packages("devtools") ## if you don't have it already
library(devtools)
install_github("rdiaz02/OncoSimul/OncoSimulR")
``` 

But sometimes the latest additions in this repo could be broken (see [Software status](#software-status)). And you
can of course clone this repo, and install from there.

## BioConductor github repository

The github repository for this package is this one:
https://github.com/rdiaz02/OncoSimul . Since mid-2017 BioConductor is
maintained using git, but since this directory contains other files and
directories in addition to the OncoSimulR package itself, I have not used
option ["Sync an existing GitHub repository with
Bioconductor"](https://www.bioconductor.org/developers/how-to/git/sync-existing-repositories). Instead,
I continue using this github repo, but then locally update a
Bioconductor-only repository of just the OncoSimulR code (as explained in
[Maintain a Bioconductor-only repository for an existing
package](https://www.bioconductor.org/developers/how-to/git/maintain-bioc-only/)).

# Documentation


As any R/BioConductor package, OncoSimulR comes with documentation for its
user-visible functions and data sets (using the help is just standard R
usage). From
[OncoSimulR's BioConductor page](https://www.bioconductor.org/packages/devel/bioc/html/OncoSimulR.html)
you have access to the standard documentation both the manual and overview
---the vignette. The best place to start is the vignette (created from the
`OncoSimulR/vignettes/OncoSimulR.Rnw` file that includes both text and
code).

<!-- you can obtain the  -->
<!-- [PDF reference manual](http://www.bioconductor.org/packages/3.2/bioc/manuals/OncoSimulR/man/OncoSimulR.pdf). A -->
<!-- better place to start, though, is the long vignette, with commented -->
<!-- examples (and created from the `OncoSimulR/vignettes/OncoSimulR.Rnw` file -->
<!-- that includes both text and code). Here is -->
<!-- [the vignette as PDF](http://www.bioconductor.org/packages/devel/bioc/vignettes/OncoSimulR/inst/doc/OncoSimulR.pdf), -->
<!-- from the BioConductor site (the development branch). -->


You can view the vignette from R itself doing


```r
browseVignettes("OncoSimulR")
```

and this gives you access to the HTML, the Rmd file (markdown + R), and the R code.

## Documentation: HTML and PDF for this repo's version


From these two links you can also
[browse the HTML vignette](https://rdiaz02.github.io/OncoSimul/OncoSimulR.html)
and [get a PDF version](https://rdiaz02.github.io/OncoSimul/pdfs/OncoSimulR.pdf).

These files correspond to the most recent, github version, of the package
(i.e., they might include changes not yet available from the BioConudctor
package).


## Further documentation

This [paper published in Bioinformatics](https://doi.org/10.1093/bioinformatics/btx077)
gives a quick overview of OncoSimulR (a former version is available as a 
[bioRxiv preprint](http://biorxiv.org/content/early/2016/08/14/069500)). You can also take a look at this
[poster presented at ECCB 2016](http://dx.doi.org/10.7490/f1000research.1112860.1).


If you use the package in publications **please cite the Bioinformatics paper**.



# Licenses and copyright


The R/BioConductor OncoSimulR package is licensed under the GPLv3
license. All of the code for the OncoSimulR BioConductor package, except
for functions `plot.stream` and `plot.stacked`, is Copyright 2012-2016 by
Ramon Diaz-Uriarte. `plot.stream` and `plot.stacked` are Copyright
2013-2016 by Marc Taylor (see also https://github.com/marchtaylor/sinkr
and
http://menugget.blogspot.com.es/2013/12/data-mountains-and-streams-stacked-area.html).


The code in `miscell-files/randutils.h` is copyright Melissa E. O'Neill,
and is licensed under "The MIT License (MIT)" in the terms explained in
the file itself. This is a license that is
[compatible with the GPL](http://directory.fsf.org/wiki/License:Expat).
The file randutils.hpp was downloaded from
https://gist.github.com/imneme/540829265469e673d045 on 2015-06-20 and is
also referenced from the main article [Ease of Use without Loss of Power]
(http://www.pcg-random.org/posts/ease-of-use-without-loss-of-power.html). I
renamed it to randutils.h to conform to R's requirements (and changed the
`auto exit_func = hash(&_Exit);` line to keep R from complaining about the
Exit function). I had to disable usage of randutils for now, since I could
not get it to compile with gcc-4.6 (since version 3.3 of R, 
the official Rtools for Windows now support C++-11, so I might change 
this in the near future).



The file under gitinfo-hooks is Copyright 2011 Brent Longborough, is
part of the
[gitinfo package](https://www.ctan.org/tex-archive/macros/latex/contrib/gitinfo?lang=en),
and is under the LaTeX Project Public License 1.3, which is
[incompatible with the GPL](http://directory.fsf.org/wiki/License:LPPLv1.3a). Note
this file is not part of the OncoSimulR BioConductor package.


The files under miscell-files/AParramon_discrete_time are copyright
Alberto Parramon, unless otherwise specified. This is an implementation of
a discrete-time version of OncoSimulR.


# Software status


|             | Bioconductor (multiple platforms)   | Travis CI  (Linux)  | Appveyor (Windows)  |
| ------------- | ------------------- | ------------- | ---------------- |
| R CMD check   | <a href="http://bioconductor.org/checkResults/release/bioc-LATEST/OncoSimulR/"><img border="0" src="http://bioconductor.org/shields/build/release/bioc/OncoSimulR.svg" alt="Build status"></a> (release)</br><a href="http://bioconductor.org/checkResults/devel/bioc-LATEST/OncoSimulR/"><img border="0" src="http://bioconductor.org/shields/build/devel/bioc/OncoSimulR.svg" alt="Build status"></a> (devel) | <a href="https://travis-ci.org/rdiaz02/OncoSimul"><img src="https://travis-ci.org/rdiaz02/OncoSimul.svg?branch=master" alt="Build status"></a> | <a href="https://ci.appveyor.com/project/rdiaz02/OncoSimul"><img src="https://ci.appveyor.com/api/projects/status/github/rdiaz02/OncoSimul?branch=master&svg=true" alt="Build status"></a> |
| Test coverage |                     | <a href="https://codecov.io/github/rdiaz02/OncoSimul?branch=master"><img src="https://codecov.io/github/rdiaz02/OncoSimul/coverage.svg?branch=master" alt="Coverage Status"/></a>   |                  |

(Note: Some of Appveyor's failures, if any, should be taken with a grain
of salt, as Appveyor can fail for reasons that have nothing to do with
the package, such as R not being downloaded correctly, etc. Look at the
details of each failure. Similarly, some of the errors in BioConductor,
specially in the development branch, can be caused, specially in Windows,
by some required packages not being yet available, often "car" and
_"igraph". Again, look at the details of each failure.)
<!-- Based on https://raw.githubusercontent.com/Bioconductor-mirror/illuminaio/master/README.md -->

<!-- [![Travis-CI Build Status](https://travis-ci.org/rdiaz02/OncoSimul.svg?branch=master)](https://travis-ci.org/rdiaz02/OncoSimul) -->
<!-- [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/rdiaz02/OncoSimul?branch=master&svg=true)](https://ci.appveyor.com/project/rdiaz02/OncoSimul) -->
<!-- [![codecov.io](https://codecov.io/github/rdiaz02/OncoSimul/coverage.svg?branch=master)](https://codecov.io/github/rdiaz02/OncoSimul?branch=master) -->

