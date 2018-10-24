#!/bin/bash
## pkill -f bioc-devel
R-3.2.0-bioc-devel CMD INSTALL OncoSimulR_1.99.1.tar.gz
R-3.2.0-bioc-devel --vanilla <  small-benchmark.R > small-benchmark-99.1.Rout

R-3.2.0-bioc-devel CMD INSTALL OncoSimulR_1.99.2.tar.gz
R-3.2.0-bioc-devel --vanilla <  small-benchmark.R > small-benchmark-99.2.Rout


