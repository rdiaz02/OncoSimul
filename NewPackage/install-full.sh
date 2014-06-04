#!/bin/bash
rm ./Full/OncoSimulR/src/OncoSimulR.so
rm ./Full/OncoSimulR/src/OncoSimul.o
rm ./Full/OncoSimulR/src/symbols.rds
R-3.1.0 CMD INSTALL ./Full/OncoSimulR
