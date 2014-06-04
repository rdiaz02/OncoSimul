#!/bin/bash
rm ./Reduced/OncoSimulR/src/OncoSimulR.so
rm ./Reduced/OncoSimulR/src/OncoSimul.o
rm ./Reduced/OncoSimulR/src/symbols.rds
R-3.1.0 CMD INSTALL ./Reduced/OncoSimulR
