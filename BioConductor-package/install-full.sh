#!/bin/bash

V_R=$1
export R_MAKEVARS_USER=/home/ramon/.R/Makevars



# export OPTIMFLAGS=-march=native -O2 -ffunction-sections -g -fpic -pipe -Wall -pedantic
# export CFLAGS=-Wall -pedantic -march=native -O2 -ffunction-sections -g -fpic -pipe 
# export CXXFLAGS=-Wall -pedantic -march=native -O2 -ffunction-sections -g -fpic -pipe
# export FFLAGS=-march=native -O2 -ffunction-sections -g -fpic -pipe -Wall -pedantic
# export FCFLAGS=-march=native -O2 -ffunction-sections -g -fpic -pipe -Wall -pedantic


rm ./Full/OncoSimulR/src/OncoSimulR.so
rm ./Full/OncoSimulR/src/OncoSimul.o
rm ./Full/OncoSimulR/src/OncoSimulR_init.o
rm ./Full/OncoSimulR/src/symbols.rds
$V_R CMD INSTALL ./Full/OncoSimulR

