#!/bin/bash

V_R=$1
export R_MAKEVARS_USER=/home/ramon/.R/Makevars


rm ./OncoSimulR/src/OncoSimulR.so
rm ./OncoSimulR/src/OncoSimul.o
rm ./OncoSimulR/src/OncoSimulR_init.o
rm ./OncoSimulR/src/symbols.rds

$V_R CMD INSTALL ./OncoSimulR

