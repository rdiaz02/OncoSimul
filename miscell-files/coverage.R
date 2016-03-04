## Run the coverage script

library(covr)
library(shiny)
library(DT)
library(OncoSimulR)
setwd("~/Proyectos/OncoSimul/OncoSimulR")

cov3 <- package_coverage(type = "all")
cov4 <- package_coverage()
save(file = "../miscell-files/coverage-results.RData", cov3, cov4)
zero_coverage(cov4) 
shine(cov3)
