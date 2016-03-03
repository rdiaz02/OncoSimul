## Run the coverage script

library(covr)
library(shiny)
library(DT)
library(OncoSimulR)
setwd("~/Proyectos/OncoSimul/OncoSimulR")

cov3 <- package_coverage(type = "all")
cov4 <- package_coverage()
zero_coverage(cov4) 
save(file = "coverage-results.RData", cov3, cov4)
shine(cov3)
