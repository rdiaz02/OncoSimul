rm(list = ls())
ff <- dir(pattern = "^test.*R$")
## No plots, as asks for input
ff <- ff[-which(ff == "test.exercise-plotting-code.R")]

system.time({
for(i in ff) {
    ## tt <- system.time(source(i))
    tt <- system.time(test_file(i))
    cat(paste("\n*************************************************  Time of ",
              i, ": ", round(tt[1], 1), "  ***** \n"))
}
})
