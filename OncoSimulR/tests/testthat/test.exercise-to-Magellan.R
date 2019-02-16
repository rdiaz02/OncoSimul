inittime <- Sys.time()
cat(paste("\n Starting to_Magella at", date(), "\n"))
test_that("Expect output", {
    r1 <- rfitness(4)
    expect_output(to_Magellan(r1, NULL))
    cs <-  data.frame(parent = c(rep("Root", 3), "a", "d", "c"),
                      child = c("a", "b", "d", "e", "c", "f"),
                      s = 0.1,
                      sh = -0.9,
                      typeDep = "MN")
    expect_output(to_Magellan(allFitnessEffects(cs), NULL))
})
cat(paste("\n Ending to_Magella at", date(), "\n"))
cat(paste("  Took ", round(difftime(Sys.time(), inittime, units = "secs"), 2), "\n\n"))
rm(inittime)
