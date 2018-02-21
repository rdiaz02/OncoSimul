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
