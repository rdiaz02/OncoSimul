test_that("Exercise LOD and POM code", {
    pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                          "TP53", "TP53", "MLL3"),
                                      child = c("KRAS","SMAD4", "CDNK2A", 
                                          "TP53", "MLL3",
                                          rep("PXDN", 3), rep("TGFBR2", 2)),
                                      s = 0.05,
                                      sh = -0.3,
                                      typeDep = "MN"))
    pancr1 <- oncoSimulIndiv(pancr, model = "Exp", keepPhylog = TRUE)
    pancr8 <- oncoSimulPop(8, pancr, model = "Exp", keepPhylog = TRUE)
    expect_true(inherits(POM(pancr1), "character"))
    require(igraph)
    expect_true(inherits(LOD(pancr1)$all_paths[[1]], "igraph.vs"))
    expect_true(inherits(LOD(pancr1)$lod_single, "igraph.vs"))
    expect_true(inherits(POM(pancr8), "list"))
    expect_true(inherits(LOD(pancr8), "list"))
    expect_true(inherits(diversityPOM(POM(pancr8)), "numeric"))
    expect_true(inherits(diversityLOD(LOD(pancr8)), "numeric"))
    expect_true(diversityPOM(POM(pancr8)) >= 0)
    expect_true(diversityLOD(LOD(pancr8)) >= 0)
    expect_error(diversityPOM(POM(pancr1)),
                 "Object must be a list", fixed = TRUE)
    expect_error(diversityLOD(LOD(pancr1)),
                 "Object must be a list", fixed = TRUE)
})
