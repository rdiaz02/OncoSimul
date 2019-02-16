
pancr <- allFitnessEffects(data.frame(parent = c("Root", rep("KRAS", 4), "SMAD4", "CDNK2A", 
                                               "TP53", "TP53", "MLL3"),
                                           child = c("KRAS","SMAD4", "CDNK2A", 
                                               "TP53", "MLL3",
                                               rep("PXDN", 3), rep("TGFBR2", 2)),
                                           s = 0.05,
                                           sh = -0.3,
                                           typeDep = "MN"))
     
     
pancr1 <- oncoSimulIndiv(pancr, model = "Exp", keepPhylog = TRUE)

## All we need to get LOD sensu stricto (i.e., identical to Szendro)
## is keep pop size of receiving, or destination, genotype.
## Then, filter those where popSize > 0


## And use first (starting from bottom) of that path
## So find, for each child, the last event with popSize == 0,
## and keep only that row.

## Probably enough to run duplicated in reverse (on the df with popSize
## child = 0)

## The indices to keep: !rev(duplicated(rev(fg3[, 2])))

fg1 <- pancr1$other$PhylogDF


fg3 <- data.frame(parent = c("", "", "A", "B", "C", "A, B"),
                  child =  c("A","B","A, B","A, B", "A, C", "A, B, C"),
                  time = 1:6,
                  pop_size_child = c(0, 0, 0, 0, 0, 0),
                  stringsAsFactors = FALSE)


fg4 <- data.frame(parent = c("", "", "A", "B", "C", "A, B"),
                  child =  c("A","B","A, B","A, B", "A, C", "A, B, C"),
                  time = 1:6,
                  pop_size_child = c(0, 0, 0, 2, 0, 0),
                  stringsAsFactors = FALSE)

## from phylogClone, key parts
fpc <- function(df) {
    tG <- unique(c(df[, 1], df[, 2]))
    g <- igraph::graph.data.frame(df[, c(1, 2)])
    nodesInP <- unique(unlist(igraph::neighborhood(g, order = 1e+09, 
                                                   nodes = tG, mode = "in")))
    allLabels <- unique(as.character(unlist(df[, c(1, 2)])))
    nodesRm <- setdiff(allLabels, V(g)$name[nodesInP])
    g <- igraph::delete.vertices(g, nodesRm)
    tmp <- list(graph = g, df = df)
    class(tmp) <- c(class(tmp), "phylogClone")
    return(tmp)
}

## Filter the PhylogDF so we obtain LOD, sensu stricto.
filter_phylog_df_LOD <- function(x) {
    x <- x[x$pop_size_child == 0, ]
    keep <- !rev(duplicated(rev(x$child)))
    return(x[keep, ])
}

all_simple_paths(fpc(filter_phylog_df_LOD(fg3))$graph, from = "",
                 to = "A, B, C",
                 mode = "out")


all_simple_paths(fpc(clean_phylog_df(fg3))$graph, to = "",
                 from = "A, B, C",
                 mode = "in")

## same, for lod
lod1 <- function(x) {
    
}





### Do it as follows:
##  Always obtain POM and LOD, the sensu stricto version
##  When running with keepPhylog, allow returning also
##  the many LODs

## LOD_all: the function that requires the "keepPhylog"

## We have two structures in C++, and the LOD one is like the phylog, w/o
## the popSize

## the LOD always kept only called in the place of the code where we
## create new species

## POM is obtained directly in the C++ code.

## Of course, change help and vignette
