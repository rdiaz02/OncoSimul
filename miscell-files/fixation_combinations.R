## Already added

## Recall all the genes in mutator must be also part of fitness effects,
## so this is enough

## namedGenes already exists
## namedGenes <- OncoSimulR:::allNamedGenes(fnme)

fixed_comb <- list(1:2, 3)

fixed_comb <- list("a1", c("a2", "a1"), c("a1", "b1"))


ng <- namedGenes
rownames(ng) <- namedGenes[, "Gene"]
fixed_comb_int <- lapply(fixed_comb, function(x) sort(ng[x, 2]))


## checks for missing
## checks for not present
