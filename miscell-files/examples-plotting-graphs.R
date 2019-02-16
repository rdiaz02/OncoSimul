## not a particularly successful thing

ig1 <- fea$graph
g1 <- igraph.to.graphNEL(ig1)


c1 <- unlist(lapply(edgeData(g1), function(x) x$color))
names(c1) <- sub("|", "~", names(c1), fixed = TRUE)
## plot(g1, edge = list(color = u1)) ## nope
eAttrs <- list()
eAttrs$color <- c1

s1 <- unlist(lapply(edgeData(g1), function(x) x$lty))
## s2 <- rep("solid", length(s1))
## s2[s1 == 2] <- "dashed"
## s2[s1 == 3] <- "dotted"
names(s1) <- names(c1)



a1 <- unlist(lapply(edgeData(g1), function(x) max(x$arrow.mode - 1, 0)))
names(a1) <- names(c1)
## nope, arrowhead does not work
plot(g1, edgeAttrs = list(arrowsize = a1, style = s2, color = c1))

lwd <- s1
lwd[lwd == 2] <- 2
lwd[lwd == 3] <- 3

plot(g1, edgeAttrs = list(arrowsize = a1, lty = s1, lwd = lwd, color = c1))




edgeRenderInfo(g1) <- list(lty = "dashed")


plot(g1, edgeAttrs = list(style = s2))




plot(pag, layout = layout.circle(pag))
plot(pag, layout = layout.auto(pag))

plot(pag, layout = layout.reingold.tilford) ## tree


plot(pag, layout = layout.graphopt(pag))
plot(pag, layout = layout.norm(layout.auto(pag)))
