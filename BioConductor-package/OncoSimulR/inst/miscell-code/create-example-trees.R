## Gerstung et al. 2009, poset 2.
p1101 <- rbind(
  c(1, 2),
  c(1, 3),
  c(3, 4),
  c(3, 5),
  c(3, 6),
  c(7, 8),
  c(7, 9),
  c(8, 10),
  c(9, 10),
  c(10, 11)
  )

p1102 <- p1101[-8, ]
p1103 <- p1101[-9, ]
p1104 <- rbind(
  c(1, 2),
  c(1, 3),
  c(3, 6),
  c(3, 8),
  c(7, 9),
  c(7, 4),
  c(9, 10),
  c(2, 11),
  c(10, 5)
  )

p901 <- rbind(
  c(1, 2),
  c(2, 3),
  c(4, 5),
  c(5, 6),
  c(7, 8),
  c(8, 9),
  c(5, 9),
  c(1, 5)
  )

p902 <- rbind(
  c(1, 2),
  c(2, 3),
  c(4, 5),
  c(5, 6),
  c(7, 8),
  c(5, 9)
  )

p903 <- rbind(
  c(1, 2),
  c(2, 3),
  c(5, 6),
  c(7, 8),
  c(8, 9),
  c(1, 5)
  )


## a transformation of p902
p904 <- rbind(
  c(1, 2),
  c(4, 5),
  c(5, 8),
  c(5, 6),
  c(1, 9),
  c(7, 3)
  )


## primary glioblastoma, Fig. 2c in Gerstung et al., 2011.
## p801 <- rbind(
##   c(2, 6),
##   c(1, 2),
##   c(3, 6),
##   c(4, 6),
##   c(5, 6),
##   c(6, 7),
##   c(7, 8))

## p802 <- rbind(
##   c(1, 2),
##   c(2, 6),
##   c(6, 7),
##   c(7, 8)
##   )

## p803 <- rbind(
##   c(1, 8),
##   c(3, 6),
##   c(6, 2),
##   c(4, 7)
##   )


## Pancreatic cancer, figure 2B inGerstung et al., 2011
p701 <- rbind(
  c(1, 2),
  c(1, 3),
  c(1, 4),
  c(1, 5),
  c(2, 6),
  c(3, 6),
  c(4, 6),
  c(4, 7),
  c(5, 7)
  )
## 02, 03, 04: collapsed lead to 701
p702 <- rbind(
  c(1, 2),
  c(1, 3),
  c(1, 4),
  c(1, 5),
  c(2, 6),
  c(4, 7)
  )

p703 <- rbind(
  c(1, 2),
  c(1, 3),
  c(1, 4),
  c(1, 5),
  c(3, 6),
  c(5, 7)
  )
p704 <- rbind(
  c(1, 2),
  c(1, 3),
  c(1, 4),
  c(1, 5),
  c(4, 6),
  c(5, 7)
  )
## moving around
p705 <- rbind(
  c(1, 2),
  c(2, 5),
  c(1, 4),
  c(1, 6),
  c(1, 7),
  c(2, 3)
  )


## Note that all the p?02, 03, etc, try
## to respect the upper level genes, those
## without restrictions.
## Also try not to revert an arrow (biologically nonplausible)
## the 02 and 03 have no conjuctions
##  The 02 and the 03, superposed, lead to the 01 in terms of conjuntcions
##  The 04 is a modification of the 02 with changes in dependencies
##      and the height of events


rt1100 <- poset.to.restrictTable(11)
rt1101 <- poset.to.restrictTable(p1101)
rt1102 <- poset.to.restrictTable(p1102)
rt1103 <- poset.to.restrictTable(p1103)
rt1104 <- poset.to.restrictTable(p1104)


rt900 <- poset.to.restrictTable(9)
rt901 <- poset.to.restrictTable(p901)
rt902 <- poset.to.restrictTable(p902)
rt903 <- poset.to.restrictTable(p903)
rt904 <- poset.to.restrictTable(p904)


## rt800 <- poset.to.restrictTable(8)
## rt801 <- poset.to.restrictTable(p801)
## rt802 <- poset.to.restrictTable(p802)
## rt803 <- poset.to.restrictTable(p803)
## rt804 <- poset.to.restrictTable(p804)


rt700 <- poset.to.restrictTable(7)
rt701 <- poset.to.restrictTable(p701)
rt702 <- poset.to.restrictTable(p702)
rt703 <- poset.to.restrictTable(p703)
rt704 <- poset.to.restrictTable(p704)
rt705 <- poset.to.restrictTable(p705)



## save(file = "./MatherR/data/example_trees.RData", list = ls())
## do not use save.image, because of RNG

## a variety of others. Not used now.

## similar to Gerstung et al., Poset 1, but with 13 and 14
## p1 <- rbind(
##   c(1, 5),
##   c(2, 5),
##   c(1, 6),
##   c(2, 6),
##   c(2, 7),
##   c(1, 4),
##   c(4, 9),
##   c(4, 10),
##   c(6, 8),
##   c(7, 8),
##   c(10, 11),
##   c(8, 11),
##   c(8, 12),
##   c(9, 12),
##   c(12, 13),
##   c(13, 14))

## p2 <- rbind(
##   c(1, 5),
##   c(2, 5),
##   c(1, 6),
##   c(2, 6),
##   c(2, 7),
##   c(1, 4),
##   c(4, 9),
##   c(4, 10),
##   c(6, 8),
##   c(7, 8),
##   c(10, 11),
##   c(8, 11),
##   c(8, 12),
##   c(9, 12),
##   c(12, 13),
##   c(13, 14),
##   c(3, 15),
##   c(3, 16),
##   c(7, 17),
##   c(7, 18),
##   c(18, 19),
##   c(18, 20))
## ##plot(poset.to.graph(p2, c("null", 1:max(p2))))

## ## no conjuntions
## p3 <- rbind(
##   c(1, 5),
##   c(1, 6),
##   c(2, 7),
##   c(1, 4),
##   c(4, 9),
##   c(4, 10),
##   c(6, 8),
##   c(10, 11),
##   c(8, 12),
##   c(12, 13),
##   c(13, 14),
##   c(3, 15),
##   c(3, 16),
##   c(7, 17),
##   c(7, 18),
##   c(18, 19),
##   c(18, 20))
## ##plot(poset.to.graph(p3, c("null", 1:max(p3))))

## # only two conjunction
## p4 <- rbind(
##   c(1, 5),
##   c(2, 5),
##   c(2, 6),
##   c(2, 7),
##   c(1, 4),
##   c(4, 9),
##   c(4, 10),
##   c(6, 8),
##   c(7, 8),
##   c(10, 11),
##   c(5, 12),
##   c(12, 13),
##   c(2, 14))

## ##plot(poset.to.graph(p4, c("null", 1:max(p4))))

## p5 <- rbind(
##   c(1, 5),
##   c(2, 5),
##   c(1, 6),
##   c(3, 6),
##   c(2, 7),
##   c(3, 7),
##   c(2, 8),
##   c(4, 8),
##   c(3, 9),
##   c(4, 9),
##   c(1, 10),
##   c(1, 11),
##   c(4, 12),
##   c(4, 13),
##   c(10, 14),
##   c(11, 15),
##   c(12, 16),
##   c(13, 17))

## ##plot(poset.to.graph(p5, c("null", 1:max(p5))))


## p6 <- rbind(
##   c(1, 5),
##   c(2, 5),
##   c(1, 6),
##   c(3, 6),
##   c(2, 7),
##   c(3, 7),
##   c(2, 8),
##   c(4, 8),
##   c(3, 9),
##   c(4, 9),
##   c(1, 10),
##   c(1, 11),
##   c(4, 12),
##   c(4, 13),
##   c(10, 14),
## #  c(11, 15),
##   c(12, 16),
##   c(13, 17),
## #  c(5, 18),
##   c(6, 19),
##   c(7, 20),
##   c(8, 15),
##   c(9, 18)
##   )

## ##plot(poset.to.graph(p6, c("null", 1:max(p6))))

## ## Gerstung et al., Poset 2.
## p7 <- rbind(
##   c(1, 2),
##   c(1, 3),
##   c(3, 4),
##   c(3, 5),
##   c(3, 6),
##   c(7, 8),
##   c(7, 9),
##   c(8, 10),
##   c(9, 10),
##   c(10, 11)
##   )
## ##plot(poset.to.graph(p7, c("null", 1:max(p7))))


## p8 <- rbind(
##   c(1, 2),
##   c(1, 3),
  
##   c(3, 4),
##   c(3, 5),
  
##   c(6, 7),
##   c(6, 8),
  
##   c(9, 10),
##   c(9, 11),
  
##   c(12, 13),
##   c(12, 14),
##   c(12, 15)
##   )

## p8 <- rbind(
##   c(1, 2),
##   c(1, 3),
  
##   c(3, 4),
##   c(3, 5),
  
##   c(6, 7),
##   c(6, 8),
  
##   c(9, 10),
##   c(9, 11),
  
##   c(12, 13),
##   c(12, 14),
##   c(12, 15)
##   )

## ## similar to 8, but with more conjunctions
## p9 <- rbind(
##   c(1, 2),
##   c(2, 3),
##   c(2, 4),

##   c(5, 6),
##   c(1, 14),
##   c(5, 14),

##   c(6, 7),
##   c(8, 9),
##   c(8, 10),

##   c(9, 11),
##   c(10, 11),
##   c(6, 15),
##   c(9, 15),
##   c(12, 13)
##   )

## ## variations on a theme with 9 drivers
## p10 <- rbind(
##   c(1, 2),
##   c(2, 3),
##   c(3, 4),
##   c(4, 5),
##   c(6, 7),
##   c(7, 8),
##   c(8, 9)
##   )

## p11 <- rbind(
##   c(1, 2),
##   c(2, 3),
##   c(3, 4),
##   c(4, 5),
##   c(5, 6),
##   c(6, 7),
##   c(8, 9)
##   )


## p12 <- rbind(
##   c(1, 2),
##   c(2, 3),
##   c(3, 4),
##   c(4, 5),
##   c(6, 7),
##   c(7, 8),
##   c(8, 9),
##   c(7, 3),
##   c(2, 8)
##   )

## p13 <- rbind(
##   c(1, 2),
##   c(2, 3),
##   c(3, 4),
##   c(4, 5),
##   c(6, 7),
##   c(7, 8),
##   c(8, 9),
##   c(7, 3),
##   c(2, 8),
##   c(3, 9),
##   c(8, 4)
##   )

## p14 <- rbind(
##   c(1, 2),
##   c(2, 3),
##   c(4, 5),
##   c(5, 6),
##   c(7, 8),
##   c(8, 9)
##   )

## p15 <- rbind(
##   c(1, 2),
##   c(2, 3),
##   c(4, 5),
##   c(5, 6),
##   c(7, 8),
##   c(8, 9),
##   c(1, 5),
##   c(2, 9)
##   )

## rt1 <- poset.to.restrictTable(p1)
## rt2 <- poset.to.restrictTable(p2)
## rt3 <- poset.to.restrictTable(p3)
## rt4 <- poset.to.restrictTable(p4)
## rt5 <- poset.to.restrictTable(p5)
## rt6 <- poset.to.restrictTable(p6)
## rt7 <- poset.to.restrictTable(p7)
## rt8 <- poset.to.restrictTable(p8)
## rt9 <- poset.to.restrictTable(p9)
## rt10 <- poset.to.restrictTable(p10)
## rt11 <- poset.to.restrictTable(p11)
## rt12 <- poset.to.restrictTable(p12)
## rt13 <- poset.to.restrictTable(p13)
## rt14 <- poset.to.restrictTable(p14)
## rt15 <- poset.to.restrictTable(p15)

## p1106 <- p1101[-8, ]
## p1106[9, 1] <- 8
## p1106[3, 1] <- 2
## p1105 <- p1102
## p1105[1, ] <- c(2, 1)
## p1105[2, 1] <- 2
## p1105[9, ] <- c(11, 10)
## p1105[8, ] <- c(9, 11)



## p905 <- p902
## p905[1, 2] <- 3
## p905[2, ] <- c(3, 2)
## p905[5, ] <- c(8, 7)
## p905[6, 1] <- 7
## p905[3, ] <- c(6, 4)
## p905[4, ] <- c(4, 5)
## p9b01 <- rbind(
##   c(1, 2),
##   c(2, 3),
##   c(4, 5),
##   c(5, 6),
##   c(7, 8),
##   c(8, 9),
##   c(4, 2),
##   c(5, 9)
##   )
## p9b02 <- p9b01
## p9b02 <- p9b01[-c(7, 8),]
## p9b04 <- p9b02
## p9b04[1, 2] <- 3
## p9b04[2, ] <- c(3, 2)
## p9b04[5, ] <- c(8, 7)
## p9b04[6, 1] <- 7
## p9b04[3, ] <- c(6, 4)
## p9b04[4, ] <- c(4, 5)


## p702 <- rbind(
##   c(1, 2),
##   c(1, 3),
##   c(2, 5),
##   c(2, 4),
##   c(3, 6),
##   c(3, 7)
##   )

## p703 <- rbind(
##   c(1, 2),
##   c(2, 5),
##   c(1, 4),
##   c(1, 6),
##   c(1, 7),
##   c(2, 3)
##   )

## p808 <- rbind(
##   c(1, 2),
##   c(3, 6),
##   c(6, 8),
##   c(5, 7)
##   )
## p805 <- rbind(
##   c(2, 6),
##   c(3, 6),
##   c(4, 6),
##   c(6, 7),
##   c(7, 8))
## p806 <- rbind(
##   c(2, 6),
##   c(3, 6),
##   c(5, 6),
##   c(6, 7),
##   c(7, 8))
## p807 <- rbind(
##   c(2, 6),
##   c(3, 6),
##   c(5, 6),
##   c(6, 8),
##   c(8, 7))



## fp <- function(x, ...) plot(poset.to.graph(x, c("null", 1:max(x))), ...)
## fp(p1)
## fp(p2)
## fp(p3)
## fp(p4)
## fp(p5)
## fp(p6)
## fp(p7)
## fp(p8)
## fp(p9)

## Old stuff used with these names for parameter searches
## rt.p11 <- poset.to.restrictTable(poset.11)
## rt.p11C <- poset.to.restrictTable(poset.11C)
## p7A <- p7
## p7B <- p7[-8,]

## p14A <- p14
## p14B <- rbind(
##   c(1, 2),
##   c(2, 3),
##   c(4, 5),
##   c(5, 6),
##   c(7, 8),
##   c(8, 9),
##   c(5, 9),
##   c(4, 2)
##   )


## ## name using number of nodes and conjunct or not

## poset.9 <- p14
## poset.9C <- p14B
## poset.11 <- p7B
## poset.11C <- p7A ## the one like in Gerstung et al.


## rt.p9 <- poset.to.restrictTable(poset.9)
## rt.p9C <- poset.to.restrictTable(poset.9C)

## rm(p14A, p14B, p7A, p7B)
