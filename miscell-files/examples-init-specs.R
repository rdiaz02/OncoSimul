source("new-restrict.R")

sd <- 0.1
sdp <- 0.15
sp <- 0.05
bauer <- data.frame(parent = c("Root", rep("p", 5)),
                    child = c("p", paste0("s", 1:5)),
                    s = c(sd, rep(sdp, 5)),
                    sh = c(0, rep(sp, 5)),
                    typeDep = "MN",
                    stringsAsFactors = FALSE)

b1 <- evalAllGenotypes(allFitnessEffects(bauer), order = FALSE)
b2 <- evalAllGenotypes(allFitnessEffects(bauer), order = TRUE, max = 2000)




m0 <- data.frame(parent = c("Root", "a", "b"),
                 child  = c("a", "b", "c"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)

epistm1 <- c("a:d" = 0.2, "d:c" = 0.3)
epistm1b <- data.frame(ids = c("a:d", "c:d"), s = c(0.2, 0.3))
oeffects1 <- c("d>a" = 0.4, "c > d" = -0.3)


epineg <- c("-a:d" = 0.2, "a:d" = 0.3, "d:c" = 0.3)
epineg2 <- c("-a:d" = 0.2, "b:c" = 0.3)
epineg3 <- c("a:-d" = 0.2, "b:c" = 0.3)


allFitnessEffects(epistasis = epineg, geneToModule = NULL)

allFitnessEffects(epistasis = epineg2, geneToModule = NULL)

allFitnessEffects(epistasis = epineg3, geneToModule = NULL)



gme <- c("Root" = "Root", "a" = "1, 2", "d" = 3, "c" = 4)

gme2 <- c("Root" = "Root", "a" = "1, 2", "d" = 3, "b" = "5, 6", "c" = 4)

allFitnessEffects(epistasis = epineg, geneToModule = gme)
allFitnessEffects(epistasis = epineg2, geneToModule = gme2)

gM3 <- c("Root" = "Root", "d" = "d9, d8",
         "a" = "1, 2", "b" = "3, 4, 5", "c" = "6")



allFitnessEffects(m0, epistasis = epineg, geneToModule = gM3)
allFitnessEffects(m0, epistasis = epineg2, geneToModule = gM3)


oo <- allFitnessEffects(m0)
oo2 <- allFitnessEffects(m0, epistm1)

evalAllGenotypes(oo2)


oa <- allFitnessEffects(m0, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2), gM3)

## this should not crash the code
ovalRGenotype(0, oa, verbose = TRUE)

ovalRGenotype(c(1, 2), oa, verbose = TRUE)

evalGenotype(c("d8", "2", "6"), oa, verbose = TRUE)


## errors are caught in genotype specification

evalGenotype(c("d8898", "2", "6"), oa, verbose = TRUE)







oa2 <- allFitnessEffects(m0, epistm1,
                        oeffects1, runif(1000), gM3)


benchmark(wrap.readFitnessEffects(m0, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2),
                                  gM3, echo = FALSE),
          replications = 100)

benchmark(allFitnessEffects(m0, epistm1,
                            oeffects1, c(0.1, 0.1, 0.2), gM3),
          replications = 100)

benchmark(readFitnessEffects(oa2, echo = FALSE),
          replications = 10000)
benchmark(readFitnessEffects(oa, echo = FALSE),
          replications = 10000)



microbenchmark(readFitnessEffects(oa, echo = FALSE), times = 1000)

microbenchmark(allFitnessEffects(m0, epistm1,
                            oeffects1, c(0.1, 0.1, 0.2), gM3),
          times = 100)

evalGenotype(c("d8", "2", "6"), oa, verbose = TRUE)





### examples here





gM2 <- c("Root" = "Root", "a" = "1, 2", "b2" = "3, 4, 5", "b" = "8",
         "c" = "7")

to.long.rt(m0, gm.to.geneModuleL(gM))

to.long.rt(m0, gm.to.geneModuleL(gM2))

m0 <- data.frame(parent = c("Root", "a", "b"),
                 child  = c("a", "b", "c"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)


gM <- c("Root" = "Root", "a" = "1, 2", "b" = "3, 4, 5", "c" = "6")
gm.to.geneModuleL(gM)


allFitnessEffects(m0)

allFitnessEffects(m0, noIntGenes = c(0.1, 0, 0.2))

allFitnessEffects(m0, noIntGenes = c("u" = 0.1, "v" = 0, "mm" = 0.2))


rtAndGeneModule(m0, gM)
rtAndGeneModule(m0)


epistm1 <- c("a:d" = 0.2, "d:c" = 0.3)
epistm1b <- data.frame(ids = c("a:d", "c:d"), s = c(0.2, 0.3))
oeffects1 <- c("d>a" = 0.4, "c > d" = -0.3)

allFitnessEffects(m0, epistasis = epistm1,
                  orderEffects = oeffects1,
                  noIntGenes = c(0.1, 0, 0.2))

wrap.readFitnessEffects(NULL,
                        NULL,
                        NULL,
                        c(0.1, 0.1),
                        NULL)

wrap.readFitnessEffects(m0,
                        NULL,
                        NULL,
                        c(0.1, 0.1),
                        NULL)

wrap.readFitnessEffects(m0,
                        NULL,
                        NULL,
                        NULL,
                        NULL)

wrap.readFitnessEffects(NULL,
                        NULL,
                        NULL,
                        NULL,
                        NULL)


wrap.readFitnessEffects(m0, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2),
                        NULL)

wrap.readFitnessEffects(m0, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2),
                        gM2)


m00 <- data.frame(parent = c("Root", "c"),
                 child  = c("c", "b"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)


m000 <- data.frame(parent = c("Root", "d", "c", "b"),
                 child  = c("c", "a", "a", "d"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)

## good to check geneModule table is OK
gM4 <- c("Root" = "Root", "d" = "d9, d8, a, z9",
         "a" = "z700, u43, 78", "b" = "2, 3, 4, 5", "c" = "1, b, 6")

wrap.readFitnessEffects(m000, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2),
                        gM4)

## a complex one with originally disordered parents and children and weird
## dep on 0 and others.
m4 <- data.frame(parent = c(rep("Root", 7), "d", "c", "b", "g", "h", "d", "c", "c", "e"),
                 child  = c(letters[2:8], rep("a", 5), rep("b", 2), rep("g", 2)),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)

wrap.readFitnessEffects(m4, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2), NULL)


wrap.readFitnessEffects(m4, epineg,
                        oeffects1, c(0.1, 0.1, 0.2), NULL)



benchmark(wrap.readFitnessEffects(m0, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2),
                                  gM3, echo = FALSE),
          replications = 100)
## FIXME make sure to test with 0 size elements: rT, epist, order


