gM3 <- c("Root" = "Root", "d" = "d9, d8",
         "a" = "1, 2", "b" = "3, 4, 5", "c" = "6")
gM2 <- c("Root" = "Root", "a" = "1, 2", "b2" = "3, 4, 5", "b" = "8",
         "c" = "7")

m0 <- data.frame(parent = c("Root", "a", "b"),
                 child  = c("a", "b", "c"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)
epistm1 <- c("a:d" = 0.2, "d:c" = 0.3)
epistm1b <- data.frame(ids = c("a:d", "c:d"), s = c(0.2, 0.3))
oeffects1 <- c("d>a" = 0.4, "c > d" = -0.3)


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
m0 <- data.frame(parent = c("Root", "a", "b"),
                 child  = c("a", "b", "c"),
                 s = 0.1, sh = -1,
                 typeDep = "MN",
                 stringsAsFactors = FALSE)

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

wrap.readFitnessEffects(m0, epistm1,
                        oeffects1, c(0.1, 0.1, 0.2),
                        gM3)
