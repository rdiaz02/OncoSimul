
create_fe <- function(a, b, c, d, e, f, g,
                        gt = c("WT", "A", "P", "C")) {
  data.frame(Genotype = gt,
             Fitness = c(
                 paste0("1 + ",
                       as.character(d), " * f_1 ",
                       "- ", as.character(c), " * f_3"),
                 paste0("1 - ", as.character(a), " + ", as.character(d), " + ",
                       as.character(f), " * f_1 ",
                      "- ", as.character(c), " * f_3"),
                 paste0("1 + ", as.character(g), " + ",
                       as.character(d), " * f_1 ",
                       "- ", as.character(c), " * (1 + ",
                       as.character(g), ") * f_3"),
                 paste0("1 - ", as.character(b), " + ",
                       as.character(e), " * f_ + ",
                       "(", as.character(d), " + ", as.character(e), ") * f_1 + ",
                       as.character(e) , " * f_2")),
             stringsAsFactors = FALSE)
}


## FIXME: later, add a different set of relationships
## e.g., WT -> P, P -> A, P -> C

create_fe2 <- function(a, b, c, d, e, f, g,
                        gt = c("WT", "A", "P", "C", "A, P", "A, C", "P, C")) {
  data.frame(Genotype = gt,
             Fitness = c(
                 paste0("1 + ",
                       as.character(d), " * f_1_2 ",
                       "- ", as.character(c), " * f_2_3"),
                 "0",
                 paste0("1 - ", as.character(a), " + ", as.character(d), " + ",
                       as.character(f), " * f_1_2 ",
                       "- ", as.character(c), " * f_2_3"),
                 "0",
                 paste0("1 + ", as.character(g), " + ",
                       as.character(d), " * f_1_2 ",
                       "- ", as.character(c), " * (1 + ",
                       as.character(g), ") * f_2_3"),
                 "0",
                 paste0("1 - ", as.character(b), " + ",
                       as.character(e), " * f_ + ",
                       "(", as.character(d), " + ", as.character(e), ") * f_1_2 + ",
                       as.character(e) , " * f_2")),
             stringsAsFactors = FALSE)
}


## Show the expressions, as such
create_fe("a", "b", "c", "d", "e", "f", "g")
create_fe2("a", "b", "c", "d", "e", "f", "g")
## WT -> A, 
