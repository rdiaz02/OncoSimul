## Copyright 2016 Ramon Diaz-Uriarte

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.



## Functions to obtain LOD and POM similar to Szendro et al., 2014, PNAS.
genot_max <- function(x) {
    x$GenotypesLabels[which.max(x$pops.by.time[nrow(x$pops.by.time), -1])]
}

LOD.internal <- function(x) {
    ## Not identical to LOD of Szendro because:
    
    ##  a) I think we want all paths, not just their single LOD, which I
    ##  think they might use out of convenience.

    ##  b) For me it is a mess, a complicated mess, to use their LOD as
    ##  they define it and there are many ambiguities in how to define it
    ##  in continuous time.

    ## This also means that single simulation might yield multiple LODs

    ## keepEvents is FALSE to make this object as small as possible.
    if(is.null(x$pops.by.time)) {
        warning("Missing needed components. This might be a failed simulation.",
                " Returning NA.")
        return(list(all_paths = NA, lod_single = NA))
    }
    pc <- phylogClone(x, keepEvents = FALSE)
    
    if((length(pc) == 1) && (is.na(pc))) {
        return(list(all_paths = NA,
                    lod_single = "No_descendants"))
    }
    pcg <- pc$graph
    end <- genot_max(x)
    all_paths <- igraph::all_simple_paths(pcg, from = "", to = end, mode = "out")
    ## the next is partially redundant
    ## graph_to_end <- igraph::make_ego_graph(pcg, order = 1e9, nodes = end,
    ##                                        mode = "in")
    ## if(length(graph_to_end) != 1) stop("length(graph_to_end) > 1")
    ## I am not sure if I should keep the last one. Redundant

    ## This gives a single path and it is the first entry into each
    ## destination. But we do not check there is no extinction afterwards.
    ## The closest to the single Szendro LOD
    if(end == "") {
        ## Max is WT
        lod_single <- "WT_is_end"
    } else {
        singlep <- pc$df
        singlep[, 1] <- as.character(singlep[, 1])
        singlep[, 2] <- as.character(singlep[, 2])
        singlep <- singlep[ do.call(order, singlep[, c(2, 3)]), ]
        singlep <- singlep[!duplicated(singlep[, 2]), ]
        gsingle <- igraph::graph_from_data_frame(singlep)
        lod_single <- igraph::all_simple_paths(gsingle, from = "", to = end, mode = "out")
        if(length(lod_single) != 1) stop("lod_single != 1")
    }
    return(list(all_paths = all_paths,
                lod_single = lod_single[[1]])) ##, graph_to_end = graph_to_end[[1]]))
                ## graph_phylog_clone = pcg))
}


POM.internal <- function(x) {
    if(is.null(x$pops.by.time)) {
        warning("Missing needed components. This might be a failed simulation.",
                " Returning NA.")
        return(NA)
    }
    x$GenotypesLabels[rle(apply(x$pops.by.time[, -1, drop = FALSE], 1, which.max))$values]
}

## First do, over a set of simulations, sim:
## l_lod_single <- mclapply(sim, LODs)
## l_poms <- mclapply(sim, POM)


diversityLOD <- function(llod) {
    nn <- names(llod[[1]])
    if( is_null_na(nn) ||
        !(nn == c("all_paths", "lod_single")))
        stop("Object must be a list of LODs")
    pathstr <- unlist(lapply(llod, function(x) paste(names(x$lod_single),
                                                     collapse = "_")))
    shannonI(table(pathstr))
}

diversityPOM <- function(lpom) {
    if(!inherits(lpom, "list"))
        stop("Object must be a list of POMs")
    ## if(!inherits(x, "oncosimul_lod_list"))
    ##     stop("This is not a list of POMs")
    pomstr <- unlist(lapply(lpom, function(x) paste(x, collapse = "_")))
    shannonI(table(pomstr))
}

## a legacy
diversity_POM <- diversityPOM
diversity_LOD <- diversityLOD

## POM.oncosimul2 <- POM.internal
## LOD2.oncosimul2 <- LOD.internal

POM <- function(x) {
    UseMethod("POM", x)
}

LOD <- function(x) {
    UseMethod("LOD", x)
}

POM.oncosimul2 <- function(x) return(POM.internal(x))
LOD.oncosimul2 <- function(x) return(LOD.internal(x))

POM.oncosimulpop <- function(x) return(lapply(x, POM.internal))
LOD.oncosimulpop <- function(x) return(lapply(x, LOD.internal))

## POM.oncosimul2 <- function(x) {
##     out <- POM.internal(x)
##     class(out) <- c(class(out), "oncosimul_pom")
##     return(out)
## }

## LOD.oncosimul2 <- function(x) {
##     out <- LOD.internal(x)
##     class(out) <- c(class(out), "oncosimul_lod")
##     return(out)
## }


## POM.oncosimulpop <- function(x) {
##     out <- lapply(x, POM.internal)
##     class(out) <- c(class(out), "oncosimul_pom_list")
##     return(out)
## }

## LOD.oncosimulpop <- function(x) {
##     out <- lapply(x, LOD.internal)
##     class(out) <- c(class(out), "oncosimul_lod_list")
##     return(out)
## }


## summary.oncosimul_lod_list <- function(x) {
##     cat("List of ", length(x), " simulations\n.")
##     cat("Shannon's diversity (entropy) = ", diversityLOD(x), "\n")
## }

## summary.oncosimul_pom_list <- function(x) {
##     cat("List of ", length(x), " simulations\n.")
##     cat("Shannon's diversity (entropy) = ", diversityPOM(x), "\n")
## }


