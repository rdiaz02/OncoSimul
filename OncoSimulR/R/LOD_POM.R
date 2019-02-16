## Copyright 2016, 2017 Ramon Diaz-Uriarte

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


## Note that, in contrast to POM, LOD is not well defined if the population
## becomes extinct.

## Functions to obtain LOD and POM similar to Szendro et al., 2014, PNAS.
genot_max <- function(x) {
    x$GenotypesLabels[which.max(x$pops.by.time[nrow(x$pops.by.time), -1])]
}



## Filter the PhylogDF so we obtain LOD, sensu stricto.

## No longer necessary with LOD from C++ removing duplicates
## ## Now this is coming from LOD_DF, which already only has
## ## implicitly pop_size_child == 0
## Only used in one function use for testing
filter_phylog_df_LOD <- function(y) {
    keep <- !rev(duplicated(rev(y$child)))
    return(y[keep, ])
}



## ## Filter the PhylogDF so we obtain LOD, sensu stricto.
## ## For the old version
## filter_phylog_df_LOD_with_n <- function(y) {
##     y <- y[y$pop_size_child == 0, , drop = FALSE]
##     keep <- !rev(duplicated(rev(y$child)))
##     return(y[keep, ])
## }


## from phylogClone, key parts for the LOD strict structure
phcl_from_lod <- function(df, x) {
    ## no longer necessary
    ## ## the !keepEvents. Which I move here to speed things up.
    ## df <- df[!duplicated(df[, c(1, 2)]), , drop = FALSE]
    
    tG <- unique(c(as.character(df[, 1]), as.character(df[, 2])))
    ## ## Do as in phylogClone. So that we have the same nodes
    ## ## in LOD all and not LOD all?
    ## z <- which_N_at_T(x, N = 1, t = "last")
    ## tG <- x$GenotypesLabels[z]
    
    ## ## FIXME: aren't these two warnings redundant or aliased?
    ## ## yes, partlt
    ##  I think this can never happen now
    ## if ((length(tG) == 1) && (tG == "")) {
    ##     warning("There never was a descendant of WT")
    ## }
    if (nrow(df) == 0) {
        warning("LOD structure has 0 rows: no descendants of initMutant ever appeared. ")
        return(NA)
    }
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

LOD.internal <- function(x) {
    if(is.null(x$pops.by.time)) {
        warning("Missing needed components. This might be a failed simulation.",
                " Returning NA.")
        return(list(all_paths = NA, lod_single = NA))
    }
    if (!inherits(x, "oncosimul2")) 
        stop("LOD information is only stored with v >= 2")
    ## y <- filter_phylog_df_LOD(x$other$LOD_DF)
    y <- x$other$LOD_DF
    pc <- phcl_from_lod(y)
    ## need eval for oncoSimulPop calls and for LOD_as_path
    initMutant <- x$InitMutant

    ## No descendants means that: never a descendant.
    ## Note the same that the init mutant be the final state.
    if((length(pc) == 1) && (is.na(pc))) {
        lod <- "No_descendants"
        ## bail out here. We do not need the rest.
        if(initMutant != "")
            attributes(lod)$initMutant <- initMutant
        return(lod)
    }
    
    pcg <- pc$graph
    end <- genot_max(x)
    
    if(end == initMutant) {
        if(initMutant == "") {
            stinitm <- "WT"
        } else {
            stinitm <- paste0("initMutant(", initMutant, ")")
        }
        lod <- paste0(stinitm, "_is_end")
    } else {
        all_paths <- igraph::all_simple_paths(pcg, from = initMutant, to = end,
                                              mode = "out")
        if(length(all_paths) > 1)
            stop("length(all_paths) > 1???")
        lod <- igraph::as_ids(all_paths[[1]])
    }
    if(initMutant != "")
        attributes(lod)$initMutant <- initMutant
    return(lod)
}


POM.internal <- function(x) {
    if(is.null(x$pops.by.time)) {
        warning("Missing needed components. This might be a failed simulation.",
                " Returning NA.")
        return(NA)
    } else {
        x$other$POM
    }
}



diversityLOD <- function(llod) {
    ## nn <- names(llod[[1]])
    nn <- llod[[1]]
    if( is_null_na(nn) ||
        !(is.list(llod)))
        stop("Object must be a list of LODs")
    pathstr <- unlist(lapply(llod, function(x) paste(x,
                                                     collapse = "_")))
    shannonI(table(pathstr))
}

## diversityLOD <- function(llod) {
##     nn <- names(llod[[1]])
##     if( is_null_na(nn) ||
##         !(is.list(llod)))
##         stop("Object must be a list of LODs")
##     pathstr <- unlist(lapply(llod, function(x) paste(names(x),
##                                                      collapse = "_")))
##     shannonI(table(pathstr))
## }

LOD_as_path <- function(llod) {
    path_l <- function(u) {
        if(length(u) == 1) {
            if(is.null(attributes(u)$initMutant))
                initMutant <- ""
            else
                initMutant <- attributes(u)$initMutant
            if(initMutant == "") initMutant <- "WT"
            if(grepl("_is_end", u))
                return(initMutant)
            if(u == "No_descendants")
                return(initMutant)
        } else {
            ## Deal with "" meaning WT
            ## the_names <- names(u)
            the_names <- u
            the_names_wt <- which(the_names == "")
            
            if(length(the_names_wt)) {
                if(length(the_names_wt) > 1) stop("more than 1 WT?!")
                if(the_names_wt > 1) stop("WT in position not 1?!")
                the_names[the_names_wt] <- "WT"
            }
            return(paste(the_names, collapse = " -> ")) 
        }
    }
    if(!is.list(llod))
        pathstr <- path_l(llod)
    else
        pathstr <- unlist(lapply(llod, path_l))
    return(pathstr)
}
## We would just need a LOD_as_DAG

diversityPOM <- function(lpom) {
    if(!inherits(lpom, "list"))
        stop("Object must be a list of POMs")
    pomstr <- unlist(lapply(lpom, function(x) paste(x, collapse = "_")))
    shannonI(table(pomstr))
}

## a legacy
diversity_POM <- diversityPOM
diversity_LOD <- diversityLOD


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



POM_pre_2.9.2 <- function(x) {
    if(is.null(x$pops.by.time)) {
        warning("Missing needed components. This might be a failed simulation.",
                " Returning NA.")
        return(NA)
    }
    x$GenotypesLabels[rle(apply(x$pops.by.time[, -1, drop = FALSE],
                                1, which.max))$values]
}



LOD.internal_pre_2.9.2 <- function(x, strict) {
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
    if(strict) {
        if (!inherits(x, "oncosimul2")) 
            stop("LOD information is only stored with v >= 2")
        y <- filter_phylog_df_LOD(x$other$LOD_DF)
        pc <- phcl_from_lod(y)
    } else {
        pc <- phylogClone(x, keepEvents = FALSE)
    }
    ## need eval for oncoSimulPop calls and for LOD_as_path
    initMutant <- x$InitMutant
    
    if((length(pc) == 1) && (is.na(pc))) {
        lodlist <- list(all_paths = NA,
                        lod_single = "No_descendants")
        ## bail out here. We do not need the rest.
        if(initMutant != "")
            attributes(lodlist)$initMutant <- initMutant
        return(lodlist)
    }
    
    pcg <- pc$graph
    end <- genot_max(x)
    
    ## if(!is.null(eval(attributes(x)$call$initMutant))) {
    ##     initMutant <- eval(attributes(s7)$call$initMutant)
    ## } else {
    ##     initMutant <- ""
    ## }
    ## browser()
    if(end == initMutant) {
        if(initMutant == "") {
            stinitm <- "WT"
        } else {
            stinitm <- paste0("initMutant(", initMutant, ")")
        }
        lod_single <- paste0(stinitm, "_is_end")
        all_paths <- list(lod_single)
    } else {
        all_paths <- igraph::all_simple_paths(pcg, from = initMutant, to = end,
                                              mode = "out")
       
        if(!strict) {
            ## the next is partially redundant
            ## graph_to_end <- igraph::make_ego_graph(pcg, order = 1e9, nodes = end,
            ##                                        mode = "in")
            ## if(length(graph_to_end) != 1) stop("length(graph_to_end) > 1")
            ## I am not sure if I should keep the last one. Redundant
            
            ## This gives a single path and it is the first entry into each
            ## destination. But we do not check there is no extinction afterwards.
            ## The closest to the single Szendro LOD
            singlep <- pc$df
            singlep[, 1] <- as.character(singlep[, 1])
            singlep[, 2] <- as.character(singlep[, 2])
            singlep <- singlep[ do.call(order, singlep[, c(2, 3)]), ]
            singlep <- singlep[!duplicated(singlep[, 2]), ]
            gsingle <- igraph::graph_from_data_frame(singlep)
            lod_single <- igraph::all_simple_paths(gsingle, from = initMutant,
                                                   to = end, mode = "out")
            if(length(lod_single) != 1) stop("lod_single != 1")
        }
    }
    if(strict) {
        if(length(all_paths) > 1)
            stop("length(all_paths) > 1???")
        lodlist <- list(all_paths = NA,
                    lod_single = all_paths[[1]])
    } else {
        lodlist <- list(all_paths = all_paths,
                    lod_single = lod_single[[1]])
    }
    if(initMutant != "")
        attributes(lodlist)$initMutant <- initMutant
    return(lodlist)
}



## LOD_as_path_pre_2.9.2 <- function(llod) {
##     path_l <- function(u) {
##         if(length(u$lod_single) == 1) {
##             if(is.null(attributes(u)$initMutant))
##                 initMutant <- ""
##             else
##                 initMutant <- attributes(u)$initMutant
##             if(initMutant == "") initMutant <- "WT"
##             if(grepl("_is_end", u$lod_single))
##                 return(initMutant)
##             if(u$lod_single == "No_descendants")
##                 return(initMutant)
##         } else {
##             ## Deal with "" meaning WT
##             the_names <- names(u$lod_single)
##             the_names_wt <- which(the_names == "")
            
##             if(length(the_names_wt)) {
##                 if(length(the_names_wt) > 1) stop("more than 1 WT?!")
##                 if(the_names_wt > 1) stop("WT in position not 1?!")
##                 the_names[the_names_wt] <- "WT"
##             }
##             return(paste(the_names, collapse = " -> ")) 
##             ## return(paste0("WT", paste(names(u$lod_single),
##             ##                           collapse = " -> ")) )
##         }
##     }
##     if(identical(names(llod), c("all_paths", "lod_single")))
##         pathstr <- path_l(llod)
##     else {
##         ## should be a list
##         pathstr <- unlist(lapply(llod, path_l))
##     }
##     return(pathstr)
##     ## pathstr <- unlist(lapply(llod, function(x) paste(names(x$lod_single),
##     ##                                                  collapse = " -> ")))
##     ## return(paste0("WT", pathstr))
## }



LOD.oncosimul2_pre_2.9.2 <- function(x, strict = TRUE)
    return(LOD.internal_pre_2.9.2(x, strict))

## LOD.oncosimulpop_pre_2.9.2 <- function(x, strict = TRUE)
##     return(lapply(x, LOD.internal_pre_2.9.2, strict))






## Note for self: we could get all the LODs per simulation in the strict
## sense of those never becoming extinct if we subset the phylogClone
## object to children in which if we arrive at the children at any two
## times t and t+k, we retain only rows where any time > t is such that
## the popsize is > 0. But this is not worth it now.
