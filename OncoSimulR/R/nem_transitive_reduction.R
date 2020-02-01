###############################################
## <code from nem package>

## Title: (Dynamic) Nested Effects Models and Deterministic Effects
##         Propagation Networks to reconstruct phenotypic hierarchies
## Version: 2.60.0
## Author: Holger Froehlich, Florian Markowetz, Achim Tresch, Theresa
## Niederberger, Christian Bender, Matthias Maneck, Claudio Lottaz, Tim
## Beissbarth
## Maintainer: Holger Froehlich <frohlich at bit.uni-bonn.de>
## https://bioconductor.org/packages/release/bioc/html/nem.html
## License: GPL (>= 2)

nem_transitive.reduction <- function(g){
	if (!any(class(g)%in%c("matrix","graphNEL"))) stop("Input must be an adjacency matrix or graphNEL object")
	if(any(class(g) == "graphNEL")){
		g = as(g, "matrix")		
	}
# 	if("Rglpk" %in% loadedNamespaces()){ # constraints müssen einzeln hinzugefügt und jedesmal das ILP gelöst werden. Danach müssen jedesmal die constraints überprüft werden und nicht mehr gebrauchte rausgeschmissen werden
# 		g = abs(g) - diag(diag(g))
# 		mat = matrix(0, ncol=sum(g), nrow=0)
# 		idx = cbind(which(g == 1, arr.ind=T), 1:sum(g))		
# 		for(y in 1:nrow(g)){
# 			for(x in 1:nrow(g)){
# 				if(g[x,y] != 0){
# 					for(j in 1:nrow(g)){
# 						if((g[y,j] != 0) && (g[x,j] != 0)){
# 							mat.tmp = double(sum(g))
# 							mat.tmp[idx[idx[,1] == x & idx[,2] == y, 3]] = 1
# 							mat.tmp[idx[idx[,1] == y & idx[,2] == j, 3]] = 1
# 							mat.tmp[idx[idx[,1] == x & idx[,2] == j, 3]] = -1
# 							mat = rbind(mat, mat.tmp)
# 						}						
# 					}
# 				}
# 			}
# 		}
# 		
# 		solve.problem = function(mat.tmp){
# 			obj = rep(1, NCOL(mat.tmp))
# 			rhs = 2
# 			dir = ">="
# 			sol = Rglpk_solve_LP(obj, mat.tmp, dir, rhs, max=TRUE, types=rep("B", NCOL(mat.tmp)))
# 			print(sol)			
# 			del = idx[which(sol$solution == 0), c(1, 2),drop=F]
# 			for(i in 1:NROW(del)){
# 				g[del[i,1], del[i,2]] = 0
# 			}
# 			g
# 		}
# 		
# 		while(NROW(mat) > 0){
# 			g = solve.problem(mat[1,,drop=F])
# 			i = 1
# 			while(i <= NROW(mat)){
# 				idx.pos = which(mat[i,,drop=F] == 1)
# 				idx.neg = which(mat[i,,drop=F] == -1)
# 				xy = idx[idx[,3] %in% idx.pos, c(1,2)]
# 				z = idx[idx[,3] == idx.neg, c(1,2)]
# 				if(!(g[xy[1,1], xy[1,2]] == 1 & g[xy[2,1], xy[2,2]] == 1 & g[z[1], z[2]] == 1)) # remove resolved constraints
# 					mat = mat[-i,,drop=F]
# 				else
# 					i = i + 1
# 			}
# 		}
# 	}
# 	else{
# 	if(class(g) == "matrix"){		
		# modified algorithm from Sedgewick book: just remove transitive edges instead of inserting them		
    g = nem_transitive.closure(g, mat=TRUE) # bug fix: algorithm only works for transitively closed graphs!
		g = g - diag(diag(g))
		type = (g > 1)*1 - (g < 0)*1	
		for(y in 1:nrow(g)){
			for(x in 1:nrow(g)){
				if(g[x,y] != 0){
					for(j in 1:nrow(g)){
						if((g[y,j] != 0) & sign(type[x,j])*sign(type[x,y])*sign(type[y,j]) != -1){ 
						    g[x,j] = 0
						}
					}
				}
			}
		}
# 	}
# 	else{		
# 		nodenames=nodes(g)		
# 		for(y in 1:length(nodes(g))){
# 			edges = edgeL(g)
# 			x = which(sapply(edges,function(l) y %in% unlist(l)))
# 			j = unlist(edges[[y]])
# 			cands = sapply(edges[x], function(e) list(intersect(unlist(e),j)))			
# 			cands = cands[sapply(cands,length) > 0]
# 			if(length(cands) > 0)
# 				for(c in 1:length(cands)){ 				
# 					jj = unlist(cands[c])					
# 					g = removeEdge(rep(names(cands)[c],length(jj)),nodenames[jj],g)
# 				}
# 		}
# 	}	
	g		
}



nem_transitive.closure <- function(g,mat=FALSE,loops=TRUE){

    if (!any(class(g)%in%c("graphNEL","matrix"))) stop("Input must be either graphNEL object or adjacency matrix")
    g <- as(g, "matrix")
    
    #-- adjacency matrix
#     if (class(g)=="matrix"){
		n <- ncol(g)
		matExpIterativ <- function(x,pow,y=x,z=x,i=1) {
		while(i < pow) {
			z <- z %*% x
			y <- y+z
			i <- i+1
		}
		return(y)
		}
	
		h <- matExpIterativ(g,n)
		h <- (h>0)*1   
		dimnames(h) <- dimnames(g)
		if (!loops) diag(h) <- rep(0,n) else diag(h) <- rep(1,n)
		if (!mat) h <- as(h,"graphNEL")	
#     }

# #     -- graphNEL object
#     if (class(g)=="graphNEL"){
#         tc <- RBGL::transitive.closure(g)    
#         if (loops) tc$edges <- unique(cbind(tc$edges,rbind(tc$nodes,tc$nodes)),MARGIN=2)
# 
#         h <- ftM2graphNEL(ft=t(tc$edges),V=tc$nodes)
#         if (mat) h <- as(h, "matrix")
#     } 
       
    return(h)
}
