## This file contains code from Youn and Simon, downloaded from
## http://linus.nci.nih.gov/Data/YounA/OrderMutation.zip

## on 2013-03-12.


## That is the code for their paper

## @article{,
## author = {Youn, Ahrim and Simon, Richard},
## doi = {10.1093/bioinformatics/bts168},
## issn = {1367-4811},
## journal = {Bioinformatics (Oxford, England)},
## month = jun,
## number = {12},
## pages = {1555--61},
## pmid = {22492649},
## title = {{Estimating the order of mutations during tumorigenesis from tumor genome sequencing data.}},
## url = {http://www.ncbi.nlm.nih.gov/pubmed/22492649},
## volume = {28},
## year = {2012}
## }


## We have included their code, contained in the file "function_library.r"
## verbatim here. The function used in our package, run.YounSimon, is like
## their function order_estimate, except we avoid using the MAX global
## variable, and pass it as an argument.
## Also modified gradient, main.function to pass MAX explicitly, and
## also hin and hin.jack and loglik (which don't use it)


## To avoid cluttering the NAMESPACE and prevent possible conflicts, we
## prepend ".YS_" to all the functions from this file.




.YS_hin <- function(x,AllP,division,sel,Y,a,aa, MAX)     # a vector function specifying inequality constraints such that hin[j] > 0 for all j
{                                                # used in the constrained optimization function auglag
    return(c(x,1-sum(x)))
}

.YS_hin.jac <- function(x,AllP,division,sel,Y,a,aa, MAX)   # Jacobian of hin, used in the constrained optimization function auglag
{
    return(rbind(diag(length(x)),rep(-1,length(x))))
}


.YS_subloglik <- function(P,N,a)
{
    temp <- rep(1,nrow(a[[N]]))
    for(j in 1:ncol(a[[N]]))
        temp <- temp*P[a[[N]][,j],j]

    return(log(sum(temp)))
}

.YS_loglik <- function(x,AllP,division,sel,Y,a,aa, MAX)     # negative log likelihood of the P_{k,i}, i=1,.. N-1. k=1,..K. where N= # of driver genes.
{
    AllP[,division[[sel]]] <- c(x,1-sum(x))     # P_{k,N}=1-\sum_{i=1}^{N-1}P_{k,i} is used in calculating the log likelihood, where N is the number of driver genes.

    loglik <- rep(0,nrow(Y))
    for(i in 1:nrow(Y))
    {
        N <- sum(Y[i,])
        P <- as.matrix(AllP[rep(1:ncol(Y),Y[i,]),1:N])
        loglik[i] <- .YS_subloglik(P,N,a)
    }
    return(-sum(loglik))
}

.YS_gradient <- function(x,AllP,division,sel,Y,a,aa, MAX)   # gradient of the negative log likelihood with respect to the P_{k,i}, i=1,.. N-1. k=1,..K. where N= # of driver genes.
{
    grad <- matrix(0,nr=ncol(Y),nc=ncol(AllP))
    AllP[,division[[sel]]] <- c(x,1-sum(x))

    loglik <- rep(0,nrow(Y))
    for(i in 1:nrow(Y))
    {
        N <- sum(Y[i,])
        P <- AllP[rep(1:ncol(Y),Y[i,]),1:N]
        if(min(division[[sel]])<=N)
        {
            if(N>1)       # for samples with more than one mutation
            {
                loglik[i] <- .YS_subloglik(P,N,a)

                count <- 0
                for(k in rep(1:ncol(Y),Y[i,]))
                {
                    count <- count+1
                    for(l in division[[sel]][division[[sel]]<=N])
                    {
                        tempP <- as.matrix(P[-count,-l])
                        if(l>MAX)
                            grad[k,l] <- grad[k,l]+exp( .YS_subloglik(tempP,N-1,a) -loglik[i])/(N-MAX)
                        else
                            grad[k,l] <- grad[k,l]+exp( .YS_subloglik(tempP,N-1,aa) -loglik[i])
                    }
                }
                grad <- t(t(grad)-grad[nrow(grad),])       # This is because P_{k,N}=1-\sum_{i=1}^{N-1}P_{k,i} is used in calculating the log likelihood
            }
            else       # for samples with one mutation
            {
                k <- which(Y[i,]==1)
                grad[k,1] <- grad[k,1]+1/P
                grad <- t(t(grad)-grad[nrow(grad),])    # This is because P_{k,N}=1-\sum_{i=1}^{N-1}P_{k,i} is used in calculating the log likelihood
            }
        }
    }

    grad <- grad[1:(nrow(grad)-1),]
    grad <- rowSums(as.matrix(grad[,division[[sel]]]))
    return(-grad)

}




.YS_main.function <- function(Y,a,aa, MAX){
    ## Modified by RDU to pass MAX explicitly
    nomut <- max(rowSums(Y))
    no.gene <- ncol(Y)

    nparam <- min(MAX+1,nomut)
    division <- vector("list",nparam)
    for(i in 1:(nparam-1)) division[[i]] <- i
    division[[nparam]] <- nparam:nomut

###############################################################

    AllPtrue <- matrix(0,nr=ncol(Y),nc=nomut)


    for(j in 1:nparam)
    {
        init=runif(no.gene); init=init/sum(init)
        AllPtrue[,division[[j]]] <-  init
    }

    prev0=Inf
    prev=Inf
    n=0
    decrease=1

    no.repeat=20
    decrease.limit=10^(-6)

    while(decrease>decrease.limit & n<=no.repeat)
    {
        n=n+1
        prev0=prev              # For each k, the length of \vec{P_{k}} optimized is N-1 where N is the number of driver genes. This is because \sum_{i=1}^{N} P_{k,i} = 1 and the value of P_{k,N} is determeined by the other N-1 P_{k,i}, thus we optimize only for P_{k,i} for i=1,..N-1. we let P_{k,N}=1-\sum_{i=1}^{N-1}{P_{k,i}}
        for(sel in 1:nparam)    # find \vec{P_{k}} maximizing the log likelihood  for each k in turn with the constraint that 0<P_{k,i} for i=1,..N-1 and 0< 1-\sum_{i=1}^{N-1}{P_{k,i}}
        {
            init=runif(no.gene); init=init/sum(init)  # initial values for \vec{P_{k}}
                                                     # auglag is a constrained optimization function in R
            x <- auglag(par=init[-no.gene],fn=.YS_loglik,gr=.YS_gradient,hin=.YS_hin,hin.jac=.YS_hin.jac,control.optim=list(maxit=10000),control.outer=list(trace=F),AllP=AllPtrue,division=division,sel=sel,Y=Y,a=a,aa=aa, MAX=MAX)
            tempAllPtrue <- AllPtrue
            nonnegative <- pmax(c(x$par,1-sum(x$par)),0)
            tempAllPtrue[,division[[sel]]] <- nonnegative/sum(nonnegative)   # This is needed since sometimes the value x$par are negative.
            fval <- .YS_loglik(tempAllPtrue[-no.gene,1],tempAllPtrue,division,1,Y,a,aa)     # negative log likelihood

            decrease=prev-fval
            if(decrease>0)          # if the negative log likelihood decreases, then change old P_{k} with new P_{k}
            {
                prev=fval
                AllPtrue <- tempAllPtrue
            }
        }
        decrease=prev0-prev
    }

    return(list(AllPtrue,prev))
}


.YS_generatea <- function(MAX,MAX.mut)
{
a <- vector("list",MAX.mut)
a[[1]] <- matrix(1,nr=1,nc=1)

for(N in 2:MAX.mut)
{
    x <- matrix(1:N,nr=N)
    for(i in 1:(min(N,MAX)-1))
    {
        y <- rep(1:N,nrow(x))[-as.vector(t(x+((1:nrow(x))-1)*N))]
        x <- t(matrix(as.vector(t(matrix(rep(as.vector(x),N-ncol(x)),nr=nrow(x)))),nr=ncol(x)))
        x <- cbind(x,y)
    }
    if(N>MAX)
    {
        y <- rep(1:N,nrow(x))[-as.vector(t(x+((1:nrow(x))-1)*N))]
        y <- t(matrix(y,nc=nrow(x)))
        x <- cbind(x,y)
    }
    a[[N]] <- x
}
return(a)
}



.YS_order_estimate <- function(sample.gene.mutation,N,parallel)
{ ## in the package, we use run.YounSimon, not this.
    a <- .YS_generatea(MAX,max(rowSums(sample.gene.mutation)))
    aa <- .YS_generatea(MAX-1,max(rowSums(sample.gene.mutation)))

## Run optimization with different inital values N times using parallel computing if the variable "parallel" is TRUE :
    if(parallel)
        tmp <- foreach (kk = 1:N) %dopar%  .YS_main.function(sample.gene.mutation,a,aa)
    else ## Otherwise, run for loop
    {
        tmp <- vector("list", N)
        for(i in 1:N) tmp[[i]] <- .YS_main.function(sample.gene.mutation,a,aa)
    }

    minusloglik=rep(Inf,N)
    for(l in 1:N)
        if(is.list(tmp[[l]])) minusloglik[l]=tmp[[l]][[2]]

    result <- tmp[[which(minusloglik==min(minusloglik))[1]]][[1]] #find the one giving the maximum likelihood

    return(result)

}



.YS_total <- function(X) # Prob that there is any mutation in a gene
    return(1-apply(1-X,1,prod))


.YS_Intersection <- function(X,N)
{
    sel <- combinations(ncol(X),N)
    res <- 0
    for(i in 1:nrow(sel))
        res <- res+apply(X[,sel[i,]],1,prod)
    return(res)
}
.YS_Union <- function(X)
{
    res <- 0
    for(i in 2:ncol(X))
        res <- res+.YS_Intersection(X,i)*(-1)^(i+1)
    return(apply(X,1,sum)+res)
}

# constructing BCa CI from bootstrapped data for each of the colon and lung data
.YS_ConfInterval <- function(x,mle,jack,sample.gene.mutation,lower,upper)
{
    no.gene=ncol(sample.gene.mutation)

    z0 <- matrix(0,nr=no.gene,nc=ncol(mle))       # bias-correction constant

    for(i in 1:no.gene)
    {
        for(j in 1:ncol(mle))
        {
            z0[i,j] <- qnorm(mean(x[i,j,]<=mle[i,j],na.rm=T))
        }
    }

    z0[z0==Inf]=1000

    theta.mean <- apply(jack,1:2,mean,na.rm=T)

    num <- denom <- 0

    for(i in 1:nrow(sample.gene.mutation))
    {
        if(sum(is.na(jack[,,i]))==0)
        {
            num <- num+(theta.mean-jack[,,i])^3
            denom <- denom+(theta.mean-jack[,,i])^2
        }
    }
    accel <- num/6/denom^1.5   # acceleration constant

    alpha1 <- pnorm(z0+(z0+qnorm(lower))/(1-accel*(z0+qnorm(lower))))
    alpha2 <- pnorm(z0+(z0+qnorm(upper))/(1-accel*(z0+qnorm(upper))))

    CI=array(0,dim=c(no.gene,2,ncol(mle)))

    for(i in 1:no.gene)
    {
        for(j in 1:ncol(mle))
        {
            CI[i,1,j] <- quantile(x[i,j,],alpha1[i,j],na.rm=T)
            CI[i,2,j] <- quantile(x[i,j,],alpha2[i,j],na.rm=T)
        }
    }
    CI=round(CI,2)
    a=round(mle,2)

    tab =NULL
    for(i in 1:ncol(a))
    {
        tab=cbind(tab,"&",a[,i],"&","(" ,CI[,1,i], ",",CI[,2,i],")")
    }
    tab=cbind(colnames(sample.gene.mutation),tab,"\\")   # for generating tables in the paper

    return(tab)
}
