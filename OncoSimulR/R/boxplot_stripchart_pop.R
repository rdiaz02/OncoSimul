## Authors: Sara Dorado Alfaro
##          Miguel Hern?ndez del Valle,
##          ?lvaro Huertas Garc?a,
##          Diego Ma?anes Cayero,
##          Alejandro Martin Mu?oz, 

## Date: 16/01/2020

## Summary: This script contains functions created for plotting the results
## from several iteration of oncoSimulPop and modified functions from 
## Oncosimul.R for placing the legend outside the plot. 

#############################################################################
############## Functions for plotting several simulations ###################

#############################################################################
## Create box-plot, title and axis parameters
## Plot box plot
simul_boxplot2 <- function(df, main,  xlab,
                           ylab, colors) {
  ## Create box plot, title and axis parameters
  e <- ggplot(df, aes(x = Genotype, y = N)) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11)
    )
  
  ## No title
  if (main == FALSE) {
    e + geom_boxplot(aes(fill = Genotype)) +
      ## Show mean
      stat_summary(fun.y = mean, geom = "point", shape = 18, 
                   size = 2.5, color = "#FC4E07") +
      ## x and y axis label
      xlab(xlab) + ylab(ylab) +
      scale_fill_manual(values = colors) +
      stat_summary(fun.y = mean, geom = "point",
                   shape = 18, size = 2.5, color = "#FC4E07") +
      xlab(xlab) + ylab(ylab) + scale_fill_manual(values = colors)
  }
  ## Title
  else {
    e + geom_boxplot(aes(fill = Genotype)) +
      ## Show mean
      stat_summary(fun.y = mean, geom = "point",
                   shape = 18, size = 2.5, color = "#FC4E07") +
      labs(title = main) +
      xlab(xlab) + ylab(ylab) + scale_fill_manual(values = colors)
  }
}

## Extract and create a data frame with results from several simulations
compositionPop2 <- function(objPop, cols = NULL,  xlab = "Genotype",
                            ylab = "N", main = FALSE) {
  ## Extract the information to create a data frame
  clon_labels <- c("WT", objPop[[1]]$geneNames)
  n_labels <- length(clon_labels)
  listPop <- vapply(objPop, function(x) tail(x[[1]], 1)[1, -1], 
                    as.double(1:n_labels))
  dfPop <- data.frame("Genotype" = rep(clon_labels, length(listPop)/n_labels),
                      "N" = c(listPop))
  
  ## Colors of boxplot
  ## Extract the colors generate by myhsvcols
  if (is.null(cols)) {
    ndr <- 1:n_labels
    y <- objPop[[1]]$pops.by.time[, 2:ncol(objPop[[1]]$pops.by.time), 
                                  drop = FALSE]
    ymax <- colSums(y)
    cols <- OncoSimulR:::myhsvcols(ndr, ymax)
    # Reorder the colors
    cols <- c(cols$colors[-1], cols$colors[1])
    ## Get the colors given by the user
  } else {
    if (length(cols) != n_labels)
      stop("The number of colors is not equal to the number of items")
  }
  
  simul_boxplot2(dfPop, colors = cols, main = main,  xlab = xlab, ylab = ylab)
}

###############################################################################
## Stripchart
stripChartPop <- function(dfPop, ylab, ...) {
  stripchart(dfPop, vertical = TRUE, ylab = ylab, ...)
  f1 <- function(x, num_genotypes) {
    ## Draw a segment between genotype 1 and genotype 2
    num_genotypes <- length(x)
    i <- 1
    while (num_genotypes > i) {
      segments(x0 = i, x1 = i+1,
               y0 = x[i],
               y1 = x[i+1],
               col = rainbow(5))
      i <- i+1
    }
  }
  ## Read data frame by rows (simulation by simulation)
  apply(dfPop, 1, f1)
}


##  Plot the data as points and join with lines the ones that come from the same 
##  simulation.
meanCompositionPop <- function(objPop, ylab = "N", ...) {
  condi <- c("WT", objPop[[1]]$geneNames)
  ## Extract the information.
  ## $pops.by.time contains all the results
  ## $pops.by.time contains the times at wich results are taken
  ## see length of times and select results from the half time to the end
  ## Calculate the mean for each Genotype in each simulation
  ## Use lapply because results can be not rectangular
  listPop <- lapply(objPop, function(x) 
    (colMeans(tail(x$pops.by.time, length(x$pops.by.time[,1])/2)[,-1])))
  ## Create data frame with the means Genotype from a list
  dfPop <- data.frame(matrix(unlist(listPop), 
                             ncol = length(condi), byrow = TRUE))
  colnames(dfPop) <- condi
  stripChartPop(dfPop, ylab = ylab, ...)
  dfPop
}