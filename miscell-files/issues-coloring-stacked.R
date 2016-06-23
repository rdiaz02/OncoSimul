## Go to epist dir, use the 56 example and do

plot(p1[[1]], type = "stacked", xlim = c(0, 2000), ylim = c(1900, 2100))

## and you will see like there is red when it shouldn't. Is just a slight
## overimposition.

## The problem is in the plot of the polygon as sometimes you get a
## vertical line where none should be. It is an R issue.

## Solutions?
##  no, you cannot pass a NA for color at that point.
##  have polygons ONLY where there is something. But that is tough.
## So use a black border line with width of 1, and this does not show.
