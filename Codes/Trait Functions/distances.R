###############Distance and other small functions

#euclidean distance measure
eucldis <- function(x) sqrt(sum(x^2))

#L1 distance
L1_inf <- function(x) max(abs(x))

#Rescaling
rescale <- function(x) {
  apply(x, MARGIN = 2, function(y) {
    mn <- min(y);mx <- max(y);n <- length(y)
    a <- (n/(n-1)) * (mn - mx/n)
    b <- (n/(n-1)) * (mx - mn/n)
    (y - a)/(b - a)
  })
}

#function to retunr the point in the vector below the mean
eucldis2 <- function(x) {
  which(sign(x) == -1)
}

# function to decide which direction to shift the vector
ass <- function(x, shift = 0.1, dims = 2) {
  # x is a list
  #tbs is a vector of points closest to the mean
  
  if(length(x) == dims) {
    
    #all had -ve sign so shift both to the left
    delta <- rep(-shift, dims)
    
  } else if(length(x) < dims & length(x) != 0) {
    
    #at least one had -ve sign so shift only that to the left and the other to the right
    ff <- x == 1:dims
    delta <- numeric(dims)
    delta[ff == TRUE] <- -shift
    delta[ff != TRUE] <- shift
    
  } else {
    
    #all had +ve sign so shift both to the right
    delta <- rep(shift, dims)
  }
  return(delta)
}


d_eq <- function(delta) {
  
  y <- numeric(2)
  
  y[1] <- delta[1]^2 + delta[2]^2 - b^2
  y[2] <- (delta[1] * (x[1] - xbar[1])) + (delta[2] * (x[2] - xbar[2]) - b*d)
  y
}

d_jb <- function(delta) {
  
  n <- length(delta)
  Df <- matrix(numeric(n*n),n,n)
  Df[1,1] <- 2*delta[1]
  Df[1,2] <- 2*delta[2]
  Df[2,1] <- x[1] - xbar[1]
  Df[2,2] <- x[2] - xbar[2]
  
  Df
}
