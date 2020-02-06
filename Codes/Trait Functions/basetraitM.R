

base_traitmatrix <- function(basetype = c("sphere", "cube"), ncols = 3, dsize = 100){
  
  #columns and rows of the trait matrix
  i <- ifelse(basetype =="cube", 1, 2)
  if(basetype == "sphere"){
    while(TRUE) {
      reference <- geozoo::sphere.solid.grid(p = ncols, n = i)$points
      if(nrow(reference) >= dsize  ) break
      i <- i + 1
    }
  } else {
    while(TRUE) {
      reference <- geozoo::cube.solid.grid(p = ncols, n = i)$points
      if(nrow(reference) >= dsize  ) break
      i <- i + 1
    }
    
  }
  
  return(reference)
}