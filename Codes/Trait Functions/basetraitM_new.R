
#source("distances.R")

base_traitmatrix_new <- function(x, gen_matrix, rescale = F, parallel = F) {
  
  #rescale to 0-1 range
  if(rescale == F) {
    x_rescale <- x 
  } else {
    x_rescale <- rescale(x)
  }
  
  #compute convex hull of the observed trait matrix
  conv_hull <- geometry::convhulln(x_rescale, option = "FA")
  #points making up the convex_hull
  mconv_hull <- as.data.frame(unique(x_rescale[conv_hull$hull, ]))
  #extract the area of the hull
  conv_area <- conv_hull$area
  
  #check which of the generated points is within the convex-hull
  if(parallel == T) {
    ncl <- parallel::makeCluster(parallel::detectCores() - 1)
    
    #parallel::clusterExport(ncl, c("mconv_hull", "gen_matrix", "conv_area"))
    parallel::clusterEvalQ(ncl, library(geometry))
    
    tests <- unlist(parallel::clusterApply(ncl, 1:nrow(gen_matrix), function(i) {
      new_area <- convhulln(rbind(mconv_hull, gen_matrix[i, ]), "FA")$area
      return(ifelse(new_area == conv_area, 1, 0))
    }))
    parallel::stopCluster(ncl)
  } else {
    tests <- sapply(1:nrow(gen_matrix), function(i) {
      
      new_area <- geometry::convhulln(rbind(mconv_hull, gen_matrix[i, ]), "FA")$area
      return(ifelse(new_area == conv_area, 1, 0))
      
    })
  }
  
  test_data <- gen_matrix[which(tests == 1), ]
  #if(nrow(test_data) <= 100) test_data <- gen_matrix else test_data <- test_data
  
  return(test_data)
}
