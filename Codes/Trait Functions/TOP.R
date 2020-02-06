#######TOP Function
TOP <- function(trait_matrix) {
  n <- nrow(trait_matrix)
  area <- c(); row_num <- c()
  while(is.null(trait_matrix) == F & 
        (nrow(trait_matrix) > dim(trait_matrix)[2])) {
    hulls <- geometry::convhulln(trait_matrix, "FA")
    area <- c(area, hulls$area)
    mat_reduce <- t(hulls$hull)
    row_num <- c(row_num, nrow(unique(trait_matrix[mat_reduce,])))
    trait_matrix <- trait_matrix[-mat_reduce, ]
    if(is.null(nrow(trait_matrix)) == T) break
  }
  
  TA <-  sum(area)
  TA2 <- TA/n
  TA3 <- mean(area)
  return(list(TA = TA, Areas = area, N_points = row_num, 
              TA2 = TA2, TA3 = TA3))
}