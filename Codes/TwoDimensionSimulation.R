library(magrittr)
library(mvtnorm)

##All needed functions are stored in the Trait functions folder

source("Trait Functions/distances.R")
source("Trait Functions/TOP.R")
source("Trait Functions/TED.R")

source("Trait Functions/basetraitM.R")
source("Trait Functions/basetraitM_new.R")



#############Different Sample Size
set.seed(1100)
n_dim <- 2
n_row <- c(50, 100, 200, seq(500, 3000, by = 500))

#base matrix for modified TED
#generate a large uniformly distributed sample space
gen_matrix <- matrix(NA, nrow = 5000, ncol = 2)
gen_matrix <- apply(gen_matrix, 2, function(x) runif(nrow(gen_matrix), 0, 1))

#setting up the parameters of the normal distribution
#snow::clusterEvalQ(ncl, library(mvtnorm))

mu_start <- c(5, 10)
delta_mu <- seq(0, 10, by = 0.5)
rho_start <- 0.007

#equal variance
sig_11_start <- 0.85
cv_start <- sig_11_start*sig_11_start*rho_start

#snow::clusterExport(ncl, c("n_dim", "n_row", "b_matrix2", "mu_start", "delta_mu", 
#                           "rho_start", "sig_11_start", "cv_start"))

Results_comm <- vector("list", length(delta_mu))
Results_samp <- vector("list", length(n_row))


# Secnario One (checking indicies reaction to changes in community structure)

for(i in 1:length(delta_mu)) {
  mu_change <- c(mu_start[1] + delta_mu[i], mu_start[2] + delta_mu[i])
  #parallel::clusterExport(ncl, c("i", "mu_change"))
  print(i)
  #5 replicates per sample size/chage in community structure
  replicates <- lapply(1:5, function(j) {
    ddata <- rmvnorm(800, mean = mu_start, sigma = matrix(c(sig_11_start, 
                                                           cv_start,cv_start, 
                                                           sig_11_start), 
                                                  nrow = 2, ncol = 2))
    if(i == 1 ) {
      ddat <- ddata
    } else {
      ddat <- ddata
      tobe_shifted <- sample(1:nrow(ddata), 0.20*nrow(ddat))
      ddat[tobe_shifted, ] <- ddat[tobe_shifted, ] - mu_start + mu_change
    }
    
    row.names(ddat) <- row.names(ddata) <- 1:nrow(ddat)
    
    #FRIC; FEve
    FDS_comm <- FD::dbFD(ddat, stand.x = F, calc.FRic = T, stand.FRic = F)
    message("Done with FDs, starting with TOP and TED")
    
    #TOP and modified TOP
    tops_comm <- TOP(ddat)
    message("Done with TOP starting with TED")
    
    ## change in community structure
    #TED
    teds_comm1 <- TED(ddat, rescale = T, base_matrix = b_matrix2,
                      modification = "None")
    #modified TED
    #base matrix for classic TED
    b_matrix2 <- base_traitmatrix(basetype = "cube", ncols = n_dim, 
                                  dsize = nrow(ddat))
    b_matrix1 <- base_traitmatrix_new(ddat, gen_matrix = gen_matrix, 
                                      rescale = T, parallel = F)
    teds_comm2 <- TED(ddat, rescale = T, base_matrix = b_matrix1,
                      modification = "None")
    
    Comm_Change <- data.frame(FRIC = FDS_comm$FRic, TOP = tops_comm$TA, 
                             TOPM = tops_comm$TA2, 
                             FEVE = FDS_comm$FEve, TED = teds_comm1, TEDM = teds_comm2, 
                             RAO = FDS_comm$RaoQ, FDIS = FDS_comm$FDis)
    message("Done with TED, next replicate commences")
    return(Comm_Change = Comm_Change)
  })
  
  Results_comm[[i]] <- replicates %>% do.call(rbind.data.frame, .)
}

#### Secnario Two
########## uniform Distribution #############
set.seed(1100)

delete <- seq(0.10, 0.50, by = 0.10)#seq(0.10, 0.50, by = 0.10)
### list for storing the results
Results_comm22 <- vector("list", length(delete))

for(i in 1:(length(delete) + 1)) {
  print(i)
  replicates_richness <- lapply(1:5, function(j) {
  Tsdata1 <- matrix(NA, ncol = 2, nrow = 1000)
  Tsdata1 <-  apply(Tsdata1, 2, function(x) runif(1000, 0, 1))
  #middle and distance to the center
  xbar <- apply(Tsdata1, 2, mean)
  dist_mean <- apply(Tsdata1 , 1, function(x) {
    eucldis(x - xbar)
  })
  
  if(i == 1) {
    
    remaining <- Tsdata1
    row.names(remaining) <- 1:nrow(remaining)
    #FRIC and FEVE
    FDS_comm <- FD::dbFD(remaining, stand.x = F, calc.FRic = T, stand.FRic = F, 
                         calc.CWM = F, calc.FGR = F, calc.FDiv = T)
    message("done with that")
    
    #TOP and modified TOP
    tops_comm <- TOP(remaining)
    message("Done with TOP starting with TED")
    
    ## change in community structure
    #TED
    #base matrix for classic TED
    b_matrix2 <- base_traitmatrix(basetype = "cube", ncols = n_dim, 
                                  dsize = nrow(remaining))
    teds_comm1 <- TED(remaining, rescale = T, base_matrix = b_matrix2,
                      modification = "None")
    #modified TED
    b_matrix1 <- base_traitmatrix_new(remaining, gen_matrix = gen_matrix, 
                                      rescale = T, parallel = F)
    teds_comm2 <- TED(remaining, rescale = F, base_matrix = b_matrix1,
                      modification = "None")
    
  } else {
    
    #middle <- mrfDepth::hdepthmedian(Tsdata1, maxdir = 500)$median
    tobe_sampled <- which(dist_mean <= quantile(dist_mean, probs = delete[i-1]))
    deltas <- sapply(tobe_sampled, dts <- function(j, b, xbar) {
      
      environment(d_eq) <- environment()
      environment(d_jb) <- environment()
      
      x <- Tsdata1[j, ]
      b <- b - dist_mean[j]
      d <- dist_mean[j]
      xbar <- xbar
      
      delta_start <- c(b^2, b^2)
      
      ds <- nleqslv::nleqslv(x = delta_start,
                             fn = d_eq,
                             jac = d_jb,
                             control=list(btol=.01, delta="newton")
      )
      ds$x
      
    }, b = min(dist_mean[-tobe_sampled]), xbar = xbar)
    
    
    Tsdata12 <- Tsdata1
    Tsdata12[tobe_sampled, ] <- Tsdata1[tobe_sampled, ] + t(deltas)
    plot(Tsdata12)
    
    row.names(Tsdata12) <- 1:nrow(Tsdata12)
    #FRIC and FEVE
    FDS_comm <- FD::dbFD(Tsdata12, stand.x = F, calc.FRic = T, stand.FRic = F, 
                         calc.CWM = F, calc.FGR = F, calc.FDiv = T)
    message("done with that")
    
    #TOP and modified TOP
    tops_comm <- TOP(Tsdata12)
    message("Done with TOP starting with TED")
    
    ## change in community structure
    #TED
    b_matrix2 <- base_traitmatrix(basetype = "cube", ncols = n_dim, 
                                  dsize = nrow(Tsdata12))
    teds_comm1 <- TED(Tsdata12, rescale = T, base_matrix = b_matrix2,
                      modification = "None")
    #modified TED
    b_matrix1 <- base_traitmatrix_new(Tsdata12, gen_matrix = gen_matrix, 
                                      rescale = T, parallel = F)
    teds_comm2 <- TED(Tsdata12, rescale = F, base_matrix = b_matrix1,
                      modification = "None")
    
  }
  
  Comm_Change <- data.frame(FRIC = FDS_comm$FRic, TOP = tops_comm$TA, 
                            TOPM = tops_comm$TA2, 
                            FEVE = FDS_comm$FEve, TED = teds_comm1, TEDM = teds_comm2, 
                            RAO = FDS_comm$RaoQ, FDIS = FDS_comm$FDis)
  message("Done with TED, next replicate commences")
  return(Comm_Change = Comm_Change)
 })
  Results_comm22[[i]] <-  replicates_richness %>% do.call(rbind.data.frame, .)#Comm_Change
}
save.image("TwoDimensionSimulationResults_results1and2.RData")


#####change in sample size
for(i in 1:length(n_row)) {
  print(i)
    replicates_samp <- lapply(1:5, function(j) {
    
      ddat2 <- rmvnorm(n_row[i], mean = mu_start, sigma = matrix(c(sig_11_start, 
                                                                 cv_start,cv_start, 
                                                                 sig_11_start), 
                                                               nrow = 2, ncol = 2))
      row.names(ddat2) <- 1:nrow(ddat2)
      FDS_sample <- FD::dbFD(ddat2, stand.x = F, calc.FRic = T, stand.FRic = F)
      message("Done with FDS, next is TOP")
      
      #TOP
      tops_sample <- TOP(ddat2)
      message("Done with TOP, next is TED")
      
      #TED
      b_matrix2 <- base_traitmatrix(basetype = "cube", ncols = n_dim, 
                                    dsize = nrow(ddat2))
      teds_sample1 <- TED(ddat2, rescale = T, base_matrix = b_matrix2,
                          modification = "None")
      #modified TED
      b_matrix3 <- base_traitmatrix_new(ddat2, rescale = T, 
                                        gen_matrix = gen_matrix, parallel = F)
      teds_sample2 <- TED(ddat2, rescale = T, base_matrix = b_matrix3,
                          modification = "None")
      message("Done with TED, next replicate commences")
      Sample_Change <- data.frame(FRIC = FDS_sample$FRic, TOP = tops_sample$TA, 
                               TOPM = tops_sample$TA2, 
                               FEVE = FDS_sample$FEve, 
                               TED = teds_sample1, TEDM = teds_sample2, 
                               RAO = FDS_sample$RaoQ, FDIS = FDS_sample$FDis)
    return(Sample_Change = Sample_Change)
  })
  Results_samp[[i]] <- replicates_samp %>% do.call(rbind.data.frame, .)
}


save.image("TwoDimensionSimulation_SampleSize_new_Results.RData")
