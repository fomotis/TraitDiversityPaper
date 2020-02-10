#vsc_home
#vsc_home <- Sys.getenv("VSC_HOME")

#UNamur CECI HOME
#ceci_home <- Sys.getenv("HOME")
#ncore <- parallel::detectCores(as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE")))
#ncl <- parallel::makeCluster(ncore)

###distances
##VSC
#source(paste0(vsc_home,"/Functions/distances.R"))

#ncore <- parallel::detectCores() - 1
#ncl <- parallel::makeCluster(ncore)
#snow::clusterEvalQ(ncl, source(paste0("Codes","/Functions/distances.R")))

##UNamur CECI
#source(paste0(ceci_home,"/Functions/distances.R"))
#snow::clusterEvalQ(ncl, source(paste0(ceci_home,"/Functions/distances.R")))

###TOP
##VSC
#source(paste0(vsc_home,"/Functions/TOP.R"))
#source(paste0(vsc_home,"/Functions/TOP_NEW.R"))

##My system

#source(paste0("Codes","/Functions/TOP_NEW.R"))
#snow::clusterEvalQ(ncl, source(paste0("Codes","/Functions/TOP.R")))
#snow::clusterEvalQ(ncl, source(paste0("Codes","/Functions/TOP_NEW.R")))

##UNamur CECI
#snow::clusterEvalQ(ncl, source(paste0(ceci_home,"/Functions/TOP.R")))
#snow::clusterEvalQ(ncl, source(paste0(ceci_home,"/Functions/TOP_NEW.R")))

###TED
##VSC
#source(paste0(vsc_home,"/Functions/TED.R"))
#source(paste0(vsc_home,"/Functions/TED_NEW.R"))
#snow::clusterEvalQ(ncl, source(paste0("Codes","/Functions/basetraitM.R")))
#snow::clusterEvalQ(ncl, source(paste0("Codes","/Functions/basetraitM_new.R")))
#snow::clusterEvalQ(ncl, source(paste0("Codes","/Functions/TED.R")))

##UNamur CECI
#snow::clusterEvalQ(ncl, source(paste0(ceci_home,"/Functions/basetraitM.R")))
#snow::clusterEvalQ(ncl, source(paste0(ceci_home,"/Functions/basetraitM_new.R")))
#snow::clusterEvalQ(ncl, source(paste0(ceci_home,"/Functions/TED.R")))

#source(paste0("Codes","/Functions/TED_NEW.R"))


#normal
# set.seed(1992)
# Tsdata1n <- rmvnorm(1000, mean = c(5, 5), sigma = matrix(c(sig_11_start, 
#                                                            cv_start,cv_start, 
#                                                            sig_11_start), 
#                                                          nrow = 2, ncol = 2)) 
# 
# xbar_n <- apply(Tsdata1n, 2, mean)
# dist_mean_n <- apply(Tsdata1n , 1, function(x) {
#   eucldis(x - xbar_n)
# })
# 
# delete <- seq(0.10, 0.50, by = 0.10)
# 
# Results_comm2 <- vector("list", length(delete))
# 
# 
# for(i in 1:(length(delete) + 1)) {
#       print(i)
#   #replicates_richness <- lapply(1:5, function(j) {
#     
#         if(i == 1) {
#           
#             remaining <- Tsdata1n
#             row.names(remaining) <- 1:nrow(remaining)
#             #FRIC and FEVE
#             FDS_comm <- FD::dbFD(remaining, stand.x = F, calc.FRic = T, stand.FRic = F, 
#                                  calc.CWM = F, calc.FGR = F, calc.FDiv = T)
#             message("done with that")
#             
#             #TOP and modified TOP
#             tops_comm <- TOP(remaining)
#             message("Done with TOP starting with TED")
#             
#             ## change in community structure
#             #TED
#             teds_comm1 <- TED(remaining, rescale = T, base_matrix = b_matrix2,
#                               modification = "None")
#             #modified TED
#             b_matrix1 <- base_traitmatrix_new(remaining, gen_matrix = gen_matrix, 
#                                               rescale = T, parallel = F)
#             teds_comm2 <- TED(remaining, rescale = F, base_matrix = b_matrix1,
#                               modification = "None")
#             
#         } else {
#           
#           tobe_sampled <- which(dist_mean_n <= quantile(dist_mean_n, probs = delete[i-1]))
#           deltas <- sapply(tobe_sampled, dts <- function(j, b, xbar_n) {
#             
#             environment(d_eq) <- environment()
#             environment(d_jb) <- environment()
#             
#             x <- Tsdata1n[j, ]
#             b <- b - dist_mean_n[j]
#             d <- dist_mean_n[j]
#             xbar <- xbar_n
#             
#             delta_start <- c(b^2, b^2)
#             
#             ds <- nleqslv::nleqslv(x = delta_start,
#                                    fn = d_eq,
#                                    jac = d_jb,
#                                    control=list(btol=.01, delta="newton")
#             )
#             ds$x
#             
#           }, b = min(dist_mean_n[-tobe_sampled]), xbar = xbar_n)
#           
#           Tsdata12n <- Tsdata1n
#           Tsdata12n[tobe_sampled, ] <- Tsdata1n[tobe_sampled, ] + t(deltas)
#           plot(Tsdata12n)
#           
#               row.names(Tsdata12n) <- 1:nrow(Tsdata12n)
#               #FRIC and FEVE
#               FDS_comm <- FD::dbFD(Tsdata12n, stand.x = F, calc.FRic = T, stand.FRic = F, 
#                                    calc.CWM = F, calc.FGR = F, calc.FDiv = T)
#               message("done with that")
#               
#               #TOP and modified TOP
#               tops_comm <- TOP(Tsdata12n)
#               message("Done with TOP starting with TED")
#               
#               ## change in community structure
#               #TED
#               teds_comm1 <- TED(Tsdata12n, rescale = T, base_matrix = b_matrix2,
#                                 modification = "None")
#               #modified TED
#               b_matrix1 <- base_traitmatrix_new(Tsdata12n, gen_matrix = gen_matrix, 
#                                                 rescale = T, parallel = F)
#               teds_comm2 <- TED(Tsdata12n, rescale = F, base_matrix = b_matrix1,
#                                 modification = "None")
#               
#         }
#      
#      Comm_Change <- data.frame(FRIC = FDS_comm$FRic, TOP = tops_comm$TA, 
#                                TOPM = tops_comm$TA2, 
#                                FEVE = FDS_comm$FEve, TED = teds_comm1, TEDM = teds_comm2, 
#                                RAO = FDS_comm$RaoQ, FDIS = FDS_comm$FDis)
#      message("Done with TED, next replicate commences")
#      #return(Comm_Change = Comm_Change)
#   #})
#   Results_comm2[[i]] <- Comm_Change #replicates_richness %>% do.call(rbind.data.frame, .)
# }
# for(i in 11:20) {
#     #TOPP <- TEDD <- data.frame()
#     Tsdata1_2and3D <- matrix(NA, ncol = n_dim, nrow = n_row[i])
#     Tsdata1_2and3D <- apply(Tsdata1_2and3D, 2, function(x) runif(n_row[i], 0, 1))
#     middle <- mrfDepth::hdepthmedian(Tsdata1_2and3D, maxdir = 500)$median
#     dist_mean <- apply(Tsdata1_2and3D - middle, 1, function(x){
#       L1_inf(x)
#     })
#     mid_points <- order(dist_mean)[1:(length(dist_mean) / 3)]#which(dist_mean < quantile(dist_mean, 0.20))
#   
#     #T1
#     T1 <- Tsdata1_2and3D[-mid_points, ]
#     T1b <- Tsdata1_2and3D
#     bb <- mid_points#sample(mid_points, length(mid_points) / 4)
#   
#     #T1b
#     T1b[bb, 1]  <- T1b[bb, 1] - min(T1b[bb, 1])# - max(T1b[bb, 1])
#     T1b[bb, 2]  <- T1b[bb, 2] + 2.5
#     
#     #Scenario T_1c
#     T1c <- T1b
#     T1c[bb,2] <- T1c[bb, 2] + 3.5
#   
#     #T2a
#     pp1 <- base::which(Tsdata1_2and3D[, 1] <= 0.3)
#     T2a <- Tsdata1_2and3D[-pp1, ]
#     pp1b <- base::which(T2a[,1] <= 0.35)
#     pp1c <- base::which(T2a[pp1b,2] == max(T2a[pp1b,2]) )
#     T2a[pp1b[pp1c], 1] <- 0.1
#     T2a[pp1b[pp1c], 2] <- max(T2a[,2])
#   
#     #T2b
#     T2b <- T2a
#     pp2 <- base::sample(base::which(T2b[,1] > 0.3),12)
#     T2b[pp2,1] <- 0.0
#     T2b[pp1b, 1] <- 0.0
#   
#     print(i)
#     parallel::clusterExport(ncl, 
#                             c("i", "Tsdata1_2and3D","T1",
#                               "T1b", "T1c", "T2a", "T2b"))
#     
#   sfc <- snow::clusterApply(cl = ncl, x = 1:10, fun = function(j) {
#     #M1
#     b_matrix <- matrix(NA, nrow = nrow(Tsdata1_2and3D), ncol = n_dim)
#     b_matrix <- apply(b_matrix, 2, function(x)runif(nrow(Tsdata1_2and3D), 0, 1))
#     #M2
#     b_matrix_new <- base_traitmatrix_new(Tsdata1_2and3D)
#     
#     TOP_T1aT2 <- TOP(Tsdata1_2and3D); TOP_T1aT2_new <- TOP_new(Tsdata1_2and3D)
#     TOP_T1 <- TOP(T1); TOP_T1_new <- TOP_new(T1)
#     TOP_T1b <- TOP(T1b); TOP_new_T1b <- TOP_new(T1b)
#     TOP_T1c <- TOP(T1c); TOP_new_T1c <- TOP_new(T1c)
#     TOP_T2a <- TOP(T2a); TOP_new_T2a <- TOP_new(T2a)
#     TOP_T2b <- TOP(T2b); TOP_new_T2b <- TOP_new(T2b)
#     
#     TOPP <- data.frame(TOP_T1aT2$TA, TOP_T1aT2_new$TA, TOP_T1aT2_new$N,
#                        TOP_T1$TA, sum(TOP_T1_new$Areas)/TOP_T1aT2_new$N, TOP_T1aT2_new$N,
#                        TOP_T1b$TA, TOP_new_T1b$TA, TOP_new_T1b$N,
#                        TOP_T1c$TA, TOP_new_T1c$TA, TOP_new_T1c$N,
#                        TOP_T2a$TA, TOP_new_T2a$TA, TOP_new_T2a$N,
#                        TOP_T2b$TA, TOP_new_T2b$TA, TOP_new_T2b$N,
#                        TOP_T1aT2$Areas[1], TOP_T1aT2_new$TA/TOP_T1aT2_new$N,
#                        nrow(T1),nrow(T1b), nrow(T1c), nrow(T2a), nrow(T2b))
#     names(TOPP) <- c("T1aT2","T1aT2_new","N_T1aT2","T1","T1_new","N_T1","T1b","T1b_new","N_T1b",
#                      "T1c","T1c_new","N_T1c","T2a","T2a_new","N_T2a","T2b","T2b_new","N_T2b",
#                      "T2c","T2c_new","n_T1","n_T1b","n_T1c","n_T2a","n_T2b")
#     
#     TEDD <- data.frame(TED(Tsdata1_2and3D, base_matrix = b_matrix), 
#                        TED(Tsdata1_2and3D, base_matrix = b_matrix_new),
#                        TED(T1, base_matrix = b_matrix), TED(T1, base_matrix = b_matrix_new),
#                        TED(T1b, base_matrix = b_matrix), TED(T1b, base_matrix = b_matrix_new),
#                        TED(T1c, base_matrix = b_matrix), TED(T1c, base_matrix = b_matrix_new),
#                        TED(T2a, base_matrix = b_matrix), TED(T2a, base_matrix = b_matrix_new),
#                        TED(T2b, base_matrix = b_matrix), TED(T2b, base_matrix = b_matrix_new))
#     
#     TEDD_classic <- data.frame(TED(Tsdata1_2and3D, base_matrix = b_matrix2),
#                        TED(T1, base_matrix = b_matrix2),
#                        TED(T1b, base_matrix = b_matrix2),
#                        TED(T1c, base_matrix = b_matrix2),
#                        TED(T2a, base_matrix = b_matrix2), 
#                        TED(T2b, base_matrix = b_matrix2))
#     names(TEDD) <- c("T1aT2","T1aT2_new","T1","T1_new","T1b","T1b_new","T1c",
#                      "T1c_new","T2a","T2a_new","T2b","T2b_new")
#     names(TEDD_classic) <- c("T1aT2", "T1", "T1b", "T1c", "T2a", "T2b")
#       
#     return(list(TOPS = TOPP, TEDS = TEDD, TEDS_Classic = TEDD_classic))
#   })
#   #names(TOPP) <- c("T1aT2","T1aT2_new","N_T1aT2","T1","T1_new","N_T1","T1b","T1b_new","N_T1b",
#   #                 "T1c","T1c_new","N_T1c","T2a","T2a_new","N_T2a","T2b","T2b_new","N_T2b",
#   #                 "T2c","T2c_new","n_T1","n_T1b","n_T1c","n_T2a","n_T2b")
#   
#   #names(TEDD) <- c("T1aT2","T1aT2_new","T1","T1_new","T1b","T1b_new","T1c",
#   #                 "T1c_new","T2a","T2a_new","T2b","T2b_new",
#   #                 "n_T1","n_T1b","n_T1c","n_T2a","n_T2b") 
#   
#   #TOPS[[i]] <-  TOPP
#   TOPS[[i]] <- lapply(sfc, "[[", 1) %>% do.call(rbind.data.frame, .)
#   #TEDS[[i]] <-  TEDD
#   TEDS[[i]] <- lapply(sfc, "[[", 2) %>% do.call(rbind.data.frame, .)
#   TEDS_Classic[[i]] <- lapply(sfc, "[[", 3) %>% do.call(rbind.data.frame, .)
# }


library(tidyverse)
library(ggpubr)
library(mvtnorm)
library(gridExtra)

x1 <- rep(1:5,each=5)
x2 <- rep(1:5,times=5)
X <- cbind(x1,x2)

#par(mfrow = c(1,2))
plot(X, pch = 19, col = 2, xlab = "Trait 1", ylab = "Trait 2", bty = "l")
points(X[13, 1], X[13, 2], col = "blue", pch = 19)

arrows(x0 = 3, y0 = 3, x1 = 4.5, y1 = 4.5, 
       angle = 45, length = 0.30, lty = 2, lwd = 2.5, col = 6)

X3 <- X[-which((X[,1]>=2 & X[,1]<=4)&(X[,2]>=2 & X[,2]<=4)),]

source("Trait Functions/TOP.R")
source("Trait Functions/TOP_New.R")
source("Trait Functions/basetraitM.R")
source("Trait Functions/basetraitM_new.R")
source("Trait Functions/distances.R")

######### Examples from Simulation scenario one ######################

mu_start <- c(5, 10)
delta_mu <- seq(0, 10, by = 0.5)
rho_start <- 0.007
#equal variance
sig_11_start <- 0.85
cv_start <- sig_11_start*sig_11_start*rho_start

#png(filename = "../../Other Presentations/Simulation_Case1_20.png", 
#    width = 800, height = 500)
#par(mfrow = c(1, 3))
scenario1_plots <- vector("list", 4)
j <- 1
for(i in 1:length(delta_mu)) {
  mu_change <- c(mu_start[1] + delta_mu[i], mu_start[2])
  #parallel::clusterExport(ncl, c("i", "mu_change"))
  print(i)
  
  ddata <- rmvnorm(800, mean = mu_start, sigma = matrix(c(sig_11_start, 
                                                          cv_start,cv_start, 
                                                          sig_11_start), 
                                                        nrow = 2, ncol = 2))
  
  if(i == 1 ) {
    ddat <- as.data.frame(ddata)
    
    p1 <- ddat %>% ggplot(aes(x = V1, y = V2)) + 
      geom_point(color = "red", size = 3) + 
      theme_pubr() + 
      labs(x = "Trait 1", y = "Trait 2") + 
      theme(axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            axis.title = element_text(size = 40), 
            title = element_text(size = 40),
            plot.title = element_text(hjust = 0.5)) + 
      ggtitle(paste0("Location Shift = ", delta_mu[i]))
    scenario1_plots[[i]] <- p1
    
    
  } else if(delta_mu[i] %in% c(3, 9)) {
    
    j <<- j + 1
    ddat <- ddata
    tobe_shifted <- sample(1:nrow(ddata), 0.20*nrow(ddat))
    ddat[tobe_shifted, 1] <- ddat[tobe_shifted, 1] - mu_start[1] + mu_change[1]
    
    ddat <- as.data.frame(ddat)
    
    p2 <- ddat %>% ggplot(aes(x = V1, y = V2)) + 
      geom_point(color = "red", size = 3) + 
      theme_pubr() + 
      labs(x = "Trait 1", y = "Trait 2") + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(), 
            axis.title = element_text(size = 40), 
            title = element_text(size = 40),
            plot.title = element_text(hjust = 0.5)) + 
      ggtitle(paste0("Location Shift = ", delta_mu[i]))
    
    scenario1_plots[[j]] <- p2
    
  }
  
}

### expectations from scenario one
demonstration <- data.frame(shift = 0:10, Value = seq(0, 1, length = 11)
)

p4 <- demonstration %>% ggplot(aes(x = shift, y = Value)) + 
  geom_line(aes(x = seq(0, 3, length.out = 11), 
                y = seq(0.6, 0.8, length.out = 11)), 
            linetype = "dashed", size = 2) + 
  #evennes
  geom_segment(aes(x = 3, y = 0.8, xend = 5.5, yend = 0.15), 
               size = 2, arrow = arrow(angle = 45), 
               linetype = "dashed") +
  geom_text(aes(x = 5.3, y = 0.1), label = "Evenness", size = 9) +
  #richness
  geom_segment(aes(x = 0, y = 0.4, xend = 4.5, yend = 0.8), 
               size = 2, arrow = arrow(angle = 45), 
               linetype = "solid") +
  geom_text(aes(x = 4.28, y = 0.83), label = "Richness", size = 9) +
  #divergence
  geom_segment(aes(x = 0, y = 0.2, xend = 6.5, yend = 0.8), 
               size = 2, arrow = arrow(angle = 45), 
               linetype = "dotted") +
  geom_text(aes(x = 5.95, y = 0.83), label = "Divergence", size = 9) +
  labs(x = "Location Shift", y = "Index value") + 
  theme_pubr() +
  scale_y_continuous(breaks = seq(0.1, 1.0, 0.1)) +
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_text(size = 40), 
        title = element_text(size = 40),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Expected trend")

scenario1_plots[[4]] <- p4

final_plot_scenario1 <- grid.arrange(scenario1_plots[[1]], scenario1_plots[[2]], 
             scenario1_plots[[3]], scenario1_plots[[4]], 
             ncol = 4)
ggsave("Images/final_plots_s1.png", plot = final_plot_scenario1, 
       device = "png", 
       height = 15, width = 34)

######### Examples from Simulation scenario Two ######################

sceneraio_2 <- matrix(NA, ncol = 2, nrow = 1000) 
sceneraio_2 <- apply(sceneraio_2, 2, function(x) runif(800, 0, 1))
delete <- seq(0.10, 0.40, by = 0.10)#c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)

#mean and distance of each vector to the center
xbar <- apply(sceneraio_2, 2, mean)
dist_mean <- apply(sceneraio_2 , 1, function(x) {
  eucldis(x - xbar)
}) 

xbar <- apply(sceneraio_2, 2, mean)
dist_mean <- apply(sceneraio_2 , 1, function(x) {
  eucldis(x - xbar)
}) 

scenario2_plots <- vector("list", 4)
j <- 1
for(i in 1:(length(delete) + 1)) {
  
  if(i == 1) {
    ddat <- as.data.frame(sceneraio_2)
    p1 <- ddat %>% ggplot(aes(x = V1, y = V2)) + 
      geom_point(color = "red", size = 3) + 
      theme_pubr() + 
      labs(x = "Trait 1", y = "Trait 2") + 
      theme(axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            axis.title = element_text(size = 40), 
            title = element_text(size = 40),
            plot.title = element_text(hjust = 0.5)) + 
      ggtitle(paste0("Shifted = ", "None"))
    scenario2_plots[[i]] <- p1
    
  } else if(i %in% c(3, 5) ) {
    
    print(i)
    j <<- j + 1
    tobe_sampled <- which(dist_mean <= quantile(dist_mean, probs = delete[i-1]))
    
    deltas <- sapply(tobe_sampled, dts <- function(j, b, xbar) {
      
      environment(d_eq) <- environment()
      environment(d_jb) <- environment()
      
      x <- sceneraio_2[j, ]
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
    
    sceneraio_2b <- sceneraio_2
    sceneraio_2b[tobe_sampled, ] <- sceneraio_2[tobe_sampled, ] + t(deltas)
    
    row.names(sceneraio_2b) <- 1:nrow(sceneraio_2b)
    ddat <- as.data.frame(sceneraio_2b)
    
    p2 <- ddat %>% ggplot(aes(x = V1, y = V2)) + 
      geom_point(color = "red", size = 3) + 
      theme_pubr() + 
      labs(x = "Trait 1", y = "Trait 2") + 
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(), 
            axis.title = element_text(size = 40), 
            title = element_text(size = 40),
            plot.title = element_text(hjust = 0.5)) + 
      ggtitle(paste0("Shifted = ", 100*delete[i-1], "%"))
      scenario2_plots[[j]] <- p2
    
  } else next
  
}

#expectation for scenario 2
p42 <- demonstration %>% ggplot(aes(x = shift, y = Value)) + 
    #evenness
  geom_segment(aes(x = 0, y = 0.8, xend = 5.5, yend = 0.15), 
               size = 2, arrow = arrow(angle = 45), 
               linetype = "dashed") +
  geom_text(aes(x = 5.6, y = 0.1), label = "Evenness", size = 9) +
  #richness
  geom_segment(aes(x = 0, y = 0.7, xend = 3.95, yend = 0.15), 
               size = 2, arrow = arrow(angle = 45), 
               linetype = "solid") +
  geom_text(aes(x = 4.05, y = 0.1), label = "Richness", size = 9) +
  #divergence
  geom_segment(aes(x = 0, y = 0.1, xend = 6.5, yend = 0.8), 
               size = 2, arrow = arrow(angle = 45), 
               linetype = "dotted") +
  geom_text(aes(x = 5.95, y = 0.83), label = "Divergence", size = 9) +
  labs(x = "% Shifted", y = "Index value") + 
  theme_pubr() +
  scale_y_continuous(breaks = seq(0.1, 1.0, 0.1)) +
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_text(size = 40), 
        title = element_text(size = 40),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Expected trend")

scenario2_plots[[4]] <- p42

final_plot_scenario2 <- grid.arrange(scenario2_plots[[1]], scenario2_plots[[2]], 
                                     scenario2_plots[[3]], scenario2_plots[[4]], 
                                     ncol = 4)
ggsave("Images/final_plots_s2.png", plot = final_plot_scenario2, 
       device = "png", 
       height = 15, width = 34)

######### Examples from Simulation scenario Three ######################

scenario3_plots <- vector("list", 4)
j <- 0
for(i in c(3, 4, 5)) {
  sceneraio_2 <- rmvnorm(n_row[i], mean = c(5, 5), sigma = matrix(c(sig_11_start, 
                                                                cv_start,cv_start, 
                                                                sig_11_start), 
                                                              nrow = 2, ncol = 2))
  j <<- j + 1
  ddat <- as.data.frame(sceneraio_2)
  
  p2 <- ddat %>% ggplot(aes(x = V1, y = V2)) + 
    geom_point(color = "red", size = 3) + 
    theme_pubr() + 
    labs(x = "Trait 1", y = "Trait 2") + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(), 
          axis.title = element_text(size = 40), 
          title = element_text(size = 40),
          plot.title = element_text(hjust = 0.5)) + 
    ggtitle(paste0("N = ", n_row[i]))
  scenario3_plots[[j]] <- p2
}
#expectation for scenario 3
p43 <- demonstration %>% ggplot(aes(x = shift, y = Value)) + 
  #evenness
  geom_segment(aes(x = 0, y = 0.2, xend = 5.5, yend = 0.2), 
               size = 2, arrow = arrow(angle = 45), 
               linetype = "dashed") +
  geom_text(aes(x = 3.5, y = 0.21), label = "Evenness", size = 9) +
  #richness
  geom_segment(aes(x = 0, y = 0.3, xend = 5.5, yend = 0.3), 
               size = 2, arrow = arrow(angle = 45), 
               linetype = "solid") +
  geom_text(aes(x = 3.5, y = 0.31), label = "Richness", size = 9) +
  #divergence
  geom_segment(aes(x = 0, y = 0.4, xend = 5.5, yend = 0.4), 
               size = 2, arrow = arrow(angle = 45), 
               linetype = "dotted") +
  geom_text(aes(x = 3.5, y = 0.41), label = "Divergence", size = 9) +
  labs(x = "Number", y = "Index value") + 
  theme_pubr() +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.1)) +
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_text(size = 40), 
        title = element_text(size = 40),
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Expected trend")

scenario3_plots[[4]] <- p43

final_plot_scenario3 <- grid.arrange(scenario3_plots[[1]], scenario3_plots[[2]], 
                                     scenario3_plots[[3]], scenario3_plots[[4]], 
                                     ncol = 4)
ggsave("Images/final_plots_s3.png", plot = final_plot_scenario3, 
       device = "png", 
       height = 15, width = 34)











################## every other thing begins here
TOP(X)
TOP_new(cbind(x1,x2))



###shifting a point until one of the convex hull
a <-  X[(X[, 1]==3 & X[, 2]==3), 2] #<- 3.5
b <- X[(X[, 1]==3 & X[, 2]==3), 1] #<- 2.5
rest2 <- rest <- c(); position <- c()
i <- 0
X2 <- X
while(a >= 1 & b >= 1) {
  X2[(X[, 1]==3 & X[, 2]==3), 2] <- a
  X2[(X[, 1]==3 & X[, 2]==3), 1] <- b
  points(c(a, b), col = 4)
  tt <- TOP(X2)
  rest <- c(rest, tt$TA); rest2 <- c(rest2, tt$TA2)
  position <- c(position, eucldis(c(3,3) - c(b, a)))
  b <- b - 0.05
  a <- a - 0.05
  i <- i + 1
  print(i)
  print(c(a, b))
}

par(mfrow=c(1, 2))
plot(X, pch = 19, col = 2, xlab = "Trait 1", ylab = "Trait 2")
points(X[13, 1], X[13, 2], col = "blue", pch = 19)
plot(position,rest,type="p",ylab="TOP",xlab="Distance to c(3,3)", main="Classical TOP")
points(position, rest, pch = 19, col = 2)


plot(position,rest2,type="l",ylab="TOP=sum(Area*no of points on hull)",xlab="Distance to c(3,3)",
     main="Thomas Suggestion")

###TED
TED(X)



############Some simulation concept
TOP(rbind(X,c(1,5)))




library(mvtnorm)
eucldis <- function(x) sqrt(sum(x^2))

#L1 distance
L1_inf <- function(x) max(abs(x))

test_data <- matrix(NA, nrow = 100, ncol = 2)
test_data <- apply(test_data, 2, function(x) runif(100, 0, 1))



#Rescaling
rescale <- function(x){
  apply(x,MARGIN = 2,function(y){
    mn <- min(y);mx <- max(y);n <- length(y)
    a <- (n/(n-1)) * (mn - mx/n)
    b <- (n/(n-1)) * (mx - mn/n)
    (y - a)/(b - a)
  })
}
#scenario 2
sceneraio_2 <- rescale(rmvnorm(800, mean = c(5, 5), sigma = matrix(c(0.85, 0.005,0.005, 0.85), 
                                                                   nrow = 2, ncol = 2))
)
conv_hull <- geometry::convhulln(sceneraio_2, option = "Fx")
mconv_hull <- as.data.frame(unique(sceneraio_2[conv_hull, ]))


plot(sceneraio_2, xlab = "Trait 1", ylab = "Trait 2", pch = 19, col = 2, type = "p", 
     xaxt = "n", yaxt = "n", bty = "l")
axis(side = 1, labels = F, col.ticks = NA)
axis(side = 2, labels = F, col.ticks = NA)


begin <- sample(1:nrow(mconv_hull), 1)
bbegin <- begin
vvec <- 1:nrow(mconv_hull)

function(mconv_hull, begin, vvec, bbegin) {
  #calculate distance between the random start on the convex-hull and the 
  #other points on the convex-hull
  
  while(length(vvec) > 1) {
    dds <- unlist(sapply(vvec, function(i) {
      y <- eucldis(mconv_hull[i, ] - mconv_hull[begin, ])
      return(y)
    }))
    
    #among the non-zeros, obtain the minimum
    ending <- vvec[which(dds == min(dds[dds != 0 ]))]
    
    #draw the line between the beginning and ending
    segments(x0 = mconv_hull$V1[vvec[vvec == begin]], y0 = mconv_hull$V2[vvec[vvec == begin]], 
             x1 = mconv_hull$V1[vvec[vvec == ending]], y1 = mconv_hull$V2[vvec[vvec == ending]], 
             col = "green4", lwd = 3, lty = 1)
    
    #reduce the convex_hull by begin
    vvec <- vvec[-which(vvec == begin)]
    begin <- ending 
    if(length(vvec) == 1) break
  }
  
  if(length(vvec) == 1) {
    
    segments(x0 = mconv_hull$V1[bbegin], y0 = mconv_hull$V2[bbegin], 
             x1 = mconv_hull$V1[ending], y1 = mconv_hull$V2[ending], 
             col = "green4", lwd = 3, lty = 1)
  }
  #legend("topright", legend = "FRIC = Polygon Volume")
  #mtext(text = "FRIC = ", line = 3, col = "black", at = 0.7)
  #mtext(text = "Polygon", line = 3, col = "green", at = 0.85)
  mtext(text = expression("FRIC=" * phantom("polygon ") * "volume"), line = 3, col = "black", at = 0.8)
  mtext(text = expression(phantom("FRIC=") * "polygon " * phantom("volume")), line = 3, col = "green4", at = 0.8)
}

###plotting minimum spanning tree
library(ape)
test_data <- matrix(NA, nrow = 20, ncol = 2)
test_data <- apply(test_data, 2, function(x) runif(20, 0, 1))
plot(mst(dist(test_data)), graph = "nsca")



######## TED Distribution Establishing
test_data_dist <- as.matrix(stats::dist(test_data, diag = F, upper = T))
test_data_dist <- as.numeric(test_data_dist[base::upper.tri(test_data_dist)])

##fitting normal distribution to the data
#test_dist_fit <- fitdist(test_data_dist, distr = "norm", method = "mle")
dbar <- mean(test_data_dist)
vdbar <- (ncol(test_data) * (2 * nrow(test_data) + 3)) / (90 * nrow(test_data) * (nrow(test_data) - 1))
vd <- nrow(test_data) * vdbar

#hist(test_data_dist, probability = T)
#dd <- function(x) dnorm(x, mean = test_dist_fit$estimate[1], sd = test_dist_fit$estimate[2])
#curve(dd, from = min(test_data_dist), to = max(test_data_dist), n = 4950)

#lines(x = test_data_dist, y = 
##        dnorm(test_data_dist, mean = test_dist_fit$estimate[1], sd = test_dist_fit$estimate[2]), 
#      col = 2)
#using ggplot2

ggplot(data = NULL, aes(x = test_data_dist)) + 
  geom_histogram(aes(y = ..density..), stat = "bin", bins = 17) + 
  geom_line(aes(x = test_data_dist, 
                y = dnorm(test_data_dist, mean = dbar, 
                          sd = sqrt(vd)))) + 
  theme_minimal()

test_data_density <- density(test_data_dist, n = 1024, from = min(test_data_dist), 
                             to = max(test_data_dist))

#approach one (using the theoretical mean and variance)
test_data_dnorm <- dnorm(test_data_density$x, mean = dbar, 
                         sd = sqrt(vd))

#plotting both density distributions
data.frame(Distances = test_data_density$x, Observed = test_data_density$y, 
           Theoretical = test_data_dnorm) %>% 
  tidyr::gather(key = "Distribution", value = "Density", -Distances) %>% 
  ggplot(aes(x = Distances, y = Density, group = Distribution, color = Distribution)) + 
  geom_line() + theme_minimal()

flexmix::KLdiv(cbind(test_data_density$y, test_data_dnorm))

#approach two(obtain the sampling distribution of d-bar)
dbar_samp <- unlist(map(1:4950, function(x) {
  ss <- sample(test_data_dist, size = length(test_data_dist), replace = T)
  return(mean(ss))
}))

flexmix::KLdiv(cbind(test_data_density$y, density(dbar_samp)$y))


###simulating from the mvnorm

png(filename = "../../Other Presentations/IntraspecificDiversity1.png", 
    width = 1800, height = 1500)
plot(x1, type = "n", xaxt = "n", yaxt = "n", bty = "l", xlab = "", ylab = ""
     , axes = T)
axis(side = 1, labels = F, col.ticks = NA, col = NA)
axis(side = 2, labels = F, col.ticks = NA, col = NA)
mtext("Trait 1", side = 1, line = 3, cex = 4.5)
mtext("Trait 2", side = 2, line = 1, cex = 4.5)

dev.off()

#scenario 1
sceneraio_1 <- matrix(NA, nrow = 800, ncol = 2)
sceneraio_1 <- apply(sceneraio_1, 2, function(x) runif(5, 0, 1))

set.seed(1992)
sceneraio_2 <- rescale(rmvnorm(20, mean = c(5, 5), sigma = matrix(c(0.85, 0.005,0.005, 0.85), 
 nrow = 2, ncol = 2))
)

############### New plots
png(filename = "Figures/Trait_Diversity_Indices_all.png", 
    width = 1800, height = 1800)
par(mfrow = c(3, 2))
# FRic
covhull <- chull(sceneraio_2)
plot(sceneraio_2, pch = 20, col = 1, cex = 10, type = "p", 
     bty = "l", xlab = "", ylab = "", ylim = c(0, 1.1),
     xaxt = "n", yaxt = "n"
    )
polygon(sceneraio_2[covhull, ], lwd = 4, border = "black")
points(sceneraio_2, pch = 20, col = 1, cex = 7.5)
axis(side = 1, labels = F, col.ticks = NA, col = NA)
axis(side = 2, labels = F, col.ticks = NA, col = NA)
mtext("Trait 1", side = 1, line = 3, cex = 2.5)
mtext("Trait 2", side = 2, line = 1, cex = 2.5)
text(x = min(sceneraio_2[,1]), y = 1.1, "a", family = "sans", 
     font = 2, cex = 4)

# TOP
plot(sceneraio_2, pch = 20, col = 1, cex = 10, type = "p", 
     bty = "l", xlab = "", ylab = "", ylim = c(0, 1.1),
     xaxt = "n", yaxt = "n"
)
axis(side = 1, labels = F, col.ticks = NA, col = NA)
axis(side = 2, labels = F, col.ticks = NA, col = NA)
mtext("Trait 1", side = 1, line = 3, cex = 2.5)
mtext("Trait 2", side = 2, line = 1, cex = 2.5)
i <- 1
sceneraio_2b <- sceneraio_2
while(TRUE) {
  i <- i+1
  covhull <- chull(sceneraio_2b)
  polygon(sceneraio_2b[covhull, ], lwd = 4, border = "black")#
  points(sceneraio_2b[covhull, ], pch = 20, col = 1, cex = 10)
  sceneraio_2b <- sceneraio_2b[-covhull, ]
  if(class(sceneraio_2b) != "matrix") {
    break
  } else if(class(sceneraio_2b) == "matrix" & nrow(sceneraio_2b) <= ncol(sceneraio_2b)) {
    break
  } else next
  
}
text(x = min(sceneraio_2[,1]), y = 1.1, "b", family = "sans", 
     font = 2, cex = 4)

# FEve
plot(ComputeMST(sceneraio_2), xaxt = "n", yaxt = "n", ylab = "", xlab = "", 
     pch = 20, ylim = c(0, 1.1),
     bty = "l", col.pts = 1, cex = 10, 
     col.segts = 1, lty = 5, lwd = 5)
axis(side = 1, labels = F, col.ticks = NA, col = NA)
axis(side = 2, labels = F, col.ticks = NA, col = NA)
mtext("Trait 1", side = 1, line = 3, cex = 2.5)
mtext("Trait 2", side = 2, line = 1, cex = 2.5)
text(x = min(sceneraio_2[,1]), y = 1.1, "c", family = "sans", 
     font = 2, cex = 3)

# TED
obs_density <- as.matrix(stats::dist(sceneraio_2, diag = F, upper = T))
obs_density <- as.numeric(obs_density[base::upper.tri(obs_density)])
obs_density <- density(obs_density, from = min(obs_density), to = max(obs_density))

x_sim1 <- base_traitmatrix("cube", ncols = 2, dsize = 200)
sim1_density <- as.matrix(stats::dist(x_sim1, diag = F, upper = T))
sim1_density <- as.numeric(sim1_density[base::upper.tri(sim1_density)])
sim1_density <- density(sim1_density, from = min(sim1_density), to = max(sim1_density))

#TEDM
gen_matrix <- matrix(NA, nrow = 500, ncol = 2)
gen_matrix <- apply(gen_matrix, 2, function(x) runif(nrow(gen_matrix), 0, 1))
x_sim2 <- base_traitmatrix_new(x = sceneraio_2, gen_matrix = gen_matrix, rescale = FALSE)
sim2_density <- as.matrix(stats::dist(x_sim2, diag = F, upper = T))
sim2_density <- as.numeric(sim2_density[base::upper.tri(sim2_density)])
sim2_density <- density(sim2_density, from = min(sim2_density), 
                        to = max(sim2_density))

plot(x = obs_density$x, y = obs_density$y, main = "",
     xaxt = "n", yaxt = "n", ylab = "", xlab = "", 
     type = "l", ylim = c(0.00, 2.5), xlim = c(-0.01, 1.0),
     bty = "l", col.pts = 1, 
     col.segts = 1, lty = 1, lwd = 2)
axis(side = 1, labels = F, col.ticks = NA, col = NA)
axis(side = 2, labels = F, col.ticks = NA, col = NA)
mtext("Distance", side = 1, line = 3, cex = 2.5)
mtext("Density", side = 2, line = 1, cex = 2.5)
lines(x = sim1_density$x, y = sim1_density$y, lty = 1, lwd = 2, 
      col = 2)
lines(x = sim2_density$x, y = sim2_density$y, lty = 1, lwd = 2,
      col = 3)
legend(x = 0.1, y = 2.6, legend = c("observed", "Discrete-Uniform (TED)", 
                                    "Continuous-Uniform (TEDM)"), 
       col = c(1, 2, 3), lty = c(1, 1, 1), horiz = TRUE, bty = "n", 
       cex = 2.5, ncol = 1,
       text.width = c(0.1, 0.1, 0.19),
       x.intersp = 0.16)
text(x = 0.01, y = 2.5, "d", family = "sans", 
     font = 2, cex = 3)

# FDis and Rao
plot(sceneraio_2, pch = 19, 
     col = 1, cex = 7.5, bg = "blue", 
     xaxt = "n", yaxt = "n", 
     ylab = "", xlab = "",
     ylim = c(0, 1.1),
     bty = "l")
midpoint <- apply(sceneraio_2, 2, mean)
#points(x = midpoint[1], y = midpoint[2], pch = 18, col = "red3", cex = 7.5)
segments(x0 = midpoint[1], y0 = midpoint[2], 
         x1 = sceneraio_2[,1], sceneraio_2[,2], lty = 3)
axis(side = 1, labels = F, col.ticks = NA, col = NA)
axis(side = 2, labels = F, col.ticks = NA, col = NA)
text(0.32, 0.87, labels = expression(d[j]), cex = 3.9, 
     font = 2, family = "sans", col = 1)
text(x = midpoint[1], y = midpoint[2], "c", font = 2, 
     family = "sans", col = 1, cex = 3.9)
mtext("Trait 1", side = 1, line = 3, cex = 2.5)
mtext("Trait 2", side = 2, line = 1, cex = 2.5)
text(x = min(sceneraio_2[,1]), y = 1.1, "e", family = "sans", 
     font = 2, cex = 3)
dev.off()

######################################## old plots
#FEve
#### plotting minimum spanning tree (FEve)
library(emstreeR)
png(filename = "../../Other Presentations/FEve_Demonstration.png", 
    width = 1800, height = 1500)
plot(ComputeMST(sceneraio_2), xaxt = "n", yaxt = "n", ylab = "", xlab = "", 
     axes = TRUE, pch = 20, ylim = c(0, 1),
     bty = "l", col.pts = 1, cex = 7.5, 
     col.segts = 1, lty = 5, lwd = 3)
axis(side = 1, labels = F, col.ticks = NA, col = NA)
axis(side = 2, labels = F, col.ticks = NA, col = NA)
mtext("Trait 1", side = 1, line = 3, cex = 4.5)
mtext("Trait 2", side = 2, line = 1, cex = 4.5)
dev.off()


#TED
png(filename = "../../Other Presentations/TED_Demonstration.png", 
    width = 1800, height = 1500)
par(mfrow = c(1, 2))

plot(sceneraio_2, pch = 20, col = 1, cex = 7.5, type = "p"
     , bty = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
axis(side = 1, labels = F, col.ticks = NA)
axis(side = 2, labels = F, col.ticks = NA)
mtext("Trait 1", side = 1, line = 3, cex = 5)
mtext("Trait 2", side = 2, line = 1, cex = 5)
mtext("Observed", side = 3, line = 1, cex = 5)

plot(x_sim, pch = 20, col = 1, cex = 7.5, type = "p"
     , bty = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
axis(side = 1, labels = F, col.ticks = NA)
axis(side = 2, labels = F, col.ticks = NA)
mtext("Trait 1", side = 1, line = 3, cex = 5)
mtext("Trait 2", side = 2, line = 1, cex = 5)
mtext("Hypothetical", side = 3, line = 1, cex = 5)
dev.off()

### demonstrating Trait Divergence
png(filename = "../../Other Presentations/TraitDivergenceDemonstration.png", 
    width = 1800, height = 1500)
ddata_div <- data.frame(x = c(3.5, 4.5, 2.5, 0.5, 6, 7), 
                        y = c(1.5, 6.5, 4, 7.5, 1, 3.75)
)
plot(sceneraio_2, pch = 19, 
     col = 1, cex = 7.5, bg = "blue", 
     xlab = "", ylab = "", axes = TRUE, 
     bty = "l")
midpoint <- apply(sceneraio_2, 2, mean)
points(x = midpoint[1], y = midpoint[2], pch = 18, col = "red3", cex = 7.5)

segments(x0 = midpoint[1], y0 = midpoint[2], 
         x1 = sceneraio_2[,1], sceneraio_2[,2], lty = 3)
dev.off()
#mtext("Trait 1", side = 1, line = 2, cex = 2.1)
#mtext("Trait 2", side = 2, line = 1, cex = 2.1)
#text(2, 6.5, labels = expression(d[i]), cex = 3, font = 2, family = "sans", col = "blue")


# FRIC
covhull <- chull(sceneraio_2)
png(filename = "../../Other Presentations/FRIC_Demonstration.png", 
    width = 1800, height = 1500)
plot(sceneraio_2, pch = 20, col = 1, cex = 7.5, type = "p", 
     axes = FALSE, bty = "l", xlab = "", ylab = "")
polygon(sceneraio_2[covhull, ], lwd = 4, border = "red")
dev.off()

#TOP
png(filename = "../../Other Presentations/TOP_Demonstration.png", 
    width = 1800, height = 1500)
i <- 1
plot(sceneraio_2, pch = 20, col = 1, cex = 7.5, type = "p", 
     axes = FALSE, bty = "l", xlab = " ", ylab = "")
#axis(side = 1, labels = F, col.ticks = NA, xpd = T, gap.axis = 0)
#axis(side = 2, labels = F, col.ticks = NA, xpd = T, gap.axis = 0)

#TOP
while(TRUE) {
    i <- i+1
    covhull <- chull(sceneraio_2)
    polygon(sceneraio_2[covhull, ], lwd = 4, border = i)#
    sceneraio_2 <- sceneraio_2[-covhull, ]
    if(class(sceneraio_2) != "matrix") {
      break
    } else if(class(sceneraio_2) == "matrix" & nrow(sceneraio_2) <= ncol(sceneraio_2)) {
      break
    } else next
       
}
dev.off()

mtext(text = expression("TOP=" * phantom("Area1") * "+" * phantom("Area2") * "+" * phantom("Area3")), 
      line = 2, col = 1, at = 0.75)
mtext(text = expression(phantom("TOP=") * "Area1" * phantom("+") * phantom("Area2") * phantom("+") * phantom("Area3")), 
      line = 2, col = 2, at = 0.75)
mtext(text = expression(phantom("TOP=") * phantom("Area1") * phantom("+") * "Area2" * phantom("+") * phantom("Area3")), 
      line = 2, col = 3, at = 0.75)
mtext(text = expression(phantom("TOP=") * phantom("Area1") * phantom("+") * phantom("Area2") * phantom("+") * "Area3"), 
      line = 2, col = 4, at = 0.75)
#mtext(text = expression(phantom("FRIC=") * "polygon " * phantom("volume")), line = 3, col = "green4", at = 0.8)
#legend("topright", , col = 1:i)



points(x = 4 - 0.5, y = 4 - 2.5, pch = 19, col = 1, cex = 3.1)
points(x = 4 + 0.5, y = 4 + 2.5, pch = 19, col = 1, cex = 3.1)
points(x = 4 - 1.5, y = 4 - 0, pch = 19, col = 1, cex = 3.1)
points(x = 4 - 3.5, y = 4 + 3.5, pch = 19, col = 1, cex = 3.1)
points(x = 4 + 2, y = 4 - 3, pch = 19, col = 1, cex = 3.1)
points(x = 4 + 3, y = 4 - 0.25, pch = 19, col = 1, cex = 3.1)



sceneraio_2 <- rescale(rmvnorm(400, mean = c(5, 5), sigma = matrix(c(0.85, 0.005,0.005, 0.85), 
                                                                  nrow = 2, ncol = 2)) )
sceneraio_2b <- sceneraio_2

#sceneraio_2b[1:100,] <- sceneraio_2b[1:100, ] - c(5, 5) + c(5.5, 5.5)


#scenario 3
sceneraio_3a <- rmvnorm(400, mean = c(5, 10), sigma = matrix(c(0.85, 0.005,0.005, 0.85), 
                                                           nrow = 2, ncol = 2))

sceneraio_3ab <- rmvnorm(400, mean = c(8, 10), sigma = matrix(c(0.85, 0.005,0.005, 0.85), 
                                                             nrow = 2, ncol = 2))
sceneraio_3 <- rescale(rbind(sceneraio_3a, sceneraio_3ab))

#scenario 4
sceneraio_4a <- rmvnorm(400, mean = c(5, 10), sigma = matrix(c(0.85, 0.005,0.005, 0.85), 
                                                             nrow = 2, ncol = 2))
sceneraio_4ab <- rmvnorm(400, mean = c(10, 10), sigma = matrix(c(0.85, 0.005,0.005, 0.85), 
                                                              nrow = 2, ncol = 2))
sceneraio_4 <- rescale(rbind(sceneraio_4a, sceneraio_4ab))

#scenario 5
sceneraio_5a <- rmvnorm(400, mean = c(5, 10), sigma = matrix(c(0.85, 0.005,0.005, 0.85), 
                                                             nrow = 2, ncol = 2))
sceneraio_5ab <- rmvnorm(400, mean = c(16, 10), sigma = matrix(c(0.85, 0.005,0.005, 0.85), 
                                                               nrow = 2, ncol = 2))
sceneraio_5 <- rescale(rbind(sceneraio_5a, sceneraio_5ab))





par(mfrow = c(1, 3))
mms = c(0, 3, 9)
map2(.x = 1:3, .y = list(sceneraio_2, 
                        sceneraio_3, sceneraio_5), function(.x, .y) {
  plot(rescale(.y), pch = 19, col = 2, main = "", 
       xlab = " ", ylab = " ", xaxt = "n", yaxt = "n", bty = "l", cex = 2.4
       )
  mtext(paste0("Location Shift = ", mms[.x]), side = 3, line = 1, cex = 1.5)
  mtext("Trait 1", side = 1, line = 3, cex = 3)
  mtext("Trait 2", side = 2, line = 1, cex = 3)
  axis(side = 1, labels = F, col.ticks = NA)
  axis(side = 2, labels = F, col.ticks = NA)
                          
})

par(mfrow = c(1, 2))
plot(sceneraio_1, frame.plot = F, pch = 19, col = 2, 
     main = paste("TOP =", round(TOP(sceneraio_1)$TA, 2), ",", "n =", nrow(sceneraio_1)), 
     xlab = "Trait 1", ylab = "Trait 2")

plot(sceneraio_2, frame.plot = F, pch = 19, col = 2, 
     main = paste("TOP =", round(TOP(sceneraio_2)$TA, 2), ",", "n =", nrow(sceneraio_2)), 
     xlab = "Trait 1", ylab = "Trait 2")

#### Drawing Expectations for each scenario
demonstration <- data.frame(shift = 0:10, Value = seq(0, 1, length = 11)
           )

#demonstration %>% ggplot(aes(x = shift, y = Value, group = 1)) + 
#  geom_line(linetype = "solid", size = 4) +
#  scale_x_continuous(breaks = seq(0, 10, by = 1))  + 
#  scale_y_continuous(breaks = seq(0, 1, by = 0.1), 
#                     labels = seq(0, 1, by = 0.1), 
#                     limits = c(0, 1), 
#                     expand = c(0, 0)) +
#  theme(axis.text.x = element_text(angle = 90, hjust = 0.9, vjust = 0.5, size = 20, color = "black"), 
#        axis.text.y = element_text(color = "black", size = 20),
#        axis.title = element_text(color = "black", size = 18, face = "bold")
#  ) + theme_pubr() + labs(x = "Location Shift", y = "Index Value")

png("../../Other Presentations/Simulation_Case2_Expectation_20.png", 
    width = 500, height = 400)

plot(demonstration, type = "n", xaxt = "n", yaxt = "n", bty = "l", xlab = "", ylab = ""
)
axis(side = 1, labels = F, col.ticks = NA)
axis(side = 2, labels = F, col.ticks = NA)
mtext("% shifted away from center", side = 1, line = 1, cex = 1.9)
mtext("Index Value", side = 2, line = 1, cex = 1.9)
#mtext("Expectations", side = 3, line = 1, cex = 3)

#scenario 1 expectations
#arrows(0, 0, 9, 0.9, length = 0.55, angle = 45, col = 1, lty = 1, lwd = 5)
#arrows(0, 0.9, 9, 0.0, length = 0.55, angle = 45, col = 1, lty = 2, lwd = 5)
#arrows(0, 0.2, 7, 0.9, length = 0.55, angle = 45, col = 1, lty = 3, lwd = 5)

#scenario 3 expectations
#arrows(0, 0.5, 9, 0.5, length = 0.55, angle = 45, col = 1, lty = 1, lwd = 5)
#arrows(0, 0.7, 9, 0.7, length = 0.55, angle = 45, col = 1, lty = 2, lwd = 5)
#arrows(0, 0.3, 9, 0.3, length = 0.55, angle = 45, col = 1, lty = 3, lwd = 5)

#scenation 2 expectations
arrows(x0 = 0, y0 = 0.8, x1 = 7, y1 = 0.10, length = 0.45, angle = 45, col = 1, lty = 1, lwd = 4)
arrows(x0 = 0, y0 = 0.3, x1 = 9, y1 = 0.95, length = 0.48, angle = 45, col = 1, lty = 3, lwd = 4)
arrows(x0 = 0, y0 = 1, x1 = 9, y1 = 0.1, length = 0.45, angle = 45, col = 1, lty = 2, lwd = 4)

dev.off()
#text(x = 9.2, y = 10, labels = "Richness \n Divergence", col = 3, face = "bold")

####Community2
sceneraio_2 <- matrix(NA, ncol = 2, nrow = 1000) #rmvnorm(1000, mean = c(5, 5), sigma = matrix(c(sig_11_start, 
              #                                     cv_start,cv_start, 
              #                                     sig_11_start), 
              #                                   nrow = 2, ncol = 2)) 

sceneraio_2 <- apply(sceneraio_2, 2, function(x) runif(800, 0, 1))
delete <- seq(0.10, 0.40, by = 0.10)#c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)

#mean and distance of each vector to the center
xbar <- apply(sceneraio_2, 2, mean)
dist_mean <- apply(sceneraio_2 , 1, function(x) {
  eucldis(x - xbar)
}) 

png("../../Other Presentations/Simulation_Case2_20.png", 
    width = 500, height = 400)
par(mfrow = c(1, 3))
#sceneraio_2b <- sceneraio_2
for(i in 1:(length(delete) + 1)) {
  
  if(i == 1) {
    plot(sceneraio_2, pch = 19, col = 2, main = "", 
                  xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "l", 
             cex = 2.4)
    axis(side = 1, labels = F, col.ticks = NA)
    axis(side = 2, labels = F, col.ticks = NA)
    mtext("Trait 1", side = 1, line = 2, cex = 2.1)
    mtext("Trait 2", side = 2, line = 1, cex = 2.1)
    mtext(paste0("Shifted = ", "None"), side = 3, line = 1, cex = 1.9)
    #plot(density(as.numeric(dist(sceneraio_2))), main = paste0("Delete = ", "None"))
  } else if(i %in% c(3, 5) ) {
      print(i)
      #middle <- mrfDepth::hdepthmedian(sceneraio_2, maxdir = 500)$median
      tobe_sampled <- which(dist_mean <= quantile(dist_mean, probs = delete[i-1]))
      
      deltas <- sapply(tobe_sampled, dts <- function(j, b, xbar) {

        environment(d_eq) <- environment()
        environment(d_jb) <- environment()

        x <- sceneraio_2[j, ]
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

      #signs <- apply(sceneraio_2[tobe_sampled,] , 1, function(x)eucldis2(x - xbar))
      #changes <- apply(sceneraio_2[tobe_sampled, ], 1, function(x) x - xbar)
      #
      #deltas <- sapply(signs, ass, shift = 0.10, dims = 2)
      
      sceneraio_2b <- sceneraio_2
      sceneraio_2b[tobe_sampled, ] <- sceneraio_2[tobe_sampled, ] + t(deltas)
    
    row.names(sceneraio_2b) <- 1:nrow(sceneraio_2b)
    plot(sceneraio_2b, pch = 19, col = 2, main = "", 
         xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
         bty = "l", 
         cex = 2.4)
    axis(side = 1, labels = F, col.ticks = NA)
    axis(side = 2, labels = F, col.ticks = NA)
    mtext("Trait 1", side = 1, line = 2, cex = 2.1)
    mtext("Trait 2", side = 2, line = 1, cex = 2.1)
    mtext( paste0("Shifted = ", 100*delete[i-1], "%"), side = 3, line = 1, cex = 1.9)
    #yyy <- density(as.numeric(dist(remaining)))
    #lines(yyy$x, yyy$y, main = paste0("Delete = ", 100*delete[i-1], "%"), col = i)
    
  } else next
  
}
dev.off()

##### Simulation 3

n_row <- c(50, 100, 200, seq(500, 3000, by = 500))
#par(mfrow = c(1, 3))
for(i in c(3, 4, 5)) {
  sceneraio_2 <- rmvnorm(1000, mean = c(5, 5), sigma = matrix(c(sig_11_start, 
                                                                cv_start,cv_start, 
                                                                sig_11_start), 
                                                              nrow = 2, ncol = 2)) 
  plot(sceneraio_2, pch = 19, col = 2, main = "", 
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "l", 
       cex = 2.2)
  axis(side = 1, labels = F, col.ticks = NA)
  axis(side = 2, labels = F, col.ticks = NA)
  mtext("Trait 1", side = 1, line = 2, cex = 2.1)
  mtext("Trait 2", side = 2, line = 1, cex = 2.1)
  mtext( paste0("N = ", n_row[i]), side = 3, line = 1, cex = 2.1)
}


# middle <- mrfDepth::hdepthmedian(sceneraio_1, maxdir = 500)$median
# d_2_mid <- apply(sceneraio_1 - middle, 1, function(x) {
#   L1_inf(x)
# })
# 
# move <- order(d_2_mid, decreasing = T)
# rel_column <- sapply(1:100, function(i) {
#   which(sceneraio_1[move[i], ] - middle == min(sceneraio_1[move[i], ] - middle))
# })
# sceneraio_1b <- sceneraio_1
# for(i in 1:100) {
#   sceneraio_1b[move[i], rel_column[i]] <- sceneraio_1b[move[i], rel_column[i]] + d_2_mid[move[100]]
# }

plot(sceneraio_2)



#### TED visualization

s1 <- base_traitmatrix("cube", ncols = 2, dsize = 200)
s12 <- matrix(NA, nrow = 200, ncol = 2)
s12 <- apply(s12, 2, function(x)runif(200, 0, 1))

par(mfrow = c(1, 2))
plot(s1, ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "l", pch = 19, col = 1, 
     main = "TED")
axis(side = 1, labels = F, col.ticks = NA)
axis(side = 2, labels = F, col.ticks = NA)
mtext("Trait 1", side = 1, line = 1)
mtext("Trait 2", side = 2, line = 1)

plot(s12, ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "l", pch = 19, col = 1, 
     main = "TEDM")
axis(side = 1, labels = F, col.ticks = NA)
axis(side = 2, labels = F, col.ticks = NA)
mtext("Trait 1", side = 1, line = 1)
mtext("Trait 2", side = 2, line = 1)

covhull2 <- chull(s12)
polygon(s12[covhull2, ], lwd = 2, border = 2)
#library(plot3D)
#scatter3D(x = s1[, 1], y = s1[, 2], z = s1[, 3], colvar = NULL, phi = 65, theta = 10, 
#          col = 2, bty = "b2")

#scatterplot3d::scatterplot3d(s1)


#### plotting minimum spanning tree (FEve)

x <- matrix(NA, nrow = 10, ncol = 2)
x <- apply(x, 2, function(x) runif(20, 0, 1))

library(ape)
mst(dist(x))
plot(mst(dist(x))
     )

plot(ComputeMST(x), ylab = "", xlab = "", xaxt = "n", pch = 19,
     yaxt = "n", bty = "l", col.pts = "red", cex = 3.1, 
     col.segts = "blue", lty = 5, lwd = 3)

mtext("Trait 1", side = 1, line = 2, cex = 2.1)
mtext("Trait 2", side = 2, line = 1, cex = 2.1)


## TED ploting
set.seed(1992)
par(mfrow = c(2, 2))
x <- matrix(NA, ncol = 2, nrow = 100)#rmvnorm(100, c(5, 5), 
    #         sigma = matrix(c(sig_11_start, 0.012, 0.012, sig_11_start), nrow = 2, ncol = 2)
    #         )

x <- apply(x, 2, function(x) runif(100, 0, 1))
plot(x, ylab = "", xlab = "", xaxt = "n", pch = 19,
     yaxt = "n", bty = "l", col = "red", cex = 3.1,)
mtext("Trait 1", side = 1, line = 2, cex = 2.1)
mtext("Trait 2", side = 2, line = 1, cex = 2.1)
mtext("Observed", side = 3, line = 2, cex = 2.1)

x_dist <- as.numeric(dist(x, diag = FALSE, upper = T))


x_sim <- base_traitmatrix("cube", ncols = 2, dsize = 100)
plot(x_sim, ylab = "", xlab = "", xaxt = "n", pch = 19,
     yaxt = "n", bty = "l", col = "blue", cex = 3.1,)
mtext("Trait 1", side = 1, line = 2, cex = 2.1)
mtext("Trait 2", side = 2, line = 1, cex = 2.1)
mtext("Hypothetical", side = 3, line = 2, cex = 2.1)

base_dist <- as.numeric(dist(x_sim, diag = FALSE, upper = T))

x_density <- density(x_dist)
base_density <- density(base_dist)

plot(x = x_density$x, y = x_density$y, col = "red", lwd = 5, lty = 1,
     ylab = "", xlab = "", bty = "l")
lines(base_density$x, base_density$y, col = "blue", lwd = 5, lty = 1)

mtext("Distance", side = 1, line = 3, cex = 2.1)
mtext("Kernel Density", side = 2, line = 2, cex = 2.1)





test <- 


 
sapply(test, ass)
 
 



x <- sceneraio_2[tobe_sampled[2],]

d <- dist_mean[tobe_sampled[2]]
b <- 0.1

delta_start <- c(0.1^2, 0.1^2)
test <- nleqslv(x = delta_start, fn = d_eq, jac = d_jb, control=list(btol=.01, delta="newton"))




 

