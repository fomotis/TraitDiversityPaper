######################## loading required packages #######
library(tidyverse)
library(data.table)
library(ggpubr)
library(mvtnorm)
library(gridExtra)

source("Codes/Trait Functions/TOP.R")
source("Codes/Trait Functions/distances.R")


############### plots in Figure 2

### scenario one
mu_start <- c(5, 10)
delta_mu <- seq(0, 10, by = 0.5)
rho_start <- 0.007
#equal variance
sig_11_start <- 0.85
cv_start <- sig_11_start*sig_11_start*rho_start

scenario1_plots <- vector("list", 4)
j <- 1
for(i in 1:length(delta_mu)) {
  mu_change <- c(mu_start[1] + delta_mu[i], mu_start[2])
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
            axis.line = element_line(size = 2),
            axis.title = element_text(size = 50, face = "bold"), 
            title = element_text(size = 40, face = "bold"),
            plot.title = element_text(hjust = -0.002)) +
      ggtitle("A")
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
            axis.line = element_line(size = 2),
            axis.title = element_text(size = 50, face = "bold"), 
            title = element_text(size = 40, face = "bold"),
            plot.title = element_text(hjust = -0.002)) +
      ggtitle(LETTERS[j])
    
    scenario1_plots[[j]] <- p2
    
  }
  
}

### expectations from scenario one
demonstration <- data.frame(shift = 0:10, Value = seq(0, 1, length = 11)
)

ltp <- c(Richness = "solid", 
         Evenness = "dashed", 
         Divergence = "dotted")

p4 <- demonstration %>% ggplot(aes(x = shift, y = Value)) + 
  geom_line(aes(x = seq(0, 3, length.out = 11), 
                y = seq(0.6, 0.8, length.out = 11)), 
            linetype = "solid", size = 3) + 
  #evennes
  geom_segment(aes(x = 3, y = 0.8, xend = 5.1, yend = 0.03), 
               size = 3, arrow = arrow(angle = 45), 
               linetype = "solid") +
  geom_text(aes(x = 4.9, y = 0.01), label = "Evenness", size = 15, angle = 360, fontface = "bold") +
  #richness
  geom_segment(aes(x = 0, y = 0.01, xend = 4.1, yend = 0.72), 
               size = 3, arrow = arrow(angle = 45), 
               linetype = "solid") +
  geom_text(aes(x = 4.5, y = 0.83), label = "Richness\nDivergence", size = 15, angle = 360, fontface = "bold") +
  #null line
  geom_segment(aes(x = 0, y = 0.2, xend = 7, yend = 0.1), 
              size = 3, arrow = arrow(angle = 45), 
              linetype = "dotted", color = "white") +
  #geom_text(aes(x = 8, y = 0.81), label = "Divergence", size = 15, angle = 360, fontface = "bold") +
  labs(x = "Location Shift", y = "Index value", linetype = "Index") + 
  theme_pubr() +
  scale_y_continuous(breaks = seq(0.1, 1.0, 0.1)) +
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  theme(axis.text = element_blank(),
        legend.position = 'top',
        axis.ticks = element_blank(), 
        axis.line = element_line(size = 2),
        axis.title = element_text(size = 50, face = "bold"), 
        title = element_text(size = 40, face = "bold"),
        plot.title = element_text(hjust = -0.002)) + 
  ggtitle("D") + scale_linetype_manual(values = ltp)

scenario1_plots[[4]] <- p4

final_plot_scenario1 <- grid.arrange(scenario1_plots[[1]], scenario1_plots[[2]], 
                                     scenario1_plots[[3]], scenario1_plots[[4]], 
                                     ncol = 4)
final_plot_scenario1 <- annotate_figure(final_plot_scenario1, 
                right = text_grob("Scenario One", face = "bold", size = 90, just = "center", rot = 90)
)
ggsave("Figures/final_plots_s1.png", plot = final_plot_scenario1, 
       device = "png", 
       height = 25, width = 35)


############### scenario two

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
            axis.line = element_line(size = 2),
            axis.title = element_text(size = 50, face = "bold"), 
            title = element_text(size = 40, face = "bold"),
            plot.title = element_text(hjust = -0.002)) + 
      ggtitle("E")
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
            axis.line = element_line(size = 2),
            axis.title = element_text(size = 50, face = "bold"), 
            title = element_text(size = 40, face = "bold"),
            plot.title = element_text(hjust = -0.002)) + 
      ggtitle(LETTERS[4+j])
    scenario2_plots[[j]] <- p2
    
  } else next
  
}

#expectation for scenario 2
p42 <- demonstration %>% ggplot(aes(x = shift, y = Value)) + 
  #evenness
  geom_segment(aes(x = 0, y = 0.8, xend = 8.5, yend = 0.15), 
               size = 3, arrow = arrow(angle = 45), 
               linetype = "solid") +
  geom_text(aes(x = 8.3, y = 0.09), label = "Evenness\nRichness", size = 15, fontface = "bold") +
  #null line
  geom_segment(aes(x = 0, y = 0.5, xend = 10, yend = 0.5), 
               size = 3, arrow = arrow(angle = 45), 
               linetype = "solid", color = "white") +
  #geom_text(aes(x = 4.05, y = 0.1), label = "Richness", size = 15, fontface = "bold") +
  #divergence
  geom_segment(aes(x = 0, y = 0.01, xend = 8.5, yend = 0.79), 
               size = 3, arrow = arrow(angle = 45), 
               linetype = "solid") +
  geom_text(aes(x = 8.3, y = 0.81), label = "Divergence", size = 15, fontface = "bold") +
  labs(x = "% Shifted", y = "Index value") + 
  theme_pubr() +
  scale_y_continuous(breaks = seq(0.1, 1.0, 0.1)) +
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  theme(axis.text = element_blank(),
        legend.position = 'top',
        axis.ticks = element_blank(), 
        axis.line = element_line(size = 2),
        axis.title = element_text(size = 50, face = "bold"), 
        title = element_text(size = 40, face = "bold"),
        plot.title = element_text(hjust = -0.002)) + 
  ggtitle("H")

scenario2_plots[[4]] <- p42

final_plot_scenario2 <- grid.arrange(scenario2_plots[[1]], scenario2_plots[[2]], 
                                     scenario2_plots[[3]], scenario2_plots[[4]], 
                                     ncol = 4)
final_plot_scenario2 <- annotate_figure(final_plot_scenario2, 
                                        right = text_grob("Scenario Two", face = "bold", size = 90, just = "center", rot = 90)
)
ggsave("Figures/final_plots_s2.png", plot = final_plot_scenario2, 
       device = "png", 
       height = 25, width = 35)


############# simulation scenario three

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
          legend.position = 'top',
          axis.ticks = element_blank(), 
          axis.line = element_line(size = 2),
          axis.title = element_text(size = 50, face = "bold"), 
          title = element_text(size = 40, face = "bold"),
          plot.title = element_text(hjust = -0.002)) + 
    ggtitle(LETTERS[8+j])
  scenario3_plots[[j]] <- p2
}
#expectation for scenario 3
p43 <- demonstration %>% ggplot(aes(x = shift, y = Value)) + 
  #evenness
  geom_segment(aes(x = 0, y = 0.9, xend = 10, yend = 0.9), 
               size = 2, arrow = arrow(angle = 45), 
               linetype = "dashed", color = "white") +
  #geom_text(aes(x = 3.5, y = 0.21), label = "Evenness", size = 9) +
  #richness
  geom_segment(aes(x = 0, y = 0.55, xend = 10, yend = 0.55), 
               size = 2, arrow = arrow(angle = 45), 
               linetype = "solid") +
  geom_text(aes(x = 5.0, y = 0.61), label = "Divergence\nEvenness", size = 15, fontface = "bold") +
  geom_text(aes(x = 5.0, y = 0.52), label = "Richness", size = 15, fontface = "bold") +
  #divergence
  geom_segment(aes(x = 0, y = 0.2, xend = 10, yend = 0.2), 
               size = 2, arrow = arrow(angle = 45), 
               linetype = "dotted", color = "white") +
  #geom_text(aes(x = 3.5, y = 0.41), label = "Divergence", size = 9) +
  labs(x = "Number", y = "Index value") + 
  theme_pubr() +
  scale_y_continuous(breaks = seq(0.0, 1.0, 0.1)) +
  scale_x_continuous(breaks = seq(1, 10, 1)) +
  theme(axis.text = element_blank(),
        legend.position = 'top',
        axis.ticks = element_blank(), 
        axis.line = element_line(size = 2),
        axis.title = element_text(size = 50, face = "bold"), 
        title = element_text(size = 40, face = "bold"),
        plot.title = element_text(hjust = -0.002)) + 
  ggtitle("M")

scenario3_plots[[4]] <- p43

final_plot_scenario3 <- grid.arrange(scenario3_plots[[1]], scenario3_plots[[2]], 
                                     scenario3_plots[[3]], scenario3_plots[[4]], 
                                     ncol = 4)

final_plot_scenario3 <- annotate_figure(final_plot_scenario3, 
                                        right = text_grob("Scenario Three", face = "bold", 
                                                          size = 90, just = "center", rot = 90)
)
ggsave("Figures/final_plots_s3.png", plot = final_plot_scenario3, 
       device = "png", 
       height = 25, width = 35)

final_plot_scenaos <- ggarrange(final_plot_scenario1, final_plot_scenario2, final_plot_scenario3, 
                                nrow = 3, align = "v")

ggsave("Figures/scenario123_plots.png", plot = final_plot_scenaos, 
       device = "png", 
       height = 47, width = 47)


########################## demonstration of the trait diversity indices
library(emstreeR)
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
     font = 2, cex = 5)

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
     font = 2, cex = 5)

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
     font = 2, cex = 5)

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
     type = "l", ylim = c(0.00, 2.7), xlim = c(-0.01, 1.2),
     bty = "l", col.pts = 1, 
     col.segts = 1, lty = 1, lwd = 6)
axis(side = 1, labels = F, col.ticks = NA, col = NA)
axis(side = 2, labels = F, col.ticks = NA, col = NA)
mtext("Distance", side = 1, line = 3, cex = 2.5)
mtext("Density", side = 2, line = 1, cex = 2.5)
lines(x = sim1_density$x, y = sim1_density$y, lty = 1, lwd = 6, 
      col = 2)
lines(x = sim2_density$x, y = sim2_density$y, lty = 1, lwd = 6,
      col = 3)
legend(x = 0.4, y = 2.7, legend = c("observed", "Discrete-Uniform (TED)", 
                                    "Continuous-Uniform (TEDM)"), 
       col = c(1, 2, 3), lty = c(1, 1, 1), lwd = rep(6, 3), horiz = FALSE, bty = "n", 
       cex = 4, ncol = 1,
       text.width = c(0.1, 0.1, 0.19),
       x.intersp = 0.16)
#text(x = 0.3, y = 1.55, "Observed", family = "sans", 
#     font = 2, cex = 5)
#text(x = 0.8, y = 1.05, "Discrete-Uniform (TED)", family = "sans", 
#     font = 2, cex = 5, col = 2)
text(x = 0.01, y = 2.7, "d", family = "sans", 
     font = 2, cex = 5)

# FDis and Rao
plot(sceneraio_2, pch = 19, 
     col = 1, cex = 10, bg = "blue", 
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
text(0.32, 0.87, labels = expression(y[i]), cex = 5, 
     font = 2, family = "sans", col = 1)
text(x = midpoint[1], y = midpoint[2], "c", font = 2, 
     family = "sans", col = 1, cex = 5)
mtext("Trait 1", side = 1, line = 3, cex = 2.5)
mtext("Trait 2", side = 2, line = 1, cex = 2.5)
text(x = min(sceneraio_2[,1]), y = 1.1, "e", family = "sans", 
     font = 2, cex = 5)
dev.off()

####################### Plots in Appendix A1

x1 <- rep(1:5,each=5)
x2 <- rep(1:5,times=5)
X <- data.frame(Trait1 = x1, Trait2 = x2)

plot1_app1 <- X  %>% ggplot(aes(x = Trait1, y = Trait2, group = 1)) + 
  geom_point(color = "black", size = 10) + 
  geom_point(aes(x = 3, y = 3), color = "red", size = 10) +
  labs( x = "Trait 1", y = "Trait 2") + 
  theme_bw() + 
  theme(  axis.title = element_text(size = 50, face = "bold"), 
          axis.text = element_text(size = 40, face = "bold"),
          title = element_text(size = 40, face = "bold")
          ) + 
  geom_segment( aes(x = 3.01, y = 3.01, xend = 4.95, yend = 4.95), 
                arrow = arrow(angle = 45, length = unit(0.4, "inches")),
                size = 5, linetype = "solid"
              ) + 
  ggtitle("A")

appendix1 <- lapply(seq(3, 5, by = 0.05), function(i) {
  
   X[13, 1] <- i
   X[13, 2] <- i
   print(X[13, ])
   position = eucldis(c(3, 3) - c(i, i))
   data.frame(Distance_to_start = position, TOP = TOP(X)$TA)  
}) %>% do.call(rbind.data.frame, .)

plot2_app1 <- appendix1 %>% ggplot(aes(x = Distance_to_start, y = TOP, group = 1)) + 
  geom_point(size = 10, color = "black") + 
  labs(x = "Distance to c(3, 3)", y = "TOP") + 
  theme_bw() + 
  theme(  axis.title = element_text(size = 50, face = "bold"),
          title = element_text(size = 40, face = "bold"),
          axis.text = element_text(size = 40, face = "bold")
  ) + scale_x_continuous(
    breaks = seq(0, 2.85, by = 0.2)
    ) + 
  ggtitle("B")
appendix1_plots <-  ggarrange(plot1_app1, plot2_app1, ncol = 2, hjust = 0.1, align = 'h')
ggsave(filename = "Figures/TOP_Appendix1.png", plot = appendix1_plots, 
       width = 40, height = 20)


################## plots in Appendix A2

Tsdata1 <- matrix(NA, ncol = 2, nrow = 50)
Tsdata1 <-  apply(Tsdata1, 2, function(x) runif(50, 0, 1))

bmat1 <- as.data.frame(base_traitmatrix(basetype = "cube", ncols = 2, 
                          dsize = 50))
bmat2 <- as.data.frame(Tsdata1)
names(bmat1) <- names(bmat2) <- c("Trait 1", "Trait 2")
bdata <- data.frame(rbind(bmat1, bmat2), Type = c(rep("Discrete-Uniform (TED)", 
                                                      nrow(bmat1)), rep("Continuous-Uniform (TEDM)", nrow(bmat2))))

tplot1 <- bdata %>%
  ggplot(aes(Trait.1, Trait.2)) + 
  geom_point(size = 10) + 
  theme_bw() +
  facet_wrap(.~Type) +
  labs(x = "Trait 1", y = "Trait 2") +
  theme(strip.text = element_text(size = 60, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_text(size = 50, face = "bold"), 
        title = element_text(size = 40, face = "bold"),
        plot.title = element_text(hjust = -0.002),
        panel.spacing = unit(5, "lines"))
  
bmat1_dist <- as.matrix(dist(bmat1, diag = F, upper = T))
bmat1_dist <- as.numeric(bmat1_dist[base::upper.tri(bmat1_dist)])
bmat1_density <- density(bmat1_dist, from = min(bmat1_dist), to = max(bmat1_dist))

bmat2_dist <- as.matrix(dist(bmat2, diag = F, upper = T))
bmat2_dist <- as.numeric(bmat2_dist[base::upper.tri(bmat2_dist)])
bmat2_density <- density(bmat2_dist, from = min(bmat2_dist), to = max(bmat2_dist))

bdata2 <- data.frame(X = c(bmat1_density$x, bmat2_density$x), 
           Y = c(bmat1_density$y, bmat2_density$y),
           Type = c(rep("Discrete-Uniform (TED)", length(bmat1_density$x)), rep("Continuous-Uniform (TEDM)", length(bmat2_density$x)))
           )
tplot2 <- bdata2 %>%
  ggplot(aes(X, Y)) + 
  geom_line(size = 10) + 
  theme_bw() +
  facet_wrap(.~Type) +
  labs(x = "Distance", y = "Density") +
  theme(strip.text = element_text(size = 60, face = "bold"),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_text(size = 50, face = "bold"), 
        title = element_text(size = 40, face = "bold"),
        plot.title = element_text(hjust = -0.002),
        panel.spacing = unit(5, "lines"))

a2plot <- ggarrange(tplot1, tplot2, nrow = 2, hjust = 0.1, align = 'h')
ggsave(filename = "Figures/tedm_motivation.png", plot = a2plot, 
       width = 30, height = 30)

#### plots in Appendix A3

#middle and distance to the center
xbar <- apply(Tsdata1, 2, mean)
dist_mean <- apply(Tsdata1 , 1, function(x) {
  eucldis(x - xbar)
})

gendata <- function(delete) {
  tobe_sampled <- which(dist_mean <= quantile(dist_mean, probs = delete))
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
  
  #distance distribution for Tsdata12
  orig_dist <- as.matrix(dist(Tsdata12, diag = F, upper = T))
  orig_dist <- as.numeric(orig_dist[base::upper.tri(orig_dist)])
  orig_density <- density(orig_dist, from = min(orig_dist), to = max(orig_dist))
  
  #distance distribution for TED
  bmat1 <- base_traitmatrix(basetype = "cube", ncols = ncol(Tsdata12), 
                            dsize = nrow(Tsdata12))
  bmat1_dist <- as.matrix(dist(bmat1, diag = F, upper = T))
  bmat1_dist <- as.numeric(bmat1_dist[base::upper.tri(bmat1_dist)])
  bmat1_density <- density(bmat1_dist, from = min(bmat1_dist), to = max(bmat1_dist))
  
  #distance distribution for TEDM
  bmat2 <- base_traitmatrix_new(Tsdata12, gen_matrix = gen_matrix, 
                                rescale = T, parallel = F)
  bmat2_dist <- as.matrix(dist(bmat2, diag = F, upper = T))
  bmat2_dist <- as.numeric(bmat2_dist[base::upper.tri(bmat2_dist)])
  bmat2_density <- density(bmat2_dist, from = min(bmat2_dist), to = max(bmat2_dist))
  
  appendix3 <- data.frame(X = c(orig_density$x, bmat1_density$x, bmat2_density$x),
                          Y = c(orig_density$y, bmat1_density$y, bmat2_density$y),
                          Distance_Distribution = c(rep("Observed", length(orig_density$x)), 
                                                    rep("TED", length(bmat1_density$x)), 
                                                    rep("TEDM", length(bmat2_density$x)) 
                                                    )
  )
  appendix3 <- cbind(appendix3, Delete = rep(paste0("Shifted = ", delete*100, "%"), nrow(appendix3)))
  return(appendix3)
}

appendix3 <- rbind(gendata(delete = 0),
            gendata(delete = 0.1),
            gendata(delete = 0.3)
)

app3 <- appendix3 %>% ggplot(aes(x = X, y = Y, group = Distance_Distribution, color = Distance_Distribution)) + 
  geom_line(size = 10) + 
  theme_bw() + 
  facet_wrap(. ~ Delete, scales = "fixed") +
  labs(x = "Distance", y = "Density", color = "Distribution") + 
  theme(  strip.text = element_text(size = 60, face = "bold"),
          axis.line.y = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text = element_text(size = 60, face = "bold"),
          axis.title = element_text(size = 70, face = "bold"), 
          legend.text = element_text(size = 60, face = "bold"),
          legend.title = element_text(size = 60, face = "bold"), 
          legend.position = "top", 
          legend.direction = "horizontal",
          panel.spacing = unit(5, "lines")
  )

ggsave(filename = "Figures/ted_explanation2.png", plot = app3, 
       width = 30, height = 30)

############ plots in Appendix A4


########################## Plotting simulation Results ###################################
load("RData/TwoDimensionSimulationResults_results1and2.RData")
load("RData/TwoDimensionSimulation_SampleSize_new_Results.RData")

max_of_maxes <- function(y) {
   
   maxes <- lapply(y, function(x) { apply(x, MARGIN = 2, max)}) %>% do.call(rbind, .)
   return(apply(maxes, 2, max))
 }

zero_one <- function(x, group_var, group) {
  
  y <- rep(NA, length(x))
  for(i in group) {
    y[group_var == i] <- x[group_var == i] / max(x[group_var == i])
  }
  
  return(y)
  
}

t_func <- function(x, rescale = FALSE, ms = NULL) {
  
  
  if(rescale  == TRUE & !is.null(ms)) {
    
    x_new <- t(apply(x, 1, function(x) x / ms))
    
  } else {
    x_new <- x
  }
  
  xm <- apply(x_new ,2, mean)
  xsd <- apply(x_new ,2, sd)
  return(data.frame(M = xm, SD = xsd))
}

t_func2 <- function(vec, tims = NULL, each = NULL) {
  
  if(is.null(each)) {
    rt <- rep(vec, tims)
  } else if(is.null(tims)) {
    rt <- rep(vec, each = each)
  } else {
    rt <- rep(vec, tims, each = each)
  }
  return(rt)
}

################# scenario one ##################
results_s1 <- lapply(Results_comm, t_func, rescale = FALSE) %>% 
  do.call(rbind.data.frame,.) %>%
  mutate(
          Index = t_func2(c("FRIC","TOP", "TOPM", "FEVE", "TED", "TEDM", "RAO", "FDIS"), tims = 21), 
          Manipulate = t_func2(delta_mu, each = 8)
         ) %>% 
  mutate(
      Index_Type = case_when(
      Index %in% c("FRIC", "TOP", "TOPM") ~ "Richness",
      Index %in% c("FEVE", "TED", "TEDM") ~ "Evenness",
      Index %in% c("RAO", "FDIS") ~ "Divergence"
    )
  ) %>% 
  mutate(
    M_rescaled = zero_one(M, group_var = Index, group = unique(Index)),
    SD_rescaled = zero_one(SD,  group_var = Index, group = unique(Index)),
    LCL = M - SD,
    UCL = M + SD
  ) %>% 
  mutate(
    
    LCL_rescaled = zero_one(LCL, group_var = Index, group = unique(Index)),
    UCL_rescaled = zero_one(UCL, group_var = Index, group = unique(Index))
    
  )

results_s1$M_rescaled[results_s1$Index == "TOP"] <- jitter(results_s1$M_rescaled[results_s1$Index == "TOP"], amount = 0.008)

#expectations for evenness
#exp_evenness <- data.frame(M = c(seq(0.8, 1, length = 7), seq(1, 0.25, length = 21-7)), SD = rep(0, 21), Index = rep("Expected", 21), 
#               Manipulate = delta_mu, Index_Type = "Evenness",  
#               M_rescaled = c(seq(0.8, 1, length = 7), seq(1, 0.25, length = 21-7)))

#expectaion for divergence
#exp_divergence <- data.frame(M = seq(0.1, 1, length = 21), SD = rep(0, 21), Index = rep("Expected", 21), 
#              Manipulate = delta_mu, Index_Type = "Divergence",  
#              M_rescaled = seq(0.1, 1, length = 21))

#expectaion for richness
#exp_richness <- data.frame(M = seq(0.1, 1, length = 21), SD = rep(0, 21), Index = rep("Expected", 21), 
#              Manipulate = delta_mu, Index_Type = "Richness",  
#              M_rescaled = seq(0.1, 1, length = 21))

#plot_data1 <- rbind(results_s1, exp_evenness, exp_divergence, exp_richness)

#bp1 <- c(FDIS = "#1B9E77", RAO = "#D95F02", 
#                 FRIC = "#7570B3", TOP = "#E7298A", TOPM = "#66A61E", 
#                 FEVE = "#E6AB02", TED = "#A6761D", TEDM = "#CC79A7", 
#         Expected = "#000000")

results_s1 %>% ggplot(aes(x = Manipulate, y = M_rescaled, group = Index, color = Index)) + 
  geom_line(size = 3) +
  facet_wrap(. ~ Index_Type) + 
  scale_color_manual(values = bp1) + 
  theme_bw()  + scale_y_continuous(breaks = seq(0.6, 1.0, by = 0.05))
  

######### scenario two ############
results_s2 <- lapply(Results_comm22, t_func, rescale = FALSE) %>% 
  do.call(rbind.data.frame, .)  %>%
  mutate(
    Index = t_func2(c("FRIC","TOP", "TOPM", "FEVE", "TED", "TEDM", "RAO", "FDIS"), tims = 6),
    Manipulate = t_func2(c(0, delete), each = 8)*100
  ) %>% 
  mutate(
    Index_Type = case_when(
      Index %in% c("FRIC", "TOP", "TOPM") ~ "Richness",
      Index %in% c("FEVE", "TED", "TEDM") ~ "Evenness",
      Index %in% c("RAO", "FDIS") ~ "Divergence"
    )
  )  %>% 
  mutate(
    M_rescaled = zero_one(M, group_var = Index, group = unique(Index)),
    LCL = M - SD,
    UCL = M + SD
  ) %>% 
  mutate(
    
    LCL_rescaled = zero_one(LCL, group_var = Index, group = unique(Index)),
    UCL_rescaled = zero_one(UCL, group_var = Index, group = unique(Index))
  )
results_s2$M_rescaled[results_s2$Index == "TOP"] <- jitter(results_s2$M_rescaled[results_s2$Index == "TOP"], amount = 0.008)

#expectations for evenness
# exp_evenness <- data.frame(M = seq(1, 0.1, length = 5), SD = rep(0, 5), Index = rep("Expected", 5), 
#                            Manipulate = delete, Index_Type = "Evenness",  
#                            M_rescaled = seq(1, 0, length = 5))
# 
# #expectaion for divergence
# exp_divergence <- data.frame(M = seq(0.1, 1, length = 5), SD = rep(0, 5), Index = rep("Expected", 5), 
#                              Manipulate = delete, Index_Type = "Divergence",  
#                              M_rescaled = seq(0, 1, length = 5))
# 
# #expectaion for richness
# exp_richness <-  data.frame(M = seq(1, 0.1, length = 5), SD = rep(0, 5), Index = rep("Expected", 5), 
#                             Manipulate = delete, Index_Type = "Richness",  
#                             M_rescaled = seq(1, 0.1, length = 5))
# 
# plot_data2 <- rbind(results_s2, exp_evenness, exp_divergence, exp_richness)

results_s2 %>% ggplot(aes(x = Manipulate, y = M_rescaled, group = Index, color = Index)) + 
  geom_line(size = 3) +
  facet_wrap(. ~ Index_Type) + 
  scale_color_manual(values = bp1) + 
  theme_bw() + scale_y_continuous(breaks = seq(0.6, 1.0, by = 0.015))

######### scenario three ############
results_s3 <- lapply(Results_samp, t_func, rescale = FALSE) %>% 
  do.call(rbind.data.frame, .)  %>%
  mutate(
    Index = t_func2(c("FRIC","TOP", "TOPM", "FEVE", "TED", "TEDM", "RAO", "FDIS"), tims = 9),
    Manipulate = t_func2(n_row/1000, each = 8)
  ) %>% 
  mutate(
    Index_Type = case_when(
      Index %in% c("FRIC", "TOP", "TOPM") ~ "Richness",
      Index %in% c("FEVE", "TED", "TEDM") ~ "Evenness",
      Index %in% c("RAO", "FDIS") ~ "Divergence"
    )
  ) %>% 
  mutate(
    M_rescaled = zero_one(M, group_var = Index, group = unique(Index)),
    LCL = M - SD,
    UCL = M + SD
  ) %>% 
  mutate(
    
    LCL_rescaled = zero_one(LCL, group_var = Index, group = unique(Index)),
    UCL_rescaled = zero_one(UCL, group_var = Index, group = unique(Index))
  )

results_s3 %>% ggplot(aes(x = Manipulate, y = M_rescaled, group = Index, color = Index)) + 
  geom_line(size = 3) +
  facet_wrap(. ~ Index_Type) + 
  scale_color_manual(values = bp1) + 
  theme_bw()


### bind the three results together

combined_result <- rbind(results_s1, results_s2, results_s3) %>% 
  mutate(
    Scenario = c(rep("Scenario One", 168), 
                 rep("Scenario Two", 48),
                 rep("Scenario Three", 72)
                 )
  )

combined_result$Scenario <- factor(combined_result$Scenario, 
                                   levels = c("Scenario One", "Scenario Two", "Scenario Three")
                                     )
combined_result$Index <- factor(combined_result$Index, 
                                levels = c("FDIS", "RAO", 
                                           "FEVE", "TED", "TEDM", 
                                           "FRIC", "TOP", "TOPM")
                                )

results_main <- combined_result %>% ggplot(aes(x = Manipulate, y = M_rescaled, group = Index, color = Index)) + 
  geom_line(size = 8) +
  facet_wrap(Scenario ~ Index_Type, scales = "free_x") + 
  scale_color_manual(values = bp1) + 
  theme_bw() + 
  labs(y = "Index Value", x = "")+
  theme(strip.text = element_text(size = 60, face = "bold"),
        axis.line.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 60, face = "bold"),
        axis.title = element_text(size = 70, face = "bold"), 
        legend.text = element_text(size = 60, face = "bold"),
        legend.title = element_text(size = 60, face = "bold"), 
        legend.position = "top", 
        legend.direction = "horizontal",
        panel.spacing = unit(5, "lines"))

ggsave(filename = "Figures/results_main.png", plot =results_main, 
       width = 30, height = 30)

######### approach two

p1 <- combined_result %>% filter(Scenario == "Scenario One") %>% ggplot(aes(x = Manipulate, y = M_rescaled, group = Index, color = Index)) + 
  geom_line(size = 8) +
  facet_wrap(Scenario ~ Index_Type, scales = "free_x") + 
  scale_color_manual(values = bp1) + 
  theme_bw() + 
  labs(y = "", x = "Location shift")+
  theme(strip.text = element_text(size = 60, face = "bold"),
        axis.line.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 60, face = "bold"),
        axis.title = element_text(size = 70, face = "bold"), 
        legend.text = element_text(size = 60, face = "bold"),
        legend.title = element_text(size = 60, face = "bold"), 
        legend.position = "top", 
        legend.direction = "horizontal") + 
  theme(panel.spacing = unit(5, "lines"))

p2 <- combined_result %>% filter(Scenario == "Scenario Two") %>% 
  ggplot(aes(x = Manipulate*100, y = M_rescaled, group = Index, color = Index)) + 
  geom_line(size = 8) +
  facet_wrap(Scenario ~ Index_Type, scales = "free_x") + 
  scale_color_manual(values = bp1) + 
  theme_bw() + 
  labs(y = "Index Value", x = "Percentage shifted outward")+
  theme(strip.text = element_text(size = 60, face = "bold"),
        axis.line.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 60, face = "bold"),
        axis.title = element_text(size = 70, face = "bold"), 
        legend.text = element_blank(),
        legend.title = element_blank(), 
        legend.position = "none") + 
  theme(panel.spacing = unit(5, "lines"))

p3 <- combined_result %>% filter(Scenario == "Scenario Three") %>% 
  ggplot(aes(x = Manipulate, y = M_rescaled, group = Index, color = Index)) + 
  geom_line(size = 8) +
  facet_wrap(Scenario ~ Index_Type, scales = "free_x") + 
  scale_color_manual(values = bp1) + 
  theme_bw() + 
  labs(y = "", x = "Number of individuals (in thousands)")+
  theme(strip.text = element_text(size = 60, face = "bold"),
        axis.line.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 60, face = "bold"),
        axis.title = element_text(size = 70, face = "bold"), 
        legend.text = element_blank(),
        legend.title = element_blank(), 
        legend.position = "none") + 
  theme(panel.spacing = unit(5, "lines"))

results_main <-  ggarrange(p1, p2, p3, nrow = 3, hjust = 0.1, align = 'hv')
ggsave(filename = "Figures/results_main.png", plot =results_main, 
       width = 30, height = 30)


######### Results with variability #################

### scenario 1
ms1 <- max_of_maxes(Results_comm)
results_s1 <- lapply(Results_comm, t_func, rescale = TRUE, ms = ms1) %>% 
  do.call(rbind.data.frame,.) %>%
  mutate(
    Index = t_func2(c("FRIC","TOP", "TOPM", "FEVE", "TED", "TEDM", "RAO", "FDIS"), tims = 21), 
    Manipulate = t_func2(delta_mu, each = 8)
  ) %>% 
  mutate(
    Index_Type = case_when(
      Index %in% c("FRIC", "TOP", "TOPM") ~ "Richness",
      Index %in% c("FEVE", "TED", "TEDM") ~ "Evenness",
      Index %in% c("RAO", "FDIS") ~ "Divergence"
    )
  ) %>% 
  mutate(
    LCL = M - SD,
    UCL = M + SD
  ) 

### scenario 2
ms2 <- max_of_maxes(Results_comm22)
results_s2 <- lapply(Results_comm22, t_func, rescale = TRUE, ms = ms2) %>% 
  do.call(rbind.data.frame, .)  %>%
  mutate(
    Index = t_func2(c("FRIC","TOP", "TOPM", "FEVE", "TED", "TEDM", "RAO", "FDIS"), tims = 6),
    Manipulate = t_func2(c(0, delete), each = 8)*100
  ) %>% 
  mutate(
    Index_Type = case_when(
      Index %in% c("FRIC", "TOP", "TOPM") ~ "Richness",
      Index %in% c("FEVE", "TED", "TEDM") ~ "Evenness",
      Index %in% c("RAO", "FDIS") ~ "Divergence"
    )
  ) %>% mutate(
    LCL = M - SD,
    UCL = M + SD
  )

#### scenario 3
ms3 <- max_of_maxes(Results_samp)
results_s3 <- lapply(Results_samp, t_func, rescale = TRUE, ms = ms3) %>% 
  do.call(rbind.data.frame, .)  %>%
  mutate(
    Index = t_func2(c("FRIC","TOP", "TOPM", "FEVE", "TED", "TEDM", "RAO", "FDIS"), tims = 9),
    Manipulate = t_func2(n_row/1000, each = 8)
  ) %>% 
  mutate(
    Index_Type = case_when(
      Index %in% c("FRIC", "TOP", "TOPM") ~ "Richness",
      Index %in% c("FEVE", "TED", "TEDM") ~ "Evenness",
      Index %in% c("RAO", "FDIS") ~ "Divergence"
    )
  ) %>% mutate(
    LCL = M - SD,
    UCL = M + SD
  )

#combine results
combined_result <- rbind(results_s1, results_s2, results_s3) %>% 
  mutate(
    Scenario = c(rep("Scenario One", 168), 
                 rep("Scenario Two", 48),
                 rep("Scenario Three", 72)
    )
  )
combined_result$Scenario <- factor(combined_result$Scenario, 
                                   levels = c("Scenario One", "Scenario Two", "Scenario Three")
)
combined_result$Index <- factor(combined_result$Index, 
                                levels = c("FDIS", "RAO", 
                                           "FEVE", "TED", "TEDM", 
                                           "FRIC", "TOP", "TOPM")
)

p1s <- combined_result %>% filter(Scenario == "Scenario One") %>% 
  ggplot(aes(x = Manipulate, y = M, ymin =  LCL, ymax =  UCL, group = Index, color = Index)) + 
  geom_errorbar(width = 0.1, size = 4, position = position_dodge(0.5)) +
  geom_point(size = 9, position = position_dodge(0.5)) +
  facet_wrap(Scenario ~ Index_Type, scales = "free") + 
  scale_color_manual(values = bp1) + 
  theme_bw() + 
  labs(y = "", x = "Location shift") +
  theme(strip.text = element_text(size = 60, face = "bold"),
        axis.line.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 60, face = "bold"),
        axis.title = element_text(size = 70, face = "bold"), 
        legend.text = element_text(size = 60, face = "bold"),
        legend.title = element_text(size = 60, face = "bold"), 
        legend.position = "top", 
        legend.direction = "horizontal") + 
  theme(panel.spacing = unit(3, "lines")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2))

p2s <- combined_result %>% filter(Scenario == "Scenario Two") %>% 
  ggplot(aes(x = Manipulate, y = M, ymin =  LCL, ymax =  UCL, group = Index, color = Index)) + 
  geom_errorbar(width = 0.1, size = 4, position = position_dodge(0.5)) +
  geom_point(size = 9, position = position_dodge(0.5)) +
  facet_wrap(Scenario ~ Index_Type, scales = "free") + 
  scale_color_manual(values = bp1) + 
  theme_bw() + 
  labs(y = "", x = "Percentage shifted outward") +
  theme(strip.text = element_text(size = 60, face = "bold"),
        axis.line.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 60, face = "bold"),
        axis.title = element_text(size = 70, face = "bold"), 
        legend.text = element_text(size = 60, face = "bold"),
        legend.title = element_text(size = 60, face = "bold"), 
        legend.position = "none", 
        legend.direction = "horizontal") + 
  theme(panel.spacing = unit(3, "lines"))

p3s <- combined_result %>% filter(Scenario == "Scenario Three") %>% 
  ggplot(aes(x = Manipulate, y = M, ymin =  LCL, ymax =  UCL, group = Index, color = Index)) + 
  geom_errorbar(width = 0.1, size = 4, position = position_dodge(0.5)) +
  geom_point(size = 9, position = position_dodge(0.5)) +
  facet_wrap(Scenario ~ Index_Type, scales = "free") + 
  scale_color_manual(values = bp1) + 
  theme_bw() + 
  labs(y = "", x = "Number of individuals (in thousands)") +
  theme(strip.text = element_text(size = 60, face = "bold"),
        axis.line.y = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 60, face = "bold"),
        axis.title = element_text(size = 70, face = "bold"), 
        legend.text = element_text(size = 60, face = "bold"),
        legend.title = element_text(size = 60, face = "bold"), 
        legend.position = "none", 
        legend.direction = "horizontal") + 
  theme(panel.spacing = unit(3, "lines")) + scale_x_continuous(breaks = seq(0, 3.5, by = 0.25))

results_mains <-  ggarrange(p1s, p2s, p3s, nrow = 3, hjust = 0.1, align = 'hv')
ggsave(filename = "Figures/results_mains.png", plot = results_mains, 
       width = 30, height = 30)
