# LC6 downstream neurons analysis

# - target color
pal_target <- pnw_palette("Bay",9)
pal_target_a <- adjustcolor(pal_target, alpha.f = 0.6)

# Figure S7, target neuron RF wt by synapses ----------------------------------------------------------------------


# PLOT,  Gaussian around each com with synap# as height,  cp data via binning
ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])
plvl <- list()
simdata <- list()
simdata_df <- list()
# LC6_tar_median <- list()
# for (j in 1:length(conn_target)) {
#   LC6_tar_median[[j]] <- (quantile(conn_target[[j]]$tofrom_glo, c(0.0)))
# }
mat_names <- c(paste("biL_", biL_skid, sep = ""), 
               paste("biR_", biR_skid, sep = ""), 
               paste("ipsi_", ipsi_skid, sep = ""))

n_lvl <- 11
breaks_3 <- seq(0,70,length.out = n_lvl)
getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
pal_tar <- getPalette(n_lvl - 1)

x_tick <- c(16, 16+45, 16+90, 16+135, 16+180)
y_tick <- c(0, 45, 90, 135, 180)

for (j in 1:length(neu_target)) {
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(neu)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- conn_target[[j]][k,"tofrom_glo"]
    # if (A >= LC6_tar_median[[j]]) { # selected neuron
      grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
      # grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/(r_xy/2)^2 - (x[2]-y0)^2/(r_xy/2)^2))
    # }
  }
  
  grid_Gaussian_cut <- grid_Gaussian %>%
    as.data.frame()
  
  # -- whole range no binning
  simdata[[j]] <- grid_Gaussian_cut
  simdata_df[[j]] <- simdata[[j]]
  colnames(simdata_df[[j]]) <- c("x","y","z")
  # simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  # simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
  # simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  # simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  simdata_df[[j]]$z <- simdata_df[[j]]$z / mean(head(sort(simdata_df[[j]]$z, decreasing = T)))
  simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  plvl[[j]] <-
    ggplot(simdata_df[[j]], aes(x, y, z = z)) +
    geom_raster(aes(fill = equalSpace), interpolate = F) +
    scale_fill_manual(values = pal_tar) +
    geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE) +
    geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
    annotate("text", x = -5, y = 185, label = "9°") +
    geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
    annotate("text", x = -15, y = 175.5, label = "9°") +
    geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
    geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
    geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
    scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
    scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
    labs(title = mat_names[j]) +
    theme_void()+
    coord_fixed(ratio = 1)
  
  plvl[[j]] <- plvl[[j]] +
    geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.7), color = "blue", alpha = 1, lwd = 1)

}

windows(record = F, width = 8, height = 8)
plvl[[1]]
windows(record = F, width = 8, height = 8)
plvl[[2]]
windows(record = F, width = 8, height = 8)
plvl[[3]]
windows(record = F, width = 8, height = 8)
plvl[[4]]
windows(record = F, width = 8, height = 8)
plvl[[5]]
windows(record = F, width = 8, height = 8)
plvl[[6]]
windows(record = F, width = 8, height = 8)
plvl[[7]]
windows(record = F, width = 8, height = 8)
plvl[[8]]
windows(record = F, width = 8, height = 8)
plvl[[9]]



#  Figure 7C, sum up cell type ------------------------------------------------------------------------------------

ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])
plvl <- list()
simdata <- list()
simdata_df <- list()
mat_names <- c("bi_all", "ipsi_all")

conn_target_agglo_sum <- list()
conn_target_agglo_sum[[1]] <- conn_target[[1]] + conn_target[[2]] + conn_target[[3]] + conn_target[[4]]
conn_target_agglo_sum[[2]] <- conn_target[[5]]+conn_target[[6]]+conn_target[[7]]+conn_target[[8]]+conn_target[[9]]

for (j in 1:length(conn_target_agglo_sum)) {  
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(neu)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- conn_target_agglo_sum[[j]][k,"tofrom_glo"]
    # if (A > 0) { # selected neuron
      # grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
      grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/(r_xy/3*2)^2 - (x[2]-y0)^2/(r_xy/3*2)^2))
    # }
  }
  simdata[[j]] <- grid_Gaussian
  simdata_df[[j]] <- simdata[[j]]
  colnames(simdata_df[[j]]) <- c("x","y","z")
  # simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  # simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
  # simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  # simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  simdata_df[[j]]$z <- simdata_df[[j]]$z / mean(head(sort(simdata_df[[j]]$z, decreasing = T)))
  simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  plvl[[j]] <-
    ggplot(simdata_df[[j]], aes(x, y, z = z)) +
    geom_raster(aes(fill = equalSpace), interpolate = F) +
    scale_fill_manual(values = pal_tar) +
    geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE) +
    geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
    annotate("text", x = -5, y = 185, label = "9°") +
    geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
    annotate("text", x = -15, y = 175.5, label = "9°") +
    geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 1) +
    geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 1) +
    geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 1) +
    scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
    scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
    labs(title = mat_names[j]) +
    theme_void() +
    coord_fixed(ratio = 1)
  
  plvl[[j]] <- plvl[[j]] +
    geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 1)

}

windows(record = F, width = 8, height = 8)
plvl[[1]] 
windows(record = F, width = 8, height = 8)
plvl[[2]] 


# target neurons syn histo glo  ---------------------------------------------------------------------------

N_gp <- 11

# combine conn data
conn_LC6_tar_bi <- rbind(conn_LC6_tar[[1]], conn_LC6_tar[[2]], conn_LC6_tar[[3]], conn_LC6_tar[[4]])
gp_x <- factor(conn_LC6_tar_bi[,"gp_x"], labels = seq(1,N_gp-1, by = 1))
conn_LC6_tar_bi[,"gp_x"] <- gp_x
conn_LC6_tar_ipsi <- rbind(conn_LC6_tar[[5]], conn_LC6_tar[[6]], conn_LC6_tar[[7]], conn_LC6_tar[[8]], conn_LC6_tar[[9]])
gp_x <- factor(conn_LC6_tar_ipsi[,"gp_x"], labels = seq(1,N_gp-1, by = 1))
conn_LC6_tar_ipsi[,"gp_x"] <- gp_x

# PLOT, target with synapse
nopen3d()
plot3d(neu_target[[5]], col = 'lightblue', alpha = 0.7)
points3d(conn_LC6_tar[[5]][,c('x','y','z')], col = 'blue')
plot3d(neu_target[[6]], col = 'pink', alpha = 0.7)
points3d(conn_LC6_tar[[6]][,c('x','y','z')], col = 'red')


# histogram of num of synp in each compartment
LC6_ipsi_dol <- matrix(ncol = 5, nrow = 10)
for (j in 1:5) {
  for (k in 1:10) {
    LC6_ipsi_dol[k,j] <- sum(LC6_wt[[4+j]][[k]][,2])
  }
}


# Figure S8B, bar plot -- obsolete
num_tot <- colSums(LC6_ipsi_dol)
for (j in 1:5) {
  LC6_ipsi_dol[,j] <- LC6_ipsi_dol[,j]/num_tot[j]
}
dev.new()
barplot(LC6_ipsi_dol[,c(3,2,4,5,1)], main = 'LC6 to ipsi in dolphin compartments', ylim = c(-0.1, 1.2), col = dolphin_col )
# legend(x = 5, y = 800, legend = c(paste(seq(1,10))), cex = 1.3, fill = dolphin_col, horiz = F)
text(x = seq(0.7,5.5, length.out = 5), y = rep(1.1,5), labels = num_tot[c(3,2,4,5,1)])
text(x = seq(0.7,5.5, length.out = 5), y = rep(-0.05,5), labels = ipsi_skid[c(3,2,4,5,1)])



# target RF in 10 compartments ------------------------------------------------------------------------------------

# separate targets
LC6_wt_biL <- list() # first 2 are biL, 2 biR, last 5 are ipsi
LC6_wt_biR <- list()
LC6_wt_bi <- list()
LC6_wt_ipsi <- list()
for (j in 1:(length(glo_div)-1)) {
  LC6_wt_biL[[j]] <- LC6_wt[[1]][[j]] + LC6_wt[[2]][[j]]
  LC6_wt_biR[[j]] <- LC6_wt[[3]][[j]] + LC6_wt[[4]][[j]]
  LC6_wt_bi[[j]] <- LC6_wt[[1]][[j]] + LC6_wt[[2]][[j]] + LC6_wt[[3]][[j]] + LC6_wt[[4]][[j]]
  LC6_wt_ipsi[[j]] <- LC6_wt[[5]][[j]] + LC6_wt[[6]][[j]] + LC6_wt[[7]][[j]] + LC6_wt[[8]][[j]] + LC6_wt[[9]][[j]]
}

# num of synapse per division
N_syn <- matrix(nrow = length(glo_div)-1, ncol = 3)
for (j in 1:(length(glo_div)-1)) {
  N_syn[j,1] <- colSums(LC6_wt_biL[[j]])[2]
  N_syn[j,2] <- colSums(LC6_wt_biR[[j]])[2]
  N_syn[j,3] <- colSums(LC6_wt_ipsi[[j]])[2]
}

# com of each division
com_syn <- cbind(rowMeans(range_syn[,1:2]), rowMeans(range_syn[,3:4]))


# Figure S8A
# -- bi all
plvl <- list()
grid_Gaussian_ls <- list()
grid_Gaussian_sum <- matrix(ncol = 4)
colnames(grid_Gaussian_sum) <- c("X","Y","Z","gp")
plvl_sum <- ggplot()
pal_dolphin <- c()
gprange <- list()
for (j in 1:(length(glo_div)-1)) {
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(conn_pre)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- LC6_wt_biR[[j]][k,2] + LC6_wt_biL[[j]][k,2]
    grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
  }
  breaks <- seq(0,180,length.out = 9)
  grid_Gaussian$equalSpace <- cut(grid_Gaussian$Z, breaks) #cut into 18 levels
  grid_Gaussian <- cbind(grid_Gaussian, rep(j,dim(grid_Gaussian)[1]))
  colnames(grid_Gaussian) <- c("X","Y","Z","equalSpace","gp")
  
  cutoffZ <- max(grid_Gaussian$Z) 
  cutofflevel <- 0.9
  grid_Gaussian <- grid_Gaussian[grid_Gaussian$Z > cutoffZ*cutofflevel, ]
  grid_Gaussian_sum <- rbind(grid_Gaussian_sum, grid_Gaussian[,c("X","Y","Z","gp")])
  grid_Gaussian_ls[[j]] <- grid_Gaussian[,c("X","Y","Z","gp")]
  gprange[[j]] <- range(grid_Gaussian_ls[[j]]$Z)
}
grid_Gaussian_sum <- grid_Gaussian_sum[-1,]
grid_Gaussian_sum$gp <- factor(grid_Gaussian_sum$gp)


xy_bd_chull_df <- as.data.frame(xy_bd_chull)
gg_cont <- ggplot(grid_Gaussian_sum, aes(X, Y, z = Z)) +
  geom_raster(aes(fill = gp)) +
  scale_fill_manual(values = dolphin_col) +
  coord_fixed(ratio = 1) +
  scale_y_reverse()

gg_cont <- gg_cont + 
  labs(title = paste("bi,", cutofflevel*100, "% cutoff"))+
  geom_polygon(data = xy_bd_chull_df, aes(x=phi_deg, y=theta_deg, z=0),colour = "black", alpha = 0) +
  theme_void() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
  annotate("text", x = -5, y = 185, label = "9°") +
  geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
  annotate("text", x = -15, y = 175.5, label = "9°") +
  geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
  geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
  geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
  annotate("text", x = 5, y = 175.5, label = "9°") 

# dolphin in colors
pglo <- ggplot(conn_LC6_tar_bi) +
  geom_point(aes(x = x, y = y, colour = gp_x), shape = 16) 

for (j in 1:(length(glo_div)-1)) { #add ahulls
  pglo <- pglo + 
    geom_segment(data = ii_ahull_ls[[j]], aes(x = x1, y = y1, xend = x2, yend = y2))
}
pglo <- pglo +
  # geom_point(aes(x = x, y = y, colour = gp_x), shape = ".") +
  scale_colour_manual(values = dolphin_col) +
  guides(col = guide_legend(title = "groupings")) +
  theme(legend.position = "bottom") +
  theme_void() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  # scale_y_reverse() +
  ylim(205000, 180000)+
  # scale_x_reverse() +
  # xlab("x") + 
  # ylab("y") +
  coord_fixed(ratio = 1) +
  geom_segment(aes(x = 380000, y = 205000, xend = 390000, yend = 205000), size=2, lineend = "round") +
  annotate("text", x = 385000, y = 204000, label = "10 Âµm")+
  labs(title = "Synapses distribution in glomerulus")
for (j in 1:(length(glo_div)-2)) {
  glo_bd <- data.frame(x1 = range_syn[j,2], x2 = range_syn[j,2], y1 = range_syn[j,3], y2 = range_syn[j,4])
  # pglo <- pglo + geom_segment(data = glo_bd, aes(x = x1, y = y1, xend = x2, yend = y2))
}
for (j in 1:(length(glo_div)-1)) {
  pglo <- pglo + 
    # annotate("text", x = com_syn[j,1], y = com_syn[j,2], label = paste(sprintf("%.1f %%", 100*N_syn[j,2]/sum(N_syn[,2]))))
    annotate("text", x = com_syn[j,1], y = com_syn[j,2], label = paste(sprintf("%.1f ", 100*N_syn[j,2]/sum(N_syn[,2]))))
}

# PLOT, glo with 3 groups
windows(record = F, width = 16, height = 10)
ggdraw() +
  draw_plot(pglo + theme(legend.justification = "bottom"), 0, 0.51, 1, 0.5) +
  draw_plot(gg_cont + theme(axis.title = element_text()), 0, 0, 0.5, 0.5) 



# -- ipsi
ID_wt <- LC6_wt_ipsi
plvl <- list()
grid_Gaussian_ls <- list()
grid_Gaussian_sum <- matrix(ncol = 4)
colnames(grid_Gaussian_sum) <- c("X","Y","Z","gp")
plvl_sum <- ggplot()
# dolphin_col <- rainbow(N_gp + 1, alpha = 1)[1:(N_gp-1)]
pal_dolphin <- c()
gprange <- list()
for (j in 1:(length(glo_div)-1)) {
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(conn_pre)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- ID_wt[[j]][k,2]
    grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
  }
  
  breaks <- seq(0,180,length.out = 9)
  grid_Gaussian$equalSpace <- cut(grid_Gaussian$Z, breaks) 
  grid_Gaussian <- cbind(grid_Gaussian, rep(j,dim(grid_Gaussian)[1]))
  colnames(grid_Gaussian) <- c("X","Y","Z","equalSpace","gp")
  
  cutoffZ <- max(grid_Gaussian$Z) # j=1 -- 179, j=4 -- 108
  cutofflevel <- 0.9
  grid_Gaussian <- grid_Gaussian[grid_Gaussian$Z > cutoffZ*cutofflevel, ]
  grid_Gaussian_sum <- rbind(grid_Gaussian_sum, grid_Gaussian[,c("X","Y","Z","gp")])
  grid_Gaussian_ls[[j]] <- grid_Gaussian[,c("X","Y","Z","gp")]
  gprange[[j]] <- range(grid_Gaussian_ls[[j]]$Z)
}
grid_Gaussian_sum <- grid_Gaussian_sum[-1,]
grid_Gaussian_sum$gp <- factor(grid_Gaussian_sum$gp)


# color contours
xy_bd_chull_df <- as.data.frame(xy_bd_chull)
gg_cont <- ggplot(grid_Gaussian_sum, aes(X, Y, z = Z)) +
  geom_raster(aes(fill = gp)) +
  scale_fill_manual(values = dolphin_col) +
  coord_fixed(ratio = 1) +
  scale_y_reverse()
# for (j in 1:(length(glo_div)-1)) {
#   gg_cont <- gg_cont + 
#     geom_contour(data = grid_Gaussian_ls[[j]], aes(X,Y,z=Z), breaks = seq(gprange[[j]][1], gprange[[j]][2], length.out = 3), color = dolphin_col[j], alpha = 0.9, lwd = 3)
# }
titletext <- paste("ipsi,", cutofflevel*100, "% cutoff")
gg_cont <- gg_cont + 
  labs(title = titletext)+
  geom_polygon(data = xy_bd_chull_df, aes(x=phi_deg, y=theta_deg, z=0),colour = "black", alpha = 0) +
  theme_void() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
  annotate("text", x = -5, y = 185, label = "9°") +
  geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
  annotate("text", x = -15, y = 175.5, label = "9°") +
  geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
  geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
  geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
  annotate("text", x = 5, y = 175.5, label = "9°") 

# dolphin in colors
pglo <- ggplot(conn_LC6_tar_ipsi) +
  geom_point(aes(x = x, y = y, colour = gp_x), shape = 16) 
# geom_segment(data = ii_ahull, aes(x = x1, y = y1, xend = x2, yend = y2)) +
for (j in 1:(length(glo_div)-1)) { #add ahulls
  pglo <- pglo + 
    geom_segment(data = ii_ahull_ls[[j]], aes(x = x1, y = y1, xend = x2, yend = y2))
}
pglo <- pglo +
  scale_colour_manual(values = dolphin_col) +
  guides(col = guide_legend(title = "groupings", override.aes = list(shape = rep(".", 10)))) +
  theme_void() +
  theme(legend.position = "bottom") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  # scale_y_reverse() +
  ylim(205000, 180000)+
  coord_fixed(ratio = 1) +
  geom_segment(aes(x = 380000, y = 205000, xend = 390000, yend = 205000), size=2, lineend = "round") +
  annotate("text", x = 385000, y = 204000, label = "10 Âµm")+
  labs(title = "Synapses distribution in glomerulus")
for (j in 1:(length(glo_div)-2)) {
  glo_bd <- data.frame(x1 = range_syn[j,2], x2 = range_syn[j,2], y1 = range_syn[j,3], y2 = range_syn[j,4])
}
for (j in 1:(length(glo_div)-1)) {
  pglo <- pglo + 
    annotate("text", x = com_syn[j,1], y = com_syn[j,2], label = paste(sprintf("%.1f ", 100*N_syn[j,3]/sum(N_syn[,3]))))
}

# PLOT, glo with 3 groups
windows(record = F, width = 16, height = 10)
ggdraw() +
  draw_plot(pglo + theme(legend.justification = "bottom"), 0, 0.51, 1, 0.5) +
  draw_plot(gg_cont + theme(axis.title = element_text()), 0, 0, 0.5, 0.5) 




# Figure 7A, targets with volume ref --------------------------------------------------------------------------

nopen3d()
par3d('windowRect' = c(100,100,1100,1100))
shade3d(v14, alpha=0.05)
shade3d(glo_vol, col= "gold", alpha = 0.5)
rgl.viewpoint(userMatrix = rotationMatrix(1*90/180*pi+pi/2,1,0,0) %*% rotationMatrix(0,0,1,0), zoom = 0.7)
tar <- neu_biL[[1]]
plot3d(tar, lwd = 2, col = "blue")
points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = "blue", size = 20)
tar <- neu_biR[[1]]
plot3d(tar, lwd = 2, col = "red")
points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = "red", size = 20)
tar <- neu_ipsi[[1]]
plot3d(tar, lwd = 2, col = "cyan")
points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = "cyan", size = 20)

# rgl.snapshot(filename = "LC6_targest.png",fmt = "png")
# rgl.snapshot(filename = "glo.png",fmt = "png")


# # -- cp LM
# nopen3d()
# par3d('windowRect' = c(100,100,1100,1100))
# # shade3d(v14, alpha=0.05)
# shade3d(glo_vol, col = 'grey90', alpha = 0.3)
# ipsi_pal <- brewer.pal(6, "RdYlBu")[c(1,2,3,5,6)]
# tar <- neu_ipsi[[2]]
# plot3d(tar, col = ipsi_pal[1])
# points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = ipsi_pal[1], size = 20)
# tar <- neu_ipsi[[3]]
# plot3d(tar, col = ipsi_pal[5])
# points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = ipsi_pal[5], size = 20)
# rgl.viewpoint(userMatrix = rotationMatrix(1*90/180*pi+pi/2,1,0,0) %*% rotationMatrix(0,0,1,0), zoom = 1)


#  Figure S7, individual targets  -------------------------------------------------------------------------------
target_names <- c(paste("BiL_", biL_skid, sep = ""), 
                  paste("BiR_", biR_skid, sep = ""), 
                  paste("Ipsi_", ipsi_skid, sep = ""))
pal_target <- pnw_palette("Bay",9)

for (j in 1:9) {
  nopen3d()
  par3d('windowRect' = c(100,100,1100,1100))
  rgl.viewpoint(userMatrix = rotationMatrix(1*90/180*pi+pi/2,1,0,0) %*% rotationMatrix(0,0,1,0), zoom = 1)
  # shade3d(v14, alpha=0.05)
  shade3d(glo_vol, col = 'grey90', alpha = 0.3)
  tar <- neu_target[[j]]
  plot3d(tar, col = pal_target[j], lwd = 3)
  points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = pal_target[j], size = 20)
  title3d(target_names[j])
  rgl.snapshot(filename = paste(target_names[j], ".png", sep = ''),fmt = "png")
}



#  Figure 7D -- time series mean + EM contour  --------------------------------------------------------------------

# -- put in contour's coord
# exp data range [-0.1, 1.4] vs [1, 280]/200
# need to reverse y-axis

# # starting locations
# tseries_x0 <- sort(unique(expBi_df[,1])) - 4.5
# tseries_y0 <- sort(unique(expBi_df[,2])) + 4.5


# -- add up the same type
mat_names <- c("bi", "ipsi")
conn_target_agglo_sum <- list()
conn_target_agglo_sum[[1]] <- conn_target[[1]] + conn_target[[2]] + conn_target[[3]] + conn_target[[4]]
conn_target_agglo_sum[[2]] <- conn_target[[5]]+conn_target[[6]]+conn_target[[7]]+conn_target[[8]]+conn_target[[9]]


exp_raw <- list(expBi3, expIpsi3)
exp_raw_mean <- list(expBi, expIpsi)
stim_end <- 250 # loom during the first 150 points, 4 sec
stim_start <- 0

# for exp data contour 
# x4 <- seq(0-4.25,117-4.25,by = 9)
# y4 <- seq(54, 108, by = 1)
# xygrid4 <- expand.grid(x = x4, y = y4)
x4 <- seq(-24, 111,by = 1) # range(loom_phi)
y4 <- seq(58, 119, by = 1) # range(loom_theta)
xygrid4 <- expand.grid(x = x4, y = y4)


ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])
plvl <- list()
simdata <- list()
simdata_df <- list()
# LC6_tar_median <- list()
fdr <- list()
ct_cut_pt <- list()

for (j in 1:length(conn_target_agglo_sum)) { 
  for (k in 1:dim(exp_raw[[j]])[3]) {
    amp_max <- quantile(na.omit(c(exp_raw[[j]][,,k])), probs = c(0.98)) # normalize to 98% 
    exp_raw[[j]][,,k] <- exp_raw[[j]][,,k] / amp_max
  }
  mai_mean <- rowMeans(exp_raw[[j]], dims = 2, na.rm = T)
  amp_max <- quantile(na.omit(c(mai_mean)), probs = c(0.98)) # normalize to 98% 
  mai_mean <- mai_mean / amp_max
  
  max(na.omit(c(mai_mean))) / amp_max
  min(na.omit(c(mai_mean))) / amp_max
  
  # t-test
  # mu_test <- sum(exp_raw[[j]][ind_mai[c(12,13,14,27,28,42,56)], stim_start:stim_end, ]) / 7 / dim(exp_raw[[j]])[3]
  mu_test <- sum(exp_raw[[j]][ind_mai[c(13,14,28,42)], stim_start:stim_end, ]) / 4 / dim(exp_raw[[j]])[3]
  pval_tmp <- c()
  for (k in 1:dim(exp_raw[[j]])[1]) {
    pval_tmp[k] <- t.test(colSums(exp_raw[[j]][k, stim_start:stim_end,]), mu = mu_test)$p.value
  }
  fdr[[j]] <- p.adjust(pval_tmp, method = 'BH') # fdr across flies CHECK
  
  exp_df <- data.frame(xygrid2, as.vector(t(exp_raw_mean[[j]])))
  colnames(exp_df) <- c("x","y","z")
  exp_loess <- loess (z ~ x * y, exp_df, degree = 2, span = 0.3)
  loess_fit <- predict (exp_loess, xygrid4, se = T)
  exp_interp <- cbind(xygrid4, melt(loess_fit$fit)$value)
  colnames(exp_interp) <- c('x','y','z')
  # exp_interp <- exp_interp[!is.na(exp_interp[,'z']), ] #why need this?
  
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  # ii <- grid_Gaussian[,1] > 16
  # grid_Gaussian[ii,2] <- slope_ex*(grid_Gaussian[ii,1] - 16) + grid_Gaussian[ii,2]
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(neu)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- conn_target_agglo_sum[[j]][k,"tofrom_glo"]
    grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
    # }
  }
  # grid_Gaussian[ii,2] <- round(- slope_ex*(grid_Gaussian[ii,1] - 16) + grid_Gaussian[ii,2])
  simdata[[j]] <- grid_Gaussian
  simdata_df[[j]] <- simdata[[j]]
  colnames(simdata_df[[j]]) <- c("x","y","z")
  # simdata_df[[j]]$z <- simdata_df[[j]]$z / mean(head(sort(simdata_df[[j]]$z, decreasing = T)))
  simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  plvl[[j]] <- ggplot()
  for (k in 1:dim(mai_mean)[1]) {
    df <- data.frame(x = seq(1,dim(mai_mean)[2])/200*6 + loom_phi[k]-4.5, y = mai_mean[ind_mai[k],]*(-6) + loom_theta[k]+4.5)
    if (fdr[[j]][ind_mai[k]] >= 0.05 ) {
      plvl[[j]] <- plvl[[j]] + geom_line(df, mapping = aes(x,y), colour = 'black', alpha = 0.3, lwd = 1)
    } else {
      plvl[[j]] <- plvl[[j]] + geom_line(df, mapping = aes(x,y), colour = 'black', lwd = 1)
    }
  }
  
  plvl[[j]] <- plvl[[j]] + labs(title = paste("target", mat_names[j], ", 70%", "N=",sum(conn_target_agglo_sum[[j]][,"tofrom_glo"]),sep = " ")) +
    geom_polygon(data = shim_xy, aes(x,y), fill = 'black', alpha = 0.3) +
    geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2) +
    geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61), color = "red", alpha = 1, lwd = 2) +
    geom_segment(aes(x = -12, y = 135, xend = -3, yend = 135), size=2, lineend = "round") +
    annotate("text", x = -3, y = 139, label = "9°") +
    geom_segment(aes(x = -12, y = 135, xend = -12, yend = 124), size=2, lineend = "round") +
    annotate("text", x = -17.5, y = 124, label = "9°") +
    # geom_segment(aes(x = -9, y = 90-lefty, xend = 16, yend = 90), size=2, lineend = "round") +
    # geom_segment(aes(x = 16, y = 90, xend = 117, yend = 90-righty), size=2, lineend = "round") +
    geom_segment(aes(x = 0, y = 50, xend = 0, yend = 120), size = 2) + # range(loom_phi)
    geom_segment(aes(x = 90, y = 50, xend = 90, yend = 120), size = 2) +
    geom_segment(aes(x = -25, y = 90, xend = 115, yend = 90), size = 2) +
    # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
    scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
    scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
    theme_void() +
    coord_fixed(ratio = 1)
  
  # area
  ct_cut_pt[[j]] <- simdata_df[[j]][simdata_df[[j]]$z > 0.7, c('x','y')]
}

windows(record = F, width = 8, height = 8)
plvl[[1]] 
windows(record = F, width = 8, height = 8)
plvl[[2]] 


# -- difference in EM contour
simdata_df_diff <- simdata_df[[1]][,1:3]
simdata_df_diff$z <- simdata_df_diff$z - simdata_df[[2]]$z
simdata_df_diff$z[simdata_df_diff$z < 0] <- simdata_df_diff$z[simdata_df_diff$z < 0] / abs(min(simdata_df_diff$z))
simdata_df_diff$z[simdata_df_diff$z > 0] <- simdata_df_diff$z[simdata_df_diff$z > 0] / max(simdata_df_diff$z)
simdata_df_diff_ib <- simdata_df_diff
simdata_df_diff_ib$z <- - simdata_df_diff_ib$z
simdata_df_diff$equalSpace <- cut(simdata_df_diff$z, seq(-1,1,length.out = n_lvl))
simdata_df_diff_ib$equalSpace <- cut(simdata_df_diff_ib$z, seq(-1,1,length.out = n_lvl))


# individual t series ---------------------------------------------------------------------------------------------

# Figure S5C

# color
# getPalette <- colorRampPalette(brewer.pal(9, "RdYlBu"))
# ct_col <- getPalette(7)
ct_col <- brewer.pal(9, "RdYlBu")[c(1,2,3,4,7,8,9)]

# -- add up the same type
mat_names <- c("bi_all", "ipsi_all")
conn_target_agglo_sum <- list()
conn_target_agglo_sum[[1]] <- conn_target[[1]] + conn_target[[2]] + conn_target[[3]] + conn_target[[4]]
conn_target_agglo_sum[[2]] <- conn_target[[5]]+conn_target[[6]]+conn_target[[7]]+conn_target[[8]]+conn_target[[9]]

exp_raw <- list(expBi3, expIpsi3)
exp_raw_mean <- list(expBi, expIpsi)
stim_end <- 250 # loom during the first 150 points, 4 sec
stim_start <- 0


plvl_ol <- list() # all exp contour
plvl_ol_em <- list() # all exp contour with EM bi - ipsi heat map
plvl_bi <- list()
plvl_ipsi <- list()
plvl_bi_ct <- list()
plvl_ipsi_ct <- list()
bi_ct <- list()
ipsi_ct <- list()
indi_max <- list() #peak value for indivial fly
fdr <- list() # False discovery rate

for (j in 1:length(conn_target_agglo_sum)) { 
  # -- fly mean
  max_tmp <- c()
  for (k in 1:dim(exp_raw[[j]])[3]) {
    amp_max <- quantile(na.omit(c(exp_raw[[j]][,,k])), probs = c(0.98)) # normalize to 98% 
    exp_raw[[j]][,,k] <- exp_raw[[j]][,,k] / amp_max
    max_tmp <- c(max_tmp, round(amp_max,1))
  }
  indi_max[[j]] <- max_tmp
  mai_mean <- rowMeans(exp_raw[[j]], dims = 2, na.rm = T)
  # amp_max <- quantile(na.omit(c(mai_mean)), probs = c(0.98)) # normalize to 98% 
  # mai_mean <- mai_mean / amp_max
  
  # t-test
  # mu_test <- sum(exp_raw[[j]][ind_mai[c(12,13,14,27,28,42,56)], stim_start:stim_end, ]) / 7 / dim(exp_raw[[j]])[3]
  mu_test <- sum(exp_raw[[j]][ind_mai[c(13,14,28,42)], stim_start:stim_end, ]) / 4 / dim(exp_raw[[j]])[3]
  pval_tmp <- c()
  for (k in 1:dim(exp_raw[[j]])[1]) {
    pval_tmp[k] <- t.test(colSums(exp_raw[[j]][k, stim_start:stim_end,]), mu = mu_test)$p.value
  }
  fdr[[j]] <- p.adjust(pval_tmp, method = 'BH')
  
  exp_df <- data.frame(xygrid2, as.vector(t(exp_raw_mean[[j]])))
  colnames(exp_df) <- c("x","y","z")
  exp_loess <- loess (z ~ x * y, exp_df, degree = 2, span = 0.3)
  loess_fit <- predict (exp_loess, xygrid4, se = T)
  exp_interp <- cbind(xygrid4, melt(loess_fit$fit)$value)
  colnames(exp_interp) <- c('x','y','z')
  exp_interp <- exp_interp[!is.na(exp_interp[,'z']), ]
  
  plvl_ol[[j]] <- ggplot()
  plvl_ol_em[[j]] <- ggplot() +
    # geom_polygon(data = shim_xy, aes(x,y), fill = 'black', alpha = 0.3) +
    geom_raster(data = simdata_df_diff, aes(x, y, z = z, fill = equalSpace), interpolate = F) +
    scale_fill_manual(values = pal_tar) +
    geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE)
  
  for (k in 1:dim(mai_mean)[1]) {
    df <- data.frame(x = seq(1,dim(mai_mean)[2])/200*6 + loom_phi[k]-4.5, y = mai_mean[ind_mai[k],]*(-6) + loom_theta[k]+4.5)
    if (fdr[[j]][ind_mai[k]] >= 0.05 ) {
      plvl_ol[[j]] <- plvl_ol[[j]] + geom_line(df, mapping = aes(x,y), colour = 'black', alpha = 0.3, lwd = 1)
      # plvl_ol_em[[j]] <- plvl_ol_em[[j]] + geom_line(df, mapping = aes(x,y), colour = 'black', alpha = 0.3, lwd = 1)
    } else {
      plvl_ol[[j]] <- plvl_ol[[j]] + geom_line(df, mapping = aes(x,y), colour = 'black', lwd = 1)
      # plvl_ol_em[[j]] <- plvl_ol_em[[j]] + geom_line(df, mapping = aes(x,y), colour = 'black', lwd = 1)
    }
  }
  plvl_ol[[j]] <- plvl_ol[[j]] + labs(title = paste("Target ipsi exp vs", mat_names[j], ", 70%", "N=",sum(conn_target_agglo_sum[[j]][,"tofrom_glo"]),sep = " ")) +
    geom_polygon(data = shim_xy, aes(x,y), fill = 'black', alpha = 0.3) +
    # geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2) +
    # geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61), color = "red", alpha = 1, lwd = 2) +
    # geom_segment(aes(x = -12, y = 115, xend = -3, yend = 115), size=2, lineend = "round") +
    # annotate("text", x = -3, y = 119, label = "9°") +
    # geom_segment(aes(x = -12, y = 115, xend = -12, yend = 106), size=2, lineend = "round") +
    # annotate("text", x = -17.5, y = 105, label = "9°") +
    geom_segment(aes(x = -12, y = 135, xend = -3, yend = 135), size=2, lineend = "round") +
    annotate("text", x = -3, y = 139, label = "9°") +
    geom_segment(aes(x = -12, y = 135, xend = -12, yend = 124), size=2, lineend = "round") +
    annotate("text", x = -17.5, y = 124, label = "9°") +
    geom_segment(aes(x = 0, y = 50, xend = 0, yend = 120), size = 2) + # range(loom_phi)
    geom_segment(aes(x = 90, y = 50, xend = 90, yend = 120), size = 2) +
    geom_segment(aes(x = -25, y = 90, xend = 115, yend = 90), size = 2) +
    # geom_segment(aes(x = 0, y = y2[1]-4.5, xend = 0, yend = tail(y2,1)+4.5), size = 1, lineend = "round") +
    # geom_segment(aes(x = 90, y = y2[1]-4.5, xend = 90, yend = tail(y2,1)+4.5), size = 1, lineend = "round") +
    # geom_segment(aes(x = x2[1], y = 90-lefty, xend = 0, yend = 90), size=1, lineend = "round") +
    # geom_segment(aes(x = 0, y = 90, xend = tail(x2,1), yend = 90-righty), size=1, lineend = "round") +
    # geom_segment(aes(x = min(loom_phi), y = 90, xend = max(loom_phi), yend = 90), size=1, lineend = "round") +
    scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
    scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
    theme_void() +
    coord_fixed(ratio = 1)
  
  # for (j in 1:length(conn_target_agglo_sum)) { 
  
  # -- fly indiv
  for (k in 1:dim(exp_raw[[j]])[3]) {
    exp_df <- data.frame(xygrid2, rowSums(exp_raw[[j]][ind_mai, stim_start:stim_end, k]))
    colnames(exp_df) <- c("x","y","z")
    exp_loess <- loess (z ~ x * y, exp_df, degree = 2, span = 0.3)
    loess_fit <- predict (exp_loess, xygrid4, se = T)
    exp_interp <- cbind(xygrid4, melt(loess_fit$fit)$value)
    colnames(exp_interp) <- c('x','y','z')
    exp_interp <- exp_interp[!is.na(exp_interp[,'z']), ]
    
    if (j == 1) {
      plvl_bi[[k]] <- ggplot()
      for (m in 1:dim(exp_raw[[j]])[1]) {
        df <- data.frame(x = seq(1,dim(exp_raw[[j]])[2])/200*6 + loom_phi[m]-4.5, y = exp_raw[[j]][ind_mai[m],,k]*(-6) + loom_theta[m]+4.5)
        plvl_bi[[k]] <- plvl_bi[[k]] + geom_line(df, mapping = aes(x,y), colour = 'black', lwd = 1)
      }
      bi_ct[[k]] <- exp_interp
      
      plvl_bi[[k]] <- plvl_bi[[k]] + labs(title = paste(k, ", bi, 60%" ,sep = " ")) +
        geom_polygon(data = shim_xy, aes(x,y), fill = 'black', alpha = 0.3) +
        # geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2) +
        geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[2*(k-1)+1], alpha = 1, lwd = 2) +
        # geom_segment(aes(x = -12, y = 115, xend = -3, yend = 115), size=2, lineend = "round") +
        # annotate("text", x = -3, y = 119, label = "9°") +
        # geom_segment(aes(x = -12, y = 115, xend = -12, yend = 106), size=2, lineend = "round") +
        # annotate("text", x = -17.5, y = 105, label = "9°") +
        geom_segment(aes(x = loom_phi[1]-4.5, y = loom_theta[1]+4.5, xend = loom_phi[1]-4.5, yend = loom_theta[1]-9+4.5), size=2, lineend = "round") +
        annotate("text", x = loom_phi[1]-4.5-4.5, y = loom_theta[1]-4.5+4.5, label = "ΔF/F") +
        annotate("text", x = loom_phi[1]-4.5-4.5, y = loom_theta[1]+4.5, label = paste(indi_max[[j]][k])) +
        geom_segment(aes(x = loom_phi[1]-4.5, y = loom_theta[1]+4.5, xend = loom_phi[1]+9-4.5, yend = loom_theta[1]+4.5), size=2, lineend = "round") +
        annotate("text", x = loom_phi[1]+4.5-4.5, y = loom_theta[1]+4+4.5, label = "5s") +
        geom_segment(aes(x = 0, y = 50, xend = 0, yend = 120), size = 2) + # range(loom_phi)
        geom_segment(aes(x = 90, y = 50, xend = 90, yend = 120), size = 2) +
        geom_segment(aes(x = -25, y = 90, xend = 115, yend = 90), size = 2) +
        # geom_segment(aes(x = 0, y = y2[1]-4.5, xend = 0, yend = tail(y2,1)+4.5), size = 1, lineend = "round") +
        # geom_segment(aes(x = 90, y = y2[1]-4.5, xend = 90, yend = tail(y2,1)+4.5), size = 1, lineend = "round") +
        # geom_segment(aes(x = x2[1], y = 90-lefty, xend = 0, yend = 90), size=1, lineend = "round") +
        # geom_segment(aes(x = 0, y = 90, xend = tail(x2,1), yend = 90-righty), size=1, lineend = "round") +
        # geom_segment(aes(x = min(loom_phi), y = 90, xend = max(loom_phi), yend = 90), size=1, lineend = "round") +
        scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
        scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
        theme_void() +
        coord_fixed(ratio = 1)
      
      # if (k == 1) {
      #   plvl_bi[[k]] <- plvl_bi[[k]] + 
      #     geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1], xend = loom_phi_mat[1,14], yend = loom_theta_mat[1,14]), size = 1, lineend = "round", col = 'blue') + # exp
      #     geom_segment(aes(x = loom_phi_mat[7,1], y = loom_theta_mat[7,1], xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]), size = 1, lineend = "round", col = 'blue') +
      #     geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1], xend = loom_phi_mat[7,1], yend = loom_theta_mat[7,1]), size = 1, lineend = "round", col = 'blue') + 
      #     geom_segment(aes(x = loom_phi_mat[1,14], y = loom_theta_mat[1,14], xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]), size = 1, lineend = "round", col = 'blue')  
      # }
      
      # -- 10 contours for ephys
      exp_interp_ct <- exp_interp
      exp_interp_ct$z <- exp_interp_ct$z / mean(head(sort(exp_interp_ct$z, decreasing = T)))
      exp_interp_ct$equalSpace <- cut(exp_interp_ct$z, seq(0,max(exp_interp_ct$z),length.out = n_lvl))
      
      plvl_bi_ct[[k]] <-
        ggplot(exp_interp_ct, aes(x, y, z = z)) +
        geom_raster(aes(fill = equalSpace), interpolate = F) +
        scale_fill_manual(values = pal_tar) +
        # geom_contour(data = exp_interp_ct, aes(x,y,z=z), breaks = c(0.7), color = "blue", alpha = 1, lwd = 1) +
        geom_contour(data = exp_interp_ct, aes(x,y,z=z), breaks = c(0.61), color = ct_col[2*(k-1)+1], alpha = 1, lwd = 2) +
        geom_segment(aes(x = 0, y = 50, xend = 0, yend = 120), size = 2) + # range(loom_phi)
        geom_segment(aes(x = 90, y = 50, xend = 90, yend = 120), size = 2) +
        geom_segment(aes(x = -25, y = 90, xend = 115, yend = 90), size = 2) +
        scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
        scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
        theme_void() +
        labs(title = paste(k, ", bi, 60%" ,sep = " ")) +
        coord_fixed(ratio = 1)
    } else {
      plvl_ipsi[[k]] <- ggplot()
      for (m in 1:dim(exp_raw[[j]])[1]) {
        df <- data.frame(x = seq(1,dim(exp_raw[[j]])[2])/200*6 + loom_phi[m]-4.5, y = exp_raw[[j]][ind_mai[m],,k]*(-6) + loom_theta[m]+4.5)
        plvl_ipsi[[k]] <- plvl_ipsi[[k]] + geom_line(df, mapping = aes(x,y), colour = 'black', lwd = 1)
      }
      ipsi_ct[[k]] <- exp_interp
      
      plvl_ipsi[[k]] <- plvl_ipsi[[k]] + labs(title = paste(k, ", ipsi, 60%" ,sep = " ")) +
        geom_polygon(data = shim_xy, aes(x,y), fill = 'black', alpha = 0.3) +
        # geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = ct_col[k], alpha = 0, lwd = 2) +
        geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[k], alpha = 1, lwd = 2) +
        geom_segment(aes(x = loom_phi[1]-4.5, y = loom_theta[1]+4.5, xend = loom_phi[1]-4.5, yend = loom_theta[1]-9+4.5), size=2, lineend = "round") +
        annotate("text", x = loom_phi[1]-4.5-4.5, y = loom_theta[1]-4.5+4.5, label = "ΔF/F") +
        annotate("text", x = loom_phi[1]-4.5-4.5, y = loom_theta[1]+4.5, label = paste(indi_max[[j]][k])) +
        geom_segment(aes(x = loom_phi[1]-4.5, y = loom_theta[1]+4.5, xend = loom_phi[1]+9-4.5, yend = loom_theta[1]+4.5), size=2, lineend = "round") +
        annotate("text", x = loom_phi[1]+4.5-4.5, y = loom_theta[1]+4+4.5, label = "5s") +
        geom_segment(aes(x = 0, y = 50, xend = 0, yend = 120), size = 2) + # range(loom_phi)
        geom_segment(aes(x = 90, y = 50, xend = 90, yend = 120), size = 2) +
        geom_segment(aes(x = -25, y = 90, xend = 115, yend = 90), size = 2) +
        # geom_segment(aes(x = 0, y = y2[1]-4.5, xend = 0, yend = tail(y2,1)+4.5), size = 1, lineend = "round") +
        # geom_segment(aes(x = 90, y = y2[1]-4.5, xend = 90, yend = tail(y2,1)+4.5), size = 1, lineend = "round") +
        # geom_segment(aes(x = x2[1], y = 90-lefty, xend = 0, yend = 90), size=1, lineend = "round") +
        # geom_segment(aes(x = 0, y = 90, xend = tail(x2,1), yend = 90-righty), size=1, lineend = "round") +
        # geom_segment(aes(x = min(loom_phi), y = 90, xend = max(loom_phi), yend = 90), size=1, lineend = "round") +
        scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
        scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
        theme_void() +
        coord_fixed(ratio = 1)
      
      # -- 10 contours for ephys
      exp_interp_ct <- exp_interp
      exp_interp_ct$z <- exp_interp_ct$z / mean(head(sort(exp_interp_ct$z, decreasing = T)))
      exp_interp_ct$equalSpace <- cut(exp_interp_ct$z, seq(0,max(exp_interp_ct$z),length.out = n_lvl))
      
      plvl_ipsi_ct[[k]] <-
        ggplot(exp_interp_ct, aes(x, y, z = z)) +
        geom_raster(aes(fill = equalSpace), interpolate = F) +
        scale_fill_manual(values = pal_tar) +
        geom_contour(data = exp_interp_ct, aes(x,y,z=z), breaks = c(0.61), color = ct_col[k], alpha = 1, lwd = 2) +
        geom_segment(aes(x = 0, y = 50, xend = 0, yend = 120), size = 2) + # range(loom_phi)
        geom_segment(aes(x = 90, y = 50, xend = 90, yend = 120), size = 2) +
        geom_segment(aes(x = -25, y = 90, xend = 115, yend = 90), size = 2) +
        scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
        scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
        theme_void() +
        labs(title = paste(k, ", ipsi, 60%" ,sep = " ")) +
        coord_fixed(ratio = 1)
    }
    if (j == 1) {
      plvl_ol[[j]] <- plvl_ol[[j]] + geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[2*(k-1)+1], alpha = 1, lwd = 2)
      plvl_ol_em[[j]] <- plvl_ol_em[[j]] + geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[2*(k-1)+1], alpha = 1, lwd = 2)
    } else {
      plvl_ol[[j]] <- plvl_ol[[j]] + geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[k], alpha = 1, lwd = 2)
      plvl_ol_em[[j]] <- plvl_ol_em[[j]] + geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[k], alpha = 1, lwd = 2)
    }
  }
}



# -- Figure S3C, ephys contour overlay t-series
windows(record = F, width = 8, height = 8)
plvl_bi[[1]] 
windows(record = F, width = 8, height = 8)
plvl_bi[[2]] 
windows(record = F, width = 8, height = 8)
plvl_bi[[3]] 
windows(record = F, width = 8, height = 8)
plvl_bi[[4]] 

windows(record = F, width = 8, height = 8)
plvl_ipsi[[1]]
windows(record = F, width = 8, height = 8)
plvl_ipsi[[2]]
windows(record = F, width = 8, height = 8)
plvl_ipsi[[3]]
windows(record = F, width = 8, height = 8)
plvl_ipsi[[4]]
windows(record = F, width = 8, height = 8)
plvl_ipsi[[5]]
windows(record = F, width = 8, height = 8)
plvl_ipsi[[6]]
windows(record = F, width = 8, height = 8)
plvl_ipsi[[7]]


# -- Figure S3 , ephys contours
for (k in 1:4) {
  pdf(file = paste("ephys_bi_contour_", k, ".pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
  print(plvl_bi_ct[[k]])
  dev.off()
}
for (k in 1:7) {
  pdf(file = paste("ephys_ipsi_contour_", k, ".pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
  print(plvl_ipsi_ct[[k]] )
  dev.off()
}

# windows(record = F, width = 8, height = 8)
# plvl_bi_ct[[1]] 

# -- Figure 3D alt
windows(record = F, width = 8, height = 8)
plvl_ol[[1]] 
windows(record = F, width = 8, height = 8)
plvl_ol[[2]] 

# -- Figure S7, bi - ipsi, indiv exp
windows(record = F, width = 8, height = 8)
plvl_ol_em[[1]] 
windows(record = F, width = 8, height = 8)
plvl_ol_em[[2]] 

# -- Figure 6D alt
# contour map + exp contour
windows(record = F, width = 8, height = 8)
plvl_ol[[1]] + geom_contour(data = simdata_df[[1]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2)

windows(record = F, width = 8, height = 8)
plvl_ol[[2]] + geom_contour(data = simdata_df[[2]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2)




#  exp and EM overlap ---------------------------------------------------------------------------------------------

# polygon from EM
simdata_ct <- list()
simdata_poly <- list()
for (j in 1:2) {
  simdata_ct[[j]] <- simdata_df[[j]][simdata_df[[j]]$z > 0.71, 1:2]
  ash <- ashape(simdata_ct[[j]]+matrix(runif(dim(simdata_ct[[j]])[1]*2, 1e-9, 2e-9), ncol = 2), alpha = 5)
  simdata_poly[[j]] <- mkpoly(ash$edges)[[1]][,3:4]
  # dev.new()
  # polygon(simdata_poly[[j]])
}


# polygon from exp
exp_raw <- list(expBi3, expIpsi3)
stim_end <- 200 # loom during the first 150 points, 4 sec
stim_start <- 50
ol_ratio <- matrix(ncol = 2, nrow = 2*(4+7)) # percentage of overlap = intersection / union, bi-biEM, ipsi-biEM, etc

NN <- 1
for (jem in 1:2) {
  for (j in 1:2) {
    for (k in 1:dim(exp_raw[[j]])[3]) {
      exp_df <- data.frame(xygrid2, rowSums(exp_raw[[j]][ind_mai, stim_start:stim_end, k]))
      colnames(exp_df) <- c("x","y","z")
      exp_loess <- loess (z ~ x * y, exp_df, degree = 2, span = 0.3)
      loess_fit <- predict (exp_loess, xygrid4, se = T)
      exp_interp <- cbind(xygrid4, melt(loess_fit$fit)$value)
      colnames(exp_interp) <- c('x','y','z')
      exp_interp <- exp_interp[!is.na(exp_interp[,'z']), ]
      
      exp_interp_ct <- exp_interp[(exp_interp$z > 0.61*max(exp_interp$z)),]
      
      N_in <- sp::point.in.polygon(exp_interp_ct[,1], exp_interp_ct[,2], simdata_poly[[jem]][,1], simdata_poly[[jem]][,2])
      N_in <- sum(N_in)
      N_out <- dim(exp_interp_ct)[1] - N_in
      N_em <- sum(simdata_ct[[jem]]$y >= 54 | simdata_ct[[jem]]$y <= 108)
      ol_ratio[NN,] <- c(j + 2*(jem %/% 2), N_in / (N_out + N_em))
      NN <- NN + 1
    }
  }
}

ol_ratio <- as.data.frame(ol_ratio)
colnames(ol_ratio) <- c('pair', 'ratio')


windows(record = F, width = 8, height = 8)
ggplot(ol_ratio, aes(x = pair, group = pair, y = ratio)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 1) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4), labels=c("bi-biEM", "ipsi-biEM", "bi-ipsiEM", "ipsi-ipsiEM"))+
  theme_bw()

# ggplot(ol_ratio, aes(x = pair, group = pair, y = ratio)) + geom_boxplot()+
#   scale_x_continuous(breaks=c(1,2,3,4), labels=c("bi-biEM", "ipsi-biEM", "bi-ipsiEM", "ipsi-ipsiEM")) +
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) 




# delta  EM  ------------------------------------------------------------------------------------------------------

simdata_df_diff <- simdata_df[[1]][,1:3]
simdata_df_diff$z <- simdata_df_diff$z - simdata_df[[2]]$z
simdata_df_diff$z <- simdata_df_diff$z / max(abs(simdata_df_diff$z))
# simdata_df_diff$equalSpace <- cut(simdata_df_diff$z, seq(0,max(abs(simdata_df_diff$z)),length.out = n_lvl))
# simdata_df_diff$equalSpace <- cut(simdata_df_diff$z, seq(min(simdata_df_diff$z),max(simdata_df_diff$z),length.out = n_lvl))
simdata_df_diff$equalSpace <- cut(simdata_df_diff$z, seq(-1,1,length.out = n_lvl))

# contour map
windows(record = F, width = 8, height = 8)
ggplot(simdata_df_diff, aes(x, y, z = z)) +
  geom_raster(aes(fill = equalSpace), interpolate = F) +
  scale_fill_manual(values = pal_tar) +
  geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE) +
  geom_segment(aes(x = 10, y = 180, xend = 19, yend = 180), size=2, lineend = "round") +
  annotate("text", x = 14.5, y = 185, label = "9°") +
  geom_segment(aes(x = 10, y = 180, xend = 10, yend = 171), size=2, lineend = "round") +
  annotate("text", x = 5, y = 175.5, label = "9°") +
  # geom_segment(aes(x = -9, y = 90-lefty, xend = 16, yend = 90), size=2, lineend = "round") +
  # geom_segment(aes(x = 16, y = 90, xend = 180, yend = 79.1), size=2, lineend = "round") + #90-slope_ex*(180-16)
  geom_segment(aes(x = min(loom_phi), y = 90, xend = 160, yend = 90), size=1, lineend = "round") +
  # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 180), size = 2, lineend = "round") +
  geom_segment(aes(x = 90, y = 0, xend = 90, yend = 180), size = 2, lineend = "round") +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  theme_void() +
  geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1], xend = loom_phi_mat[1,14], yend = loom_theta_mat[1,14]), size = 1, lineend = "round", col = 'red') + # exp
  geom_segment(aes(x = loom_phi_mat[7,1], y = loom_theta_mat[7,1], xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]), size = 1, lineend = "round", col = 'red') +
  geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1], xend = loom_phi_mat[7,1], yend = loom_theta_mat[7,1]), size = 1, lineend = "round", col = 'red') +
  geom_segment(aes(x = loom_phi_mat[1,14], y = loom_theta_mat[1,14], xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]), size = 1, lineend = "round", col = 'red')
labs(title = "bi - ipsi") +
  coord_fixed(ratio = 1)

# Figure delta EM
# bi - ipsi
windows(record = F, width = 8, height = 8)
ggplot(simdata_df_diff, aes(x, y, z = z)) +
  geom_raster(aes(fill = equalSpace), interpolate = F) +
  scale_fill_manual(values = pal_tar) +
  guides(fill = guide_legend(reverse=T)) +
  geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE) +
  geom_segment(aes(x = 10, y = 190, xend = 19, yend = 190), size=2, lineend = "round") +
  annotate("text", x = 14.5, y = 195, label = "9°") +
  geom_segment(aes(x = 10, y = 190, xend = 10, yend = 181), size=2, lineend = "round") +
  annotate("text", x = 5, y = 185.5, label = "9°") +
  geom_segment(aes(x = min(loom_phi), y = 90, xend = 160, yend = 90), size=1, lineend = "round") +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 180), size = 1, lineend = "round") +
  geom_segment(aes(x = 90, y = 0, xend = 90, yend = 180), size = 1, lineend = "round") +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  theme_void() +
  labs(title = "bi - ipsi and bi contour") +
  geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1]-4.5, xend = loom_phi_mat[1,14], yend = loom_theta_mat[1,14]-4.5), size = 1, lineend = "round", col = 'red') + # exp
  geom_segment(aes(x = loom_phi_mat[7,1], y = loom_theta_mat[7,1]+4.5, xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]+4.5), size = 1, lineend = "round", col = 'red') +
  geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1]-4.5, xend = loom_phi_mat[7,1], yend = loom_theta_mat[7,1]+4.5), size = 1, lineend = "round", col = 'red') +
  geom_segment(aes(x = loom_phi_mat[1,14], y = loom_theta_mat[1,14]-4.5, xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]+4.5), size = 1, lineend = "round", col = 'red') +
  geom_contour(data = bi_ct[[1]], aes(x,y,z=z), breaks = c(0.61 * max(bi_ct[[1]][,3])), color = ct_col[1], alpha = 1, lwd = 2) +
  geom_contour(data = bi_ct[[2]], aes(x,y,z=z), breaks = c(0.61 * max(bi_ct[[2]][,3])), color = ct_col[3], alpha = 1, lwd = 2) +
  geom_contour(data = bi_ct[[3]], aes(x,y,z=z), breaks = c(0.61 * max(bi_ct[[3]][,3])), color = ct_col[5], alpha = 1, lwd = 2) +
  geom_contour(data = bi_ct[[4]], aes(x,y,z=z), breaks = c(0.61 * max(bi_ct[[4]][,3])), color = ct_col[7], alpha = 1, lwd = 2) +
  coord_fixed(ratio = 1)

# ipsi - bi
windows(record = F, width = 8, height = 8)
# pdf(file = "ipsi_bi.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(simdata_df_diff_ib, aes(x, y, z = z)) +
  geom_raster(aes(fill = equalSpace), interpolate = F) +
  scale_fill_manual(values = pal_tar ) +
  guides(fill = guide_legend(reverse=T)) +
  geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE) +
  geom_segment(aes(x = 10, y = 190, xend = 19, yend = 190), size=2, lineend = "round") +
  annotate("text", x = 14.5, y = 195, label = "9°") +
  geom_segment(aes(x = 10, y = 190, xend = 10, yend = 181), size=2, lineend = "round") +
  annotate("text", x = 5, y = 185.5, label = "9°") +
  geom_segment(aes(x = min(loom_phi), y = 90, xend = 160, yend = 90), size=1, lineend = "round") +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 180), size = 1, lineend = "round") +
  geom_segment(aes(x = 90, y = 0, xend = 90, yend = 180), size = 1, lineend = "round") +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  theme_void() +
  labs(title = "ipsi - bi and ipsi contour") +
  geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1]-4.5, xend = loom_phi_mat[1,14], yend = loom_theta_mat[1,14]-4.5), size = 1, lineend = "round", col = 'red') + # exp
  geom_segment(aes(x = loom_phi_mat[7,1], y = loom_theta_mat[7,1]+4.5, xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]+4.5), size = 1, lineend = "round", col = 'red') +
  geom_segment(aes(x = loom_phi_mat[1,1], y = loom_theta_mat[1,1]-4.5, xend = loom_phi_mat[7,1], yend = loom_theta_mat[7,1]+4.5), size = 1, lineend = "round", col = 'red') +
  geom_segment(aes(x = loom_phi_mat[1,14], y = loom_theta_mat[1,14]-4.5, xend = loom_phi_mat[7,14], yend = loom_theta_mat[7,14]+4.5), size = 1, lineend = "round", col = 'red') +
  geom_contour(data = ipsi_ct[[1]], aes(x,y,z=z), breaks = c(0.61 * max(ipsi_ct[[1]][,3])), color = ct_col[1], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[2]], aes(x,y,z=z), breaks = c(0.61 * max(ipsi_ct[[2]][,3])), color = ct_col[2], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[3]], aes(x,y,z=z), breaks = c(0.61 * max(ipsi_ct[[3]][,3])), color = ct_col[3], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[4]], aes(x,y,z=z), breaks = c(0.61 * max(ipsi_ct[[4]][,3])), color = ct_col[4], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[5]], aes(x,y,z=z), breaks = c(0.61 * max(ipsi_ct[[5]][,3])), color = ct_col[5], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[6]], aes(x,y,z=z), breaks = c(0.61 * max(ipsi_ct[[6]][,3])), color = ct_col[6], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[7]], aes(x,y,z=z), breaks = c(0.61 * max(ipsi_ct[[7]][,3])), color = ct_col[7], alpha = 1, lwd = 2) +
  coord_fixed(ratio = 1)

dev.off()


# cross section occoupancy  delta-EM ---------------------------------------------------------------------------------

# resample glo LC6
LC6_sp_glo_rs <- LC6 # initiate
for (ni in 1:length(LC6)) {
  tar <- LC6[[ni]]
  tar_xyz <- xyzmatrix(tar$d)
  ind <- which.min(rowSums(sweep(tar_xyz, MARGIN = 2, STATS = start_xyz, '-')^2))
  targ <- as.ngraph(tar)
  distal_points <- igraph::graph.dfs(targ, root=ind, unreachable=FALSE, neimode='out')$order
  tar_distal <- subset(tar, distal_points)
  sp <- spine(tar_distal, UseStartPoint = T)
  sp_glo <- subset(sp, pointsinside(sp, surf=glo.msh, rval='distance') > -500)
  LC6_sp_glo_rs[[ni]] <- resample(sp_glo, stepsize = 400)
}

# choose 8 LC, frontal (delta RF), lateral or random
ind_bi_ipsi_front <- c(48,1,4,39,44,14,19,8) #front
ind_bi_ipsi_lat <- c(24,25,56,55,57,9,26,2) #lateral

ind_bi_ipsi <- ind_bi_ipsi_front
# ind_bi_ipsi <- ind_bi_ipsi_lat

# ind_bi_ipsi <- sample(65, 8) # random

# -- ref neu
ind_sep <- c(44,46)

o_xyz_m <- xyzmatrix(LC6_sp_glo_rs[[ind_sep[1]]]$d)
z_xyz_m <- xyzmatrix(LC6_sp_glo_rs[[ind_sep[2]]]$d)

neu_target_rs <- resample(neu_target, stepsize = 400)

# nopen3d()
# plot3d(LC6[[i_AP[1]]], WithNodes = F)
# points3d(xyzmatrix(LC6[[i_AP[1]]]$d), size = 10)
# plot3d(LC6)
# identify3d(xyzmatrix(LC6[[i_AP[1]]]$d))

dL <- 500 # half thickness of cross section
range_x <- c(358000, 438000)
N_cs <- diff(range_x) / 2 / dL

xp <- c(1, 0, 0)

yz_ls <- list()
yz_ls_tar <- list()
dd_ls_tar <- list()
ref_d <- c()

for (j in 1:N_cs) {
  ii <- o_xyz_m[,1] > (range_x[1] + (j-1)*2*dL) & o_xyz_m[,1] < (range_x[1] + (j)*2*dL)
  o_xyz <- colMeans(o_xyz_m[ii,])
  d1 <- (o_xyz + xp*dL) %*% xp #distance to plane
  d2 <- (o_xyz - xp*dL) %*% xp #distance to plane
  
  # use 2 planes
  z_ind <- apply(z_xyz_m, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
  if (sum(z_ind) == 1 ) {
    z_xyz <- z_xyz_m[z_ind,] 
  } else {
    z_xyz <- colMeans(z_xyz_m[z_ind,]) 
  }
  z_xyz[1] <- o_xyz[1] # use this plane as y-z
  zp <- z_xyz - o_xyz
  
  ref_d <- c(ref_d, sqrt(sum(zp^2))) #unit length
  
  pt_yz <- matrix(NA, ncol = 2, nrow = length(LC6))
  for (k in 1:nrow(pt_yz)) {
    xyz_m <- xyzmatrix(LC6_sp_glo_rs[[k]]$d)
    # xyz_m <- xyz_m[xyz_m[,1] > glo_div[j] & xyz_m[,1] < glo_div[j+1],]
    xyz_m_ind <- apply(xyz_m, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
    if (sum(xyz_m_ind) > 0) {
      if (sum(xyz_m_ind) == 1) { #if single pt
        xyz <- xyz_m[xyz_m_ind,] 
      } else {
        xyz <- colMeans(xyz_m[xyz_m_ind,]) 
      }
      v <- xyz - o_xyz
      v <- v - c(v %*% xp)*xp 
      ang <- acos(v %*% zp / sqrt(sum(v^2)) / sqrt(sum(zp^2)))
      ang <- ang * sign(xp %*% cross3D(zp, v) )
      mag <- sqrt(sum(v^2)) / sqrt(sum(zp^2)) #normalize to 1
      pt_yz[k,] <- c(sin(ang)*mag, cos(ang)*mag)
    }
  }
  pt_yz[ind_sep[1],] <- c(0,0)
  pt_yz[ind_sep[2],] <- c(0,1)
  yz_ls[[j]] <- pt_yz
  
  # target
  dd_tar <- matrix(NA, ncol = 1, nrow = 9)
  pt_yz_tar <- matrix(NA, ncol = 2, nrow = 9)
  for (k in 1:nrow(pt_yz_tar)) {
    ii <- pointsinside(xyzmatrix(neu_target_rs[[k]]$d), surf=glo.msh, rval='logical') 
    xyz_m <- xyzmatrix(neu_target_rs[[k]]$d)[ii,]
    xyz_m_ind <- apply(xyz_m, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
    xyz_m_in <- xyz_m[xyz_m_ind,,drop=F]
    if (sum(xyz_m_ind) > 0) {
      if (sum(xyz_m_ind) == 1) { #if single pt
        xyz <- xyz_m[xyz_m_ind,] 
      } else {
        xyz <- colMeans(xyz_m[xyz_m_ind,]) 
      }
      v <- xyz - o_xyz
      v <- v - c(v %*% xp)*xp 
      ang <- acos(v %*% zp / sqrt(sum(v^2)) / sqrt(sum(zp^2)))
      ang <- ang * sign(xp %*% cross3D(zp, v) )
      mag <- sqrt(sum(v^2)) / sqrt(sum(zp^2)) #normalize to 1
      pt_yz_tar[k,] <- c(sin(ang)*mag, cos(ang)*mag)
      
      dd <- c()
      for (itar in 1:nrow( xyz_m_in )) {
        # for (iLC in 1:65) {
        for (iLC in ind_bi_ipsi) {
          xyz_m_LC <- xyzmatrix(LC6_sp_glo_rs[[iLC]]$d)
          xyz_m_ind <- apply(xyz_m_LC, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
          if (sum(xyz_m_ind) > 0) {
            if (sum(xyz_m_ind) == 1) { #if single pt
              xyz_LC <- xyz_m_LC[xyz_m_ind,] 
            } else {
              xyz_LC <- colMeans(xyz_m_LC[xyz_m_ind,]) 
            }
            dd <- c(dd, sqrt(sum((xyz_LC - xyz_m_in[itar,])^2)))
          }
        }
      }
      dd_tar[k] <- mean(dd)
    }
  }
  yz_ls_tar[[j]] <- pt_yz_tar
  dd_ls_tar[[j]] <- dd_tar
}


# dd_ls_tar_front <- dd_ls_tar
# dd_ls_tar_lateral <- dd_ls_tar

# - control, max diameter in each cross section
x <- seq(range_x[1], range_x[2], by = dL)
y <- seq(182500, 204500, by = dL)
z <- seq(131000, 161000, by = dL)
glo_xyzgrid <- expand.grid(x, y, z)
glo_xyzgrid <- subset(glo_xyzgrid, pointsinside(glo_xyzgrid, surf=glo.msh, rval='distance') > -dL/2)
colnames(glo_xyzgrid) <- c('x','y','z')

cs_dia <- matrix(ncol = 1, nrow = N_cs)
for (j in 1:N_cs) {
  # ii <- glo_xyzgrid[,1] < (range_x[2] - (j-1)*2*dL) & glo_xyzgrid[,1] > (range_x[2] - (j)*2*dL)
  ii <- glo_xyzgrid[,1] > (range_x[1] + (j-1)*2*dL) & glo_xyzgrid[,1] < (range_x[1] + (j)*2*dL)
  xyz <- glo_xyzgrid[ii,]
  cs_dia[j] <- max(dist(xyz))
}


# PLOT, Figure S7E
dd_tar <- matrix(unlist(dd_ls_tar), ncol = 9, byrow = T)
# dd_tar <- dd_tar[seq(nrow(dd_tar),1), ] #reverse

# --- choose one
ii_target <- c(3,7) #2
# ii_target <- c(4,5,6,8,9) #rest

df <- data.frame(dd_tar[1:80, ii_target])
df <- cbind(seq(1,80), df)
colnames(df) <- c('pos', ii_target)


df$dia <- cs_dia[1:80]/2
dfm <- melt(df, id = 'pos')

tick_x <- rep("", 81)
tick_x[seq(0,80,by = 20)+1] <- seq(0,80,by = 20)


dev.new()
# pdf(file = paste("avg dist to lat LC rest.pdf",sep = ''), width = 8, height = 4,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(x=pos, y=value) ) + 
  # geom_point( ) +
  geom_line(aes(colour = factor(variable), group = variable), lwd = 3) +
  # coord_cartesian(ylim = c(-.5, 1)) +
  # theme_minimal_grid() +
  # theme_bw() +
  scale_colour_manual(values = c(pal_target[ii_target], 'grey'),
                      labels = c(target_names[ii_target], "radius")) +
  scale_x_continuous(breaks = seq(0,80, by=1), labels = tick_x, expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,8000, by=2000), labels = seq(0,8, by=2)) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color='black'),
        # panel.grid.minor.x = element_line(colour = 'black', size = 0.2),
        axis.ticks.length.x = unit(-.25, "cm"),
        axis.text.x = element_text(margin = margin(t = .3, unit = "cm"))) +
  labs(title = "avg dist to lat LC", x='glomerulus axis [um]', y='avg dist [um]') 

dev.off()


# Figure S7D, glo with 80 cuts
nopen3d()
par3d('windowRect' = c(100,100,1700,1700))
shade3d(glo.msh, col='grey',alpha = 0.5, lit=F)
rgl.viewpoint(userMatrix = rotationMatrix(0/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0))
par3d("FOV" = 0)
segments3d(rbind(c(range_x[1],210000,150000), c(range_x[1]+80*2*dL,210000,150000)), lwd = 5, col = 'black')
segments3d(rbind(c(range_x[1]+80*2*dL,210000,150000),
                 c(range_x[1]+80*2*dL,210000-cos(30/180*pi)*10000,150000-sin(30/180*pi)*10000)), lwd = 5)
segments3d(rbind(c(range_x[1]+20*2*dL,230000-cos(30/180*pi)*12000,130000-sin(30/180*pi)*12000),
                 c(range_x[1]+20*2*dL,230000-cos(30/180*pi)*25000,130000-sin(30/180*pi)*25000)),
           lwd = 5, col = 'black')
# axes3d(c('x','y','z')); title3d('','','x','y','z')

# rgl.snapshot(filename = "80 sectors.png",fmt = "png")


# # PLOT cable
# m <- 30
# 
# df_plt <- rbind(yz_ls[[m]][c(ind_sep,ind_bi_ipsi), ], yz_ls_tar[[m]][c(3,4,7),]) %>% 
#   as.data.frame()
# df_plt$colcol <- factor(c(rep(1,2),
#                           rep(2,length(ind_bi_ipsi)),
#                           rep(3,2),
#                           rep(4,1)))
# colnames(df_plt) <- c('y','z','colcol')
# gpl <- ggplot() +
#   geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 20 ) +
#   scale_colour_manual(values = col4,
#                       breaks = c("1", "2", "3","4"),
#                       labels = c("ref", "front", "bi", "ipsi") ) +
#   guides(colour = guide_legend("groups") )+
#   # scale_y_reverse() +
#   # xlim(xrange) +
#   # ylim(yrange) +
#   coord_fixed(ratio = 1) +
#   theme_minimal() +
#   # theme_void() +
#   # geom_segment(aes(x = 0, y = -0.5, xend = 1/ref_d[m]*1e3, yend = -0.5), size=2, lineend = "round") +
#   # annotate("text", x = 0.05, y = -0.8, label = '1 um') +
#   labs(title = m)
# 
# windows(record = F, width = 8, height = 8)
# # pdf(file = paste(m, ".pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
# gpl
# 
# dev.off()



# Figure 7E, cross section, 

j <- 20

ii <- o_xyz_m[,1] > (range_x[1] + (j-1)*2*dL) & o_xyz_m[,1] < (range_x[1] + (j)*2*dL)
# ii <- o_xyz_m[,1] < (range_x[2] - (j-1)*2*dL) & o_xyz_m[,1] > (range_x[2] - (j)*2*dL)
o_xyz <- colMeans(o_xyz_m[ii,])
d1 <- (o_xyz + xp*dL) %*% xp #distance to plane
d2 <- (o_xyz - xp*dL) %*% xp #distance to plane
# use 2 planes
z_ind <- apply(z_xyz_m, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
z_xyz <- colMeans(z_xyz_m[z_ind,,drop=F]) 
zp <- z_xyz - o_xyz
zp <- zp - c(zp %*% xp)*xp # +z

ref_d <- c(ref_d, sqrt(sum(zp^2))) #unit length

pt_yz <- matrix(NA, ncol = 2, nrow = length(LC6))
for (k in 1:nrow(pt_yz)) {
  xyz_m <- xyzmatrix(LC6_sp_glo_rs[[k]]$d)
  xyz_m_ind_1 <- apply(xyz_m, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
  xyz_m_ind <- xyz_m_ind_1 #& xyz_m_ind_2 # within 2 planes and 20 um within ref pt
  if (sum(xyz_m_ind) > 0) {
    xyz <- colMeans(xyz_m[xyz_m_ind,,drop=F]) 
    v <- xyz - o_xyz
    v <- v - c(v %*% xp)*xp 
    ang <- acos(v %*% zp / sqrt(sum(v^2)) / sqrt(sum(zp^2)))
    ang <- ang * sign(xp %*% cross3D(zp, v) )
    mag <- sqrt(sum(v^2)) / sqrt(sum(zp^2)) #normalize to 1
    pt_yz[k,] <- c(sin(ang)*mag, cos(ang)*mag)
  }
}
pt_yz[ind_sep[1],] <- c(0,0)
pt_yz[ind_sep[2],] <- c(0,1)

it <- 3
ii <- pointsinside(xyzmatrix(neu_target_rs[[it]]$d), surf=glo.msh, rval='logical') 
xyz_m <- xyzmatrix(neu_target_rs[[it]]$d)[ii,]
xyz_m_ind <- apply(xyz_m, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
xyz_m_in <- xyz_m[xyz_m_ind,,drop=F]
# xyz <- colMeans(xyz_m[xyz_m_ind,]) 
pt_yz_k <- matrix(NA, ncol = 2, nrow = nrow(xyz_m_in) )
for (k in 1:nrow(pt_yz_k)) {
  xyz <- xyz_m_in[k,]
  v <- xyz - o_xyz
  v <- v - c(v %*% xp)*xp 
  ang <- acos(v %*% zp / sqrt(sum(v^2)) / sqrt(sum(zp^2)))
  ang <- ang * sign(xp %*% cross3D(zp, v) )
  mag <- sqrt(sum(v^2)) / sqrt(sum(zp^2)) #normalize to 1
  pt_yz_k[k,] <- c(sin(ang)*mag, cos(ang)*mag)
}

df <- rbind(pt_yz, pt_yz_k)
col_ind <- rep(1,65)
col_ind[ind_bi_ipsi_front] <- 2
df <- cbind(df, c(col_ind, rep(3,nrow(pt_yz_k))))
df <- as.data.frame(df)
colnames(df) <- c('x','y','type')

gpl <- ggplot(df) +
  geom_point(aes(x=x, y=y, colour = factor(type), size = type), shape = 16, alpha=0.8 ) +
  scale_colour_manual(values = c('grey', col4_a1[4], pal_target[it]),
                      breaks = c("1", "2", "3"),
                      labels = c("LC6", "LC6 front", neu_target[[it]]$NeuronName) ) +
  scale_size(range = c(4,6)) +
  coord_fixed(ratio = 1) +
  theme_void() +
  geom_segment(aes(x = 0, y = -0.5, xend = 1/ref_d[it]*1e3, yend = -0.5), size=2, lineend = "round") +
  annotate("text", x = 0.05, y = -0.8, label = '1 um') +
  labs(title = "cross section occupancy")

windows(record = F, width = 8, height = 8)
# pdf(file = paste("cross section occ.pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl

dev.off()



# - Figure S7, 8 LC choice

windows(record = F, width = 8, height = 8)

# pdf(file = "2 groups of 8 LO.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
bd_theta <- seq(1, 180, by = 1)
xy_bd <- matrix(ncol = 2)
bd_grid <- expand.grid(bd_phi, bd_theta)
plot(bd_grid, ylim = rev(range(bd_grid$Var2)), type = "n", axes = FALSE, ann = F)
i1i <- 1; i2i <- 1; i3i <- 1; i4i <- 1
for (j in 1:length(xy_poly)) {
  xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
  if (j %in% ind_bi_ipsi_front) {
    # geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 18 )
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4[4], cex = 3, pch = 16)
    # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = i1i, pos = 1, offset = 0.3)
    # i1i <- i1i + 1
  }
  else if (j %in% ind_bi_ipsi_lat) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4[2], cex = 3, pch = 16)
    # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = i3i, pos = 1, offset = 0.3)
    # i3i <- i3i + 1
  }
  else {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="grey", cex = 2, pch = 16) 
  }
}
xy_bd <- xy_bd[-1,]
hpts_2 <- chull(xy_bd)
hpts_2 <- c(hpts_2, hpts_2[1])
xy_bd_chull <- xy_bd[hpts_2,] # hull edge points
polygon(xy_bd_chull)
lines(rbind(c(-11,180), c(-2,180)), lwd = 1)
text(-5, 180, labels = expression(paste("9",degree)), pos = 1, offset = 0.3)
lines(rbind(c(-11,180), c(-11,171)), lwd = 1)
text(-11, 175, labels = expression(paste("9",degree)), pos = 2, offset = 0.2)
lines(rbind(c(0,0), c(0,180)), lwd = 1) 
text(0, -5, labels = "front", pos = 1, offset = 0)
lines(rbind(c(90,0), c(90,180)), lwd = 1)
text(90, -5, labels = expression(paste("side 90",degree)), pos = 1, offset = 0)
lines(rbind(c(-12,90), c(162,90)), lwd = 1)
text(-17, 90, labels = "equator", pos = 1, offset = 0, srt = 90)

dev.off()



# target neuite vs syn --------------------------------------------------------------------------------------------

# break_x <- seq(340000, 450000, length.out = 23)

nmb <- 16
break_x <- seq(range_x[1] - 2*dL*nmb, range_x[2] + 2*dL*nmb, by = 2*dL*4)

# break_x <- glo_div

tar_syn_xbin <- matrix(ncol = 9, nrow = length(break_x)-1)
for (j in 1:9) {
  tar_syn_xbin[,j] <- hist(conn_LC6_tar[[j]]$x, breaks = break_x, plot = F)$counts
}

# moving 3-avg
tar_syn_xbin3 <- matrix(ncol = 9, nrow = length(break_x)-3)
for (j in 1:nrow(tar_syn_xbin3)) {
  tar_syn_xbin3[j,] <- colMeans(tar_syn_xbin[j:(j+2),]) %>% round(.,1)
}

# tar_syn_xbin <- cbind(hist(conn_LC6_tar[[1]]$x, breaks = break_x, plot = F)$mids, tar_syn_xbin)
xx <- hist(conn_LC6_tar[[1]]$x, breaks = break_x, plot = F)$mids[2:(length(break_x)-2)]
xx <- (xx - range_x[1]) / 2/dL
xx <- xx[3:(length(xx)-2)]
yy <- tar_syn_xbin3
yy <- yy[3:(2+length(xx)), ]
df <- as.data.frame(cbind(xx, yy))
colnames(df) <- c('mid', seq(1,9))
dfm <- melt(df, id = 'mid')

dev.new()
# pdf(file = "#syn per 12um 3bin avg.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(x=mid, y=value) ) + 
  geom_point(size = 3 ) +
  geom_line(aes(colour = factor(variable), group = variable), lwd = 2) +
  # coord_cartesian(ylim = c(-.5, 1)) +
  # theme_minimal_grid() +
  scale_colour_manual(values = pal_target, labels = target_names ) +
  scale_x_continuous(breaks = seq(0,80, by=1), labels = tick_x, expand = c(0,0),limits = c(0,80)) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color='black'),
        axis.ticks.length.x = unit(-.25, "cm"),
        axis.text.x = element_text(margin = margin(t = .3, unit = "cm"))) +
  labs(title = "target syn histo", x='glomerulus axis[um]', y='count') 

dev.off()


# -- node hist
neu_target_rs <- resample(neu_target, stepsize = 400)

tar_node_xbin <- matrix(ncol = 9, nrow = length(break_x)-1)
for (j in 1:9) {
  i1 <- pointsinside(xyzmatrix(neu_target_rs[[j]]$d), surf=glo.msh, rval='distance') > -500
  i2 <- xyzmatrix(neu_target_rs[[j]]$d)[,1] > break_x[1] & xyzmatrix(neu_target_rs[[j]]$d)[,1] < tail(range_x,1)
  xx <- xyzmatrix(neu_target_rs[[j]]$d)[i1 & i2,'X']
  tar_node_xbin[,j] <- hist(xx, breaks = break_x, plot = F)$counts
}


# moving 3-avg
tar_node_xbin3 <- matrix(ncol = 9, nrow = length(break_x)-3)
for (j in 1:nrow(tar_node_xbin3)) {
  tar_node_xbin3[j,] <- colMeans(tar_node_xbin[j:(j+2),]) %>% round(.,1)
}

xx <- hist(conn_LC6_tar[[1]]$x, breaks = break_x, plot = F)$mids[2:(length(break_x)-2)]
xx <- (xx - range_x[1]) / 2/dL
xx <- xx[3:(length(xx)-2)]
yy <- tar_node_xbin3
yy <- yy[3:(2+length(xx)), ] * 0.4 # [um]
df <- as.data.frame(cbind(xx, yy))
colnames(df) <- c('mid', seq(1,9))
dfm <- melt(df, id = 'mid')


dev.new()
# pdf(file = "neurite in um 3bin avg.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(x=mid, y=value) ) + 
  geom_point(size =3 ) +
  geom_line(aes(colour = factor(variable), group = variable), lwd = 2) +
  # coord_cartesian(ylim = c(-.5, 1)) +
  # theme_minimal_grid() +
  scale_colour_manual(values = pal_target, labels = target_names ) +
  scale_x_continuous(breaks = seq(0,80, by=1), labels = tick_x, expand = c(0,0),limits = c(0,80)) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color='black'),
        axis.ticks.length.x = unit(-.25, "cm"),
        axis.text.x = element_text(margin = margin(t = .3, unit = "cm"))) +
  labs(title = "target node histo", x='glomerulus axis[um]', y='cable length [um]') 

dev.off()

# -- ratio
xx <- hist(conn_LC6_tar[[1]]$x, breaks = break_x, plot = F)$mids[2:(length(break_x)-2)]
xx <- (xx - range_x[1]) / 2/dL
xx <- xx[3:(length(xx)-3)]
yy <- tar_syn_xbin3 / tar_node_xbin3
yy <- yy[3:(2+length(xx)), ] * 2.5 #per um
df <- as.data.frame(cbind(xx, yy))
colnames(df) <- c('mid', seq(1,9))
dfm <- melt(df, id = 'mid')


# colnames(df) <- c('mid', seq(1,9))
# df$mid <- hist(xx, breaks = break_x, plot = F)$mids
# dfm <- melt(df, id = 'mid')

dev.new()

# pdf(file = paste("syn density.pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(x=mid, y=value) ) + 
  geom_point(size = 3 ) +
  geom_line(aes(colour = factor(variable), group = variable), lwd = 2) +
  # coord_cartesian(ylim = c(-.5, 1)) +
  # theme_minimal_grid() +
  scale_colour_manual(values = pal_target, labels = target_names ) +
  # scale_x_continuous(breaks = seq(range_x[1],range_x[1]+80*2*dL, by=2*dL), limits = c(340000,450000),labels = tick_x, expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,1, by=0.25), labels = seq(0,1, by=0.25), limits = c(0,1)) +
  scale_x_continuous(breaks = seq(0,80, by=1), labels = tick_x, expand = c(0,0), limits = c(0,80)) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color='black'),
        axis.ticks.length.x = unit(-.25, "cm"),
        axis.text.x = element_text(margin = margin(t = .3, unit = "cm"))) +
  labs(title = "target syn density", x='glomerulus axis[um]', y='syn count per um') 

dev.off()


