# codes used for analyze LC6 and downstream neurons
# !!! open with encoding UTF-8
# All data is saved in the .RData file
# For a paritcular figure, search for, eg. "5A" (main figures) or "S5A" (supplimentary)


# load libraries 
library(natverse)
library(rjson)
library(alphashape3d)
library(tidyverse)
library(RColorBrewer)
library(PNWColors)
library(ggplot2)
library(sf)
library(ggExtra)
library(cowplot)
library(alphahull)
library(reshape2)

setwd("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/LC6_code_R")

# setwd("C:/Users/zhaoa/Dropbox (HHMI)/LC6 downstream paper new/Code")
setwd("./") # set to this code's directory
# source("someFunc.R") # load some useful functions

# clean everythign up.  
rm(list=ls())

#close any open rgl windows
while (rgl.cur() > 0) { rgl.close() }

# set up for 3d plots based on rgl package
rgl::setupKnitr()



# Buchner 71 eye map from Andrew Straw----------------------------------------------------------------------------------------------

Npt <- 699
buchner <- read.csv("data/buchner71_tp.csv", header = FALSE)
buchner <- buchner[1:Npt,]
buchner <- buchner / pi * 180

dev.new()
plot(buchner)

range(buchner[buchner[,2] > 70 & buchner[,2] < 110, 1]) # [-7, 160] as horizontal range
# range(buchner[buchner,2]) 

buchner_phi <- c(-10, 160) # for longitude



# func to generate polygon from set of points ---------------------------------------------------------------------

mkpoly <- function(xy) {
  xy <- as.data.frame(xy)[,1:6]
  xyset <- list() # vertices for a list of polygons
  if (dim(xy)[1] > 2) {
    N <- 1
    xyset[[N]] <- xy[1,]
    xy <- xy[-1, ]
    ii <- c()
    while (dim(xy)[1] >= 1) {
      ii[1] <- match(tail(xyset[[N]], 1)[2], xy[, 1])
      ii[2] <- match(tail(xyset[[N]], 1)[2], xy[, 2])
      if (!is.na(ii[1])) {
        xyset[[N]] <- rbind(xyset[[N]], xy[ii[1], ])
        xy <- xy[-ii[1], ]
      } else if (!is.na(ii[2])){
        xytmp <- xy[ii[2], c(2,1,5,6,3,4)]
        colnames(xytmp) <- colnames(xyset[[N]])
        xyset[[N]] <- rbind(xyset[[N]], xytmp)
        xy <- xy[-ii[2], ]
      } else {
        N <- N + 1
        xyset[[N]] <- xy[1, ]
        xy <- xy[-1, ]
      }
    }
  }
  return(xyset)
}

# load neurons ----------------------------------------------------------------------------------------------------

# LC6 neurons
anno_LC6 <- catmaid_query_by_annotation("^LC6 neuron$")
anno_RHS <- catmaid_query_by_annotation("^LO-R")
anno_rightLC6 <- anno_LC6[anno_LC6[,1] %in% anno_RHS[,1],]
LC6_skid <- anno_rightLC6[,"skid"]
neu <-  read.neurons.catmaid(LC6_skid, .progress='text')
LC6 <- neu

altTract <- c(8,14,19,39)

# target neurons
anno_ipsi <- catmaid_query_by_annotation("^putative ss2036$")
ipsi_skid <- anno_ipsi[,"skid"]
neu_ipsi <-  read.neurons.catmaid(ipsi_skid, .progress='text') 

# anno_biL <- catmaid_query_by_annotation("^LC6 target - bilateral$")
# biL_skid <- anno_biL$skid # LHS x4 , 
# biL_skid <- biL_skid[1:2] #select 2
biL_skid <- c(3149978, 3262999)
neu_biL <- read.neurons.catmaid(biL_skid, .progress='text') 

biR_skid <- c(3154013, 3155729) #RHS x2
neu_biR <- read.neurons.catmaid(biR_skid, .progress='text')

neu_target <- c(neu_biL, neu_biR, neu_ipsi)

# TM5
neu_JSON <- fromJSON(file = "Tm5_LC6 mapping.json")
neu_skid <- c()
for (j in 1:length(neu_JSON)) {
  neu_skid[j] <- neu_JSON[[j]]$skeleton_id
}
TM5 = read.neurons.catmaid(neu_skid, .progress='text')

# load glomerulus volumne mesh 
glo_vol <- catmaid_get_volume("v14.LC6_glomerulus_preSyn_R")

# whole brain mesh
v14 <- catmaid_get_volume(439, rval = 'mesh3d')  

# make a layer using LC6 ----------------------------------------------------------------------------------------

# define a plane separating lobula portion of the LC6
a = -0.6; b = 1; c = 1.3; d = -230000

# get all end points and branch points
for (j in 1:length(neu)) {
  tar <- neu[[j]]
  xyz_ep <-  tar$d[tar$EndPoints, ] %>% xyzmatrix()
  xyz_bp = tar$d[tar$BranchPoints, ] %>% xyzmatrix()
  xyz_LO <- rbind(xyz_ep, xyz_bp) %>% 
    as_tibble() %>% 
    mutate(LO = a*X + b*Y + c*Z + d) %>%
    filter(LO > 0) %>%
    select(X,Y,Z)
  if (j == 1) {
    xyz_node <- xyz_LO
  } else {
    xyz_node <- bind_rows(xyz_node, xyz_LO)
  }
}

# fit curved surface
polyfitorder <- 2 # 2nd order surface
gridMargin <- 20000 #add margin on the edge
xyz_node <- data.matrix(xyz_node)
X <- xyz_node[, 1];  Y <- xyz_node[, 2]; Z <- xyz_node[, 3]
fitlm <- lm(Z ~ poly(X, Y, degree = polyfitorder, raw = TRUE)) #linear model fit

# make an underlying grid to interpret the fit as a point set
dx2 <- 500 
dy2 <- 500
xx <- seq(range(X)[1] - gridMargin, range(X)[2] + gridMargin, by = dx2)
yy <- seq(range(Y)[1] - gridMargin, range(Y)[2] + gridMargin, by = dy2)
xygrid <- expand.grid(xx, yy)
xygrid <- setNames(data.frame(xygrid), c('X', 'Y'));
valfit <- predict(fitlm, xygrid) #generate values from the fit
xyz_lm <- cbind(xygrid, valfit)
dist_min <- apply(xyz_lm, 1, function(pt) {min(rowSums(sweep(xyz_node, 2, pt) ^ 2))}) #calculate min distance
ii <- dist_min > 20e6 # old, for visulization
xyz_layer_25e6 <- xyz_lm[!ii,] # pts
valfit_sup <- valfit
valfit_sup[ii] <- NA
m_surf <-  matrix(valfit_sup, nrow = length(xx), ncol = length(yy)) # surface

ii <- dist_min > 5e6
xyz_layer <- xyz_lm[!ii,] # pts

# coarse version for plotting
dx2c <- 1000
dy2c <- 1000
xxc <- seq(range(X)[1] - gridMargin, range(X)[2] + gridMargin, by = dx2c)
yyc <- seq(range(Y)[1] - gridMargin, range(Y)[2] + gridMargin, by = dy2c)
xygridc <- expand.grid(xxc, yyc)
xygridc <- setNames(data.frame(xygridc), c('X', 'Y'));
# valfitc <- predict(fitlm, xygridc) #generate values from the fit
# xyz_lmc <- cbind(xygridc, valfitc)
dist_min <- apply(xygridc, 1, function(pt) {min(rowSums(sweep(xyz_node[,c(1,2)], 2, pt) ^ 2))}) #calculate min distance
ii <- dist_min > 20e6
xy_layer_coarse <- xygridc[!ii,] # pts

# make alpha mesh
msh.a <- ashape3d(xyz_node, alpha = 20000) # 20000 look ok
msh <- as.mesh3d(msh.a)
neu_lo <- nlapply(neu, subset, function(x) pointsinside(x, msh))

# Figure 5B
nopen3d()
par3d('windowRect' = c(100,100,1700,1700))
plot3d(neu[[10]],  col= "#d7191c", lwd = 5, soma=T, WithNodes = F)
plot3d(neu[[51]],  col= "#2c7bb6", lwd = 5, soma=T, WithNodes = F)
plot3d(neu,  col= 'grey90', soma=T, WithNodes = F)
shade3d(glo_vol, col= "cyan", alpha = 1)
surface3d(xx,yy,m_surf, color = "#fdae61", alpha = 1)
plot3d(TM5[[1]], col = 'gold4', lwd = 5) #TM5
plot3d(TM5[[2]], col = 'brown', lwd = 5)
rgl.viewpoint(userMatrix = rotationMatrix(30/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0))
rgl.pop(type = "light")
rgl.light(theta = 30, phi = 30)

# # save
# rgl.snapshot(filename = "LC6_3d.png",fmt = "png")


# 2d projections --------------------------------------------------------------------------------------------------
# projection dendrites onto the layer, both contour and center-of-mass

# xyz_layer <- xyz_layer
row.names(xyz_layer) <-  seq(1, dim(xyz_layer)[1])

ind_pj <- list() #index of projected grid points
xyz_com <- list() # center-of-mass
xyz_pj_com <- list() # xyz of com projecting on grid
for (j in 1:length(neu)){
  tar <- neu[[j]]
  xyz_ep <-  tar$d[tar$EndPoints, ] %>% 
    xyzmatrix() %>%
    as_tibble() %>% 
    mutate(LO = a*X + b*Y + c*Z + d) %>%
    filter(LO > 0) %>%
    select(X,Y,Z) %>%
    data.matrix()
  xyz_bp <-  tar$d[tar$BranchPoints, ] %>% 
    xyzmatrix() %>%
    as_tibble() %>% 
    mutate(LO = a*X + b*Y + c*Z + d) %>%
    filter(LO > 0) %>%
    select(X,Y,Z) %>%
    data.matrix()
  
  # center-of-mass
  xyz_eb <- rbind(xyz_ep,xyz_bp)
  xyz_com[[j]] <- colSums(xyz_eb)/dim(xyz_eb)[1]
  
  # project dendrite end points and com to the fitted grid by shortest distance
  xyz_dend_pj <- rbind(xyz_com[[j]], xyz_ep) #append com at the beginning
  Nsigma <- 5 #exclude > 5 sigma
  com_dist <- as.matrix(dist(xyz_dend_pj))
  thrhd_dist <- sd(com_dist[,1])*Nsigma
  xyz_dend_pj <- cbind(xyz_dend_pj, com_dist[,1]) %>%
    as_tibble() %>%
    filter(V4 < thrhd_dist)%>%
    select(X,Y,Z) %>%
    data.matrix()
  
  ind_min <- c()
  ind_min <-apply(xyz_dend_pj, 1, function(pt) {which.min(rowSums(sweep(xyz_layer, 2, pt) ^ 2))}) # index in xyz_layer with min distance
  ind_com <- ind_min[1] #index of grid point that's closest to com
  ind_min <- unique(ind_min[-1]) #index of grid points
  xyz_pj_com[[j]] <- xyz_layer[ind_com,]
  ind_pj[[j]] <- row.names(xyz_layer[ind_min,])
}

# projection down onto (x,y) plane
xy_pj <- list()
for (j in 1:length(ind_pj)){
  xy_tmp <- list()
  xy_ashape <- list()
  ii <- sort(as.integer(ind_pj[[j]]))
  xy_tmp <- xyz_layer[ii, 1:2]
  xy_ashape <- ashape(xy_tmp + matrix(runif(dim(xy_tmp)[1]*2, 1e-9, 2e-9), ncol = 2), alpha = 6000)
  xy_edge <- xy_ashape$edges[,1:6]
  xy_pj[[j]] <- list(xy=xy_tmp, ashape=xy_ashape, edge=xy_edge)
}


# ahull
xy_ashape_grid <- ashape(xyz_layer[,1:2] + matrix(runif(dim(xyz_layer)[1]*2, 1e-9, 2e-9), ncol = 2), alpha = 6000)


# # use chull for the grid projection
# xy_grid <- xy_ashape_grid$edges[,c("x1","y1")]
# hpts <- chull(xy_grid)
# hpts <- c(hpts, hpts[1])
# xy_edge_grid <- xy_grid[hpts,] # hull edge points

# ### ### use alpha-hull for the grid projection
xy_grid_ahull <- mkpoly(xy_ashape_grid$edges)[[1]][,3:4]
xy_edge_grid <- xy_grid_ahull # hull edge points

# get edge points for each LC6 projection
xy_poly <- list()
for (j in 1:length(ind_pj)) {
  ls_poly <- mkpoly(xy_pj[[j]]$ashape$edges)
  xy_poly[[j]] <- ls_poly[[1]][,3:4]
}


# use TM5 to determine the center and meridian  -------------------------------------------------------------------

tar <- TM5[[1]]
ng <- as.ngraph(tar)
distal_points <- igraph::graph.dfs(ng, root=2137, unreachable=FALSE, neimode='out')$order
proximal_points <- setdiff(igraph::V(ng), distal_points)
# points3d(xyzmatrix(tar$d[proximal_points,]), col = 'red', size = 5)
TM5_u <- colMeans(xyzmatrix(tar$d[proximal_points,])) #upper pt

tar <- TM5[[2]]
ng <- as.ngraph(tar)
distal_points <- igraph::graph.dfs(ng, root=873, unreachable=FALSE, neimode='out')$order
proximal_points <- setdiff(igraph::V(ng), distal_points)
TM5_c <- colMeans(xyzmatrix(tar$d[proximal_points,])) # central pt

grid_c <- which.min(rowSums((sweep(xyz_layer, 2, c(TM5_c), "-"))^2))
grid_u <- which.min(rowSums((sweep(xyz_layer, 2, c(TM5_u), "-"))^2))


# --- map equator and coord system 
ymax <- max(xy_edge_grid[,2])
ymin <- min(xy_edge_grid[,2])
center_new <- xyz_layer[grid_c,1:2]
x_med_new <- as.numeric(center_new[1])
y_eq_new <- as.numeric(center_new[2])
angR_new <- acos((x_med_new - xyz_layer[grid_u,1])/sqrt(sum((xyz_layer[grid_c,1:2] - xyz_layer[grid_u,1:2])^2)))

# PLOT,PAPER, 2d projection with dendrites
ang_2 <- pi/2 - angR_new
rot_2 <- matrix(c(cos(ang_2), sin(ang_2), -sin(ang_2), cos(ang_2)), ncol = 2)
xy_layer_coarse_rot <- sweep(xy_layer_coarse, MARGIN = 2, STATS = c(x_med_new, y_eq_new))
xy_layer_coarse_rot <- sweep(t(rot_2 %*% t(xy_layer_coarse_rot)), 2, c(x_med_new, y_eq_new), '+')

# boundary
xy_edge_grid_rot <- sweep(xy_edge_grid, MARGIN = 2, STATS = c(x_med_new, y_eq_new))
xy_edge_grid_rot <- sweep(t(rot_2 %*% t(xy_edge_grid_rot)), 2, c(x_med_new, y_eq_new), '+')

# Figure 5C
windows(record = F, width = 8, height = 8)
# pdf('LC6_2d_lo.pdf', width = 8, height = 8, family = "Courier")
plot(xy_layer_coarse_rot, col="#fdae61", cex = 1, pch = ".", ylim = rev(range(xy_layer_coarse_rot[,2])), asp = 1)
for (j in 1:length(neu_lo)) {
  neu_lo_xy <- neu_lo[[j]]$d[,c("X","Y")]
  pp <- sweep(neu_lo_xy, MARGIN = 2, STATS = c(x_med_new, y_eq_new))
  pp <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new), '+')
  # points(neu_lo_xy, col = "grey80", pch = ".", cex = 1)
  points(pp, col = "grey80", pch = ".", cex = 1)
  points(colMeans(pp)[1], colMeans(pp)[2], pch = 20, col = "blue", cex = 1.5)
}
ng=as.ngraph(neu[[10]])
distal_points=igraph::graph.dfs(ng, root=449, unreachable=FALSE, neimode='out')$order
distal_tree=subset(neu[[10]], distal_points)
pp <- sweep(distal_tree$d[,c("X","Y")], MARGIN = 2, STATS = c(x_med_new, y_eq_new))
distal_tree$d[,c("X","Y")] <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new), '+')
plot(distal_tree,  col= "#d7191c", lwd = 2, soma=T, WithNodes = F, add = T)
ng=as.ngraph(neu[[51]])
distal_points=igraph::graph.dfs(ng, root=636, unreachable=FALSE, neimode='out')$order
distal_tree=subset(neu[[51]], distal_points)
pp <- sweep(distal_tree$d[,c("X","Y")], MARGIN = 2, STATS = c(x_med_new, y_eq_new))
distal_tree$d[,c("X","Y")] <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new), '+')
plot(distal_tree,  col= "#2c7bb6", lwd = 2, soma=T, WithNodes = F, add = T)
# lines(rbind(xyz_layer[grid_c,1:2], xyz_layer[grid_c,1:2]+(xyz_layer[grid_u,1:2]-xyz_layer[grid_c,1:2])*2), lwd = 3, col = 'cyan')
pp <- sweep(xyz_layer[grid_u,1:2], MARGIN = 2, STATS = c(x_med_new, y_eq_new))
pp <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new), '+')
lines(rbind(c(x_med_new,y_eq_new), c(x_med_new,y_eq_new)+(pp-c(x_med_new,y_eq_new))*2), lwd = 3, col = 'cyan')
points(matrix(c(x_med_new,y_eq_new), ncol =2), pch = 18, col = 'brown', cex = 2)
points(c(x_med_new,y_eq_new)+(pp-c(x_med_new,y_eq_new)), pch = 18, col = 'gold4', cex = 2)
lines(c(300000,310000), c(300000, 300000), col = 'black', lwd = 3)
text(x = 305000, 295000, labels = "10 µm")
# polygon(xy_edge_grid_rot)


# wrap the lobula layer onto a hemisphere, use line-polygon boundary to calculate radial distance --------------

# get edge points for each LC6 projection
xy_poly <- list()
for (j in 1:length(ind_pj)) {
  ls_poly <- mkpoly(xy_pj[[j]]$ashape$edges)
  xy_poly[[j]] <- ls_poly[[1]][,3:4]
}

grid_bdpt <- xy_edge_grid
grid_bdpt <- rbind(grid_bdpt, grid_bdpt[1,])

library(sf)

poly_st <- st_polygon(list(data.matrix(grid_bdpt)))
xy_ori <- c(x_med_new, y_eq_new)
R <- (ymax-ymin)*2
line_st <- st_linestring(t(matrix(c(x_med_new, y_eq_new, x_med_new+100000, y_eq_new+100000), nrow = 2)))
int_st = st_intersection(line_st, poly_st) # the 2nd point is the intersection
xy_com <- list()
for (j in 1:length(xy_poly)) {
  xy_pj_com <- xyz_pj_com[[j]][c('X','Y')]
  colnames(xy_pj_com) <- c("x1", "y1")
  xy_poly[[j]] <- rbind(xy_pj_com, xy_poly[[j]])
  xy_poly[[j]] %<>% 
    as_tibble() %>%
    mutate(phiC = acos((x1-x_med_new)/sqrt((x1-x_med_new)^2+(y1-y_eq_new)^2))*(-1)^(y1<y_eq_new) + 2*pi*(y1<y_eq_new)) %>%  #angle
    mutate(thetaC = NA) %>% #radius
    transmute(x1, y1, phiC, thetaC) %>%
    data.matrix()
  for (k in 1:dim(xy_poly[[j]])[1]) {
    alpha <- xy_poly[[j]][k,'phiC']
    line_st <- st_linestring(rbind(xy_ori, c(xy_ori[1] + R*cos(alpha), xy_ori[2] + R*sin(alpha))))
    int <- data.matrix(st_intersection(line_st, poly_st))[2,]
    xy_poly[[j]][k,'thetaC'] <- pi/2 * dist(rbind(xy_poly[[j]][k,1:2], xy_ori)) / dist(rbind(int,xy_ori))
  }
  # now turn inside out and a 90-rotation about x-axis to have the front edge on the x-z plane, 
  # angle wrt looking from behiind the eye
  xy_poly[[j]] %<>%
    as_tibble() %>%
    mutate(thetaC = pi - thetaC) %>% # turn inside out
    mutate(phiC = phiC - pi/2 - angR_new) %>% # align front edge to x-z plane
    mutate(x = sin(thetaC)*cos(phiC), y = sin(thetaC)*sin(phiC), z = cos(thetaC)) %>% #(x,y,z) coord
    mutate(xr = x, yr = -z, zr = y) %>% # +90 rotation about x-axis
    mutate(xrr = xr, yrr = yr, zrr = zr) %>% # do nothing
    mutate(theta = acos(zrr/sqrt(xrr^2+yrr^2+zrr^2)), phi = acos(xrr/sqrt(xrr^2+yrr^2))) %>%
    mutate(theta_deg = theta/pi*180, phi_deg = phi/pi*180/180*diff(buchner_phi)+buchner_phi[1]) %>% # longitude use buchner_phi
    select(x1, y1, theta_deg, phi_deg) %>%
    data.matrix()
  xy_com[[j]] <- xy_poly[[j]][1,]
  xy_poly[[j]] <- xy_poly[[j]][-1,]
}


# Figure 5D
windows(record = F, width = 8, height = 8)
bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
bd_theta <- seq(1, 180, by = 1)
xy_bd <- matrix(ncol = 2)
bd_grid <- expand.grid(bd_phi, bd_theta)
plot(bd_grid, ylim = rev(range(bd_grid$Var2)), type = "n", axes = FALSE, ann = F)
for (j in 1:length(xy_poly)) {
  xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
  if (j == 10) {
    polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "#d7191c", density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 3, pch = 3)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
  }
  else if (j == 51) {
    polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "#2c7bb6", density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 3, pch = 3)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
  }
  else {
    # polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "cyan", border = 'black', density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
    text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = j, pos = 1, offset = 0.3)
  }
}
xy_bd <- xy_bd[-1,]
hpts_2 <- chull(xy_bd)
hpts_2 <- c(hpts_2, hpts_2[1])
xy_bd_chull <- xy_bd[hpts_2,] # hull edge points
polygon(xy_bd_chull)
lines(rbind(c(-11,180), c(-2,180)), lwd = 3)
text(-5, 180, labels = expression(paste("9",degree)), pos = 1, offset = 0.3)
lines(rbind(c(-11,180), c(-11,171)), lwd = 3)
text(-11, 175, labels = expression(paste("9",degree)), pos = 2, offset = 0.2)
lines(rbind(c(0,0), c(0,180)), lwd = 3) 
text(0, -5, labels = "front", pos = 1, offset = 0)
lines(rbind(c(90,0), c(90,180)), lwd = 3)
text(90, -5, labels = expression(paste("side 90",degree)), pos = 1, offset = 0)
lines(rbind(c(-12,90), c(162,90)), lwd = 3)
text(-17, 90, labels = "equator", pos = 1, offset = 0, srt = 90)


# Figure S6A
bd_phi2 <- seq(buchner_phi[1], buchner_phi[2], by = 2)
bd_theta2 <- seq(1, 180, by = 2)
pt_grid <- expand.grid(bd_phi2, bd_theta2) 
pt_grid <- cbind(pt_grid, rep(0, dim(pt_grid)[1]))
colnames(pt_grid) <- c('x','y','paint')
for (j in 1:length(xy_poly)) {
  ii <- sp::point.in.polygon(pt_grid[,"x"], pt_grid[,"y"], xy_poly[[j]][,"phi_deg"], xy_poly[[j]][,"theta_deg"])
  pt_grid[,"paint"] <- pt_grid[,"paint"] + ii
}
pt_grid <- pt_grid[pt_grid[,"paint"] != 0, ] 
pt_grid[,"paint"] <- factor(pt_grid[,"paint"], labels = seq(1,5, by = 1))

# PLOT
windows(record = F, width = 8, height = 8)
xy_bd_chull_df <- as.data.frame(xy_bd_chull)
x_tick <- c(0, 45, 90, 135, 180) + buchner_phi[1]
y_tick <- c(0, 45, 90, 135, 180)
dev.new()
ggplot() + 
  geom_point(data = pt_grid, aes(x = x, y = y, colour = paint)) +
  scale_color_brewer(palette =  "Blues") +
  scale_y_reverse() +
  guides(col = guide_legend(title = "groupings")) +
  geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg)) +
  geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
  annotate("text", x = -5, y = 185, label = "9°") +
  geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
  annotate("text", x = -15, y = 175.5, label = "9°") +
  geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
  geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
  geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
  scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
  scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
  theme(axis.title = element_blank(), axis.text = element_blank()) +
  theme_void() +
  labs(title = "2d projection on LO")




# LC6 glomerulus ------------------------------------------------------------------------------------------------

# all the connectors on LC6
LC6_pre <- data.frame()
LC6_post <- data.frame()
conn_glu <- list()
conn_pre <- list()
conn_post <- list()
for (j in 1:length(neu)) {
  tar <- neu[[j]]
  conn_glu[[j]] <- tar$connectors %>% 
    as_tibble() %>%
    mutate(glu = a*x + b*y + c*z + d) %>%
    filter(glu < 0)
  conn_pre[[j]] <- conn_glu[[j]] %>% 
    filter(prepost == 0) %>%
    select(treenode_id, connector_id, x, y, z) %>%
    as.data.frame()
  conn_post[[j]] <- conn_glu[[j]] %>% 
    filter(prepost == 1) %>%
    select(treenode_id, connector_id, x, y, z) %>%
    as.data.frame()
  LC6_pre <- rbind(LC6_pre, cbind(rep(j, dim(conn_pre[[j]])[1]),conn_pre[[j]]))
  LC6_post <- rbind(LC6_post, cbind(rep(j, dim(conn_post[[j]])[1]),conn_post[[j]]))
}
colnames(LC6_pre)[1] <- "ID"
colnames(LC6_post)[1] <- "ID"

if (exists("LC6LC6")) {
  rm(LC6LC6)
}
for (j in 1:length(neu)) {
  # tar_pre <- neu[[j]]
  for (k in 1:length(neu)) {
    # tar_post <- neu[[k]]
    cft <- catmaid_get_connectors_between(pre_skids = neu[[j]]$skid, post_skids = neu[[k]]$skid)
    if (!is.null(cft)) {
      cft_nrow <- dim(cft)[1]
      cft_comb <- cbind(rep(j,cft_nrow), rep(k,cft_nrow), cft[,1:8])
      if (exists("LC6LC6")) {
        LC6LC6 <- rbind(LC6LC6,cft_comb)
      } else {
        LC6LC6 <- cft_comb
      }
    }
  }
}
LC6LC6 <- as.data.frame(LC6LC6)
colnames(LC6LC6) <- c('pre_ind','post_ind','pre_skid','post_skid','conn_id','pre_node_id','post_node_id','x','y','z')
rownames(LC6LC6) <- seq(1,dim(LC6LC6)[1])


# -- look at pairwise distance vs pairwise connection
a1 = -0.6; b1 = 1; c1 = 1.3; d1 = -200000 #separate LO and glomerulus
dist_Nconn <- matrix(ncol = 7) #pairwise dist and connections
# colnames(dist_Nconn) <- c("from","to", "dist", "fromto_LO", "fromto_glu", "tofrom_LO", "tofrom_glu")
for (j in 1:(length(conn_pre)-1)) {
  for (k in (j+1):length(conn_pre)) {
    dist_d <- dist(rbind(xyz_com[[j]], xyz_com[[k]]))
    conn_fromto <- catmaid_get_connectors_between(pre_skids = neu[[j]]$skid, post_skids = neu[[k]]$skid)
    if (is.null(conn_fromto)) {
      fromto_LO <- 0
      fromto_glu <- 0
    } else {
      mat_conn <- data.matrix(conn_fromto[,c("connector_x", "connector_y", "connector_z")])
      fromto_LO <- sum(mat_conn %*% c(a1,b1,c1) + d1 > 0)
      fromto_glu <- dim(conn_fromto)[1] - fromto_LO
    }
    conn_tofrom <- catmaid_get_connectors_between(pre_skids = neu[[k]]$skid, post_skids = neu[[j]]$skid)
    if (is.null(conn_tofrom)) {
      tofrom_LO <- 0
      tofrom_glu <- 0
    } else {
      mat_conn <- data.matrix(conn_tofrom[,c("connector_x", "connector_y", "connector_z")])
      tofrom_LO <- sum(mat_conn %*% c(a1,b1,c1) + d1 > 0)
      tofrom_glu <- dim(conn_tofrom)[1] - tofrom_LO
    }
    dist_Nconn_tmp <- matrix(c(j,k,dist_d, fromto_LO, fromto_glu, tofrom_LO, tofrom_glu), nrow = 1)
    dist_Nconn <- rbind(dist_Nconn, dist_Nconn_tmp)
  }
}
dist_Nconn <- dist_Nconn[-1,]
colnames(dist_Nconn) <- c("from","to","dist_com","fromto_LO","fromto_glu","tofrom_LO","tofrom_glu")
dist_Nconn %<>% 
  as_tibble() %>%
  mutate(Nconn_glu = tofrom_glu + fromto_glu) %>%
  as.data.frame()
colSums(dist_Nconn)
dist_com_mean <- mean(dist_Nconn$dist_com)
Nconn_glu_tot <- sum(dist_Nconn$Nconn_glu)


# dist to a neuron avg over all partners weighted by num of conn
dist_Nconn %<>% as_tibble() %>%
  mutate(distxN = dist_com*Nconn_glu) %>%
  as.data.frame()
LC6_avgdist_wt <- matrix(ncol = 2, nrow = length(conn_pre))
LC6_avgdist <- c()
LC6_nnbdist <- c()
# N_nb <- c()
for (j in 1:length(conn_pre)) {
  tib_tmp <- dist_Nconn %>%
    as_tibble() %>%
    mutate(bool = from == j | to == j) %>%
    filter(bool == TRUE) %>%
    # transmute(Nconn_glu, distxN) %>%
    data.matrix()
  # N_nb <- c(N_nb, sum(tib_tmp[,"distxN"] != 0))
  LC6_avgdist[j] <- sum(tib_tmp[,"dist_com"])/dim(tib_tmp)[1]
  LC6_nnbdist[j] <- min(tib_tmp[,"dist_com"])
  LC6_avgdist_wt[j,] <- c(j, sum(tib_tmp[,"distxN"])/sum(tib_tmp[,"Nconn_glu"]))
}
# avg_nb <- mean(N_nb)
# randomize connection partners
dist_nonzero <- dist_Nconn %>%
  filter(distxN != 0) %>%
  transmute(Nconn_glu, distxN) %>%
  data.matrix()
mat_avgdist_rand <- matrix(ncol = 2, nrow = length(conn_pre))
for (j in 1:length(conn_pre)) {
  # -- maintain connection num for given LC6
  ii_rand <- sample.int(length(conn_pre)-1)
  tib_tmp <- dist_Nconn %>%
    as_tibble() %>%
    mutate(bool = from == j | to == j) %>%
    filter(bool == TRUE) %>%
    select(dist_com, Nconn_glu) %>%
    mutate(dist_com = dist_com[ii_rand]) %>%
    mutate(distxN = dist_com*Nconn_glu) %>%
    as.data.frame()
    mat_avgdist_rand[j,] <- c(j, sum(tib_tmp[,"distxN"])/sum(tib_tmp[,"Nconn_glu"]))
}


# Figure 5E
dat_ggplot <- data.frame(rbind(LC6_avgdist_wt, mat_avgdist_rand))
gp  <- factor(c(rep(1,length(conn_pre)),rep(2,length(conn_pre))), labels = c("EM data","Randomized"))
dat_ggplot <- cbind(dat_ggplot, gp)
colnames(dat_ggplot) <- c("neu", "dist", "group")
dat_ggplot$dist <- dat_ggplot$dist/1000

windows(record = F, width = 8, height = 8)
p <- ggplot(dat_ggplot, aes(x = neu, y = dist, colour = gp)) + 
  geom_point(shape = 16, size = 3) +
  theme_void() +
  theme(legend.position = "top", plot.margin = unit(c(3,-5.5,4,3), "mm")) + 
  labs(title = "Avg distance over connection num") + 
  xlim(0, 65) + 
  ylim(30, 80) +
  xlab("neurons index") + 
  ylab("Avg dist [um]") +
  coord_fixed(ratio = 4) +
  geom_hline(yintercept = mean(LC6_nnbdist)/1000, linetype = 2) +
  geom_hline(yintercept = mean(LC6_avgdist)/1000, linetype = 2) +
  scale_colour_discrete(name="Groups")
ggMarginal(p, margins = "y", size = 2, type = "boxplot", outlier.size =3, groupColour = TRUE, groupFill = TRUE)

# KS test
ks.test(LC6_avgdist_wt, mat_avgdist_rand)

# Mann-Whitney
wilcox.test(LC6_avgdist_wt, mat_avgdist_rand, alternative = 'less')


# --- glomerulus dolphin compartment

N_gp <- 11

# for making ahull
conn_LC6 <- LC6_pre # only pre-synapses to divide the dolphin
glu_div <- quantile(conn_LC6$x, probs = seq(0,1,length.out = N_gp)) #separate into divisions of equal pts based on x
# glu_div <- quantile(LC6_pre_glo_all[,1], probs = seq(0,1,length.out = 11)) 
glu_div[1] <- glu_div[1]-1
conn_LC6 <- conn_LC6 %>%
  as_tibble() %>%
  mutate(gp_x = 0)
for (j in 1:length(glu_div)) {
  conn_LC6 %<>% as_tibble() %>%
    mutate(gp_x = gp_x + (x>glu_div[j]))
}

# assign group
conn_LC6LC6 <- LC6LC6
conn_LC6LC6 <- conn_LC6LC6 %>%
  as_tibble() %>%
  mutate(gp_x = 0)
for (j in 1:length(glu_div)) {
  conn_LC6LC6 %<>% as_tibble() %>%
    mutate(gp_x = gp_x + (x>glu_div[j]))
}

ID_wt <- list() #neuron skid with synapse weight in each division
for (j in 1:length(glu_div)) {
  df_tmp <- conn_LC6LC6 %>%
    filter(gp_x == j) %>%
    # transmute(ID) %>%
    data.frame()
  mat_tmp <- matrix(ncol = 2, nrow = length(neu))
  for (k in 1:length(neu)) {
    mat_tmp[k,] <- c(k, sum(df_tmp$pre_ind == k | df_tmp$post_ind == k))  
  }
  ID_wt[[j]] <- mat_tmp
}

# avg radius of projection
avg_size_xy <- c()
for (j in 1:length(neu)) {
  tmp <- as.data.frame(xy_poly[[j]])
  avg_size_xy <- c(avg_size_xy, (max(tmp$theta_deg)-min(tmp$theta_deg)+max(tmp$phi_deg)-min(tmp$phi_deg))/4)
}
r_xy <- mean(avg_size_xy)
# half-width sqrt(-log(0.5))*r_xy*2
# plot Gaussian around each com with synap# as height
ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])

n_lvl <- 11 # 10 compartments
plvl <- list() # ggplot
grid_Gaussian_ls <- list()
grid_Gaussian_sum <- matrix(ncol = 4)
colnames(grid_Gaussian_sum) <- c("X","Y","Z","gp")
plvl_sum <- ggplot()
getPalette <- colorRampPalette(brewer.pal(9, "RdYlBu"))
dolphin_col <- getPalette(n_lvl - 1)
pal_dolphin <- c()
gprange <- list()
for (j in 1:(length(glu_div)-1)) {
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(conn_pre)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- ID_wt[[j]][k,2]
    # grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
    grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
  }
  
  breaks <- seq(0,180,length.out = 9)
  grid_Gaussian$equalSpace <- cut(grid_Gaussian$Z, breaks) #cut into 18 levels
  grid_Gaussian <- cbind(grid_Gaussian, rep(j,dim(grid_Gaussian)[1]))
  colnames(grid_Gaussian) <- c("X","Y","Z","equalSpace","gp")

  plvl[[j]] <- ggplot(grid_Gaussian, aes(X, Y, z= Z)) + 
    geom_raster(aes(fill = equalSpace), interpolate = T) +
    geom_contour(color = "black", alpha = 0.5) +
    theme_void() +
    theme(axis.title = element_blank()) +
    scale_fill_manual(values = (brewer.pal(8, 'Blues')), breaks = seq(0, 190, length.out = 8)) +
    scale_y_reverse() +
    xlim(0, 180) +
    ylim(0, 180) +
    labs(title = paste("Group", j))
  
  cutoffZ <- max(grid_Gaussian$Z) 
  cutoff_per <- 0.5
  grid_Gaussian <- grid_Gaussian[grid_Gaussian$Z > cutoffZ*cutoff_per, ]
  grid_Gaussian_sum <- rbind(grid_Gaussian_sum, grid_Gaussian[,c("X","Y","Z","gp")])
  grid_Gaussian_ls[[j]] <- grid_Gaussian[,c("X","Y","Z","gp")]
  gprange[[j]] <- range(grid_Gaussian_ls[[j]]$Z)
}
grid_Gaussian_sum <- grid_Gaussian_sum[-1,]
grid_Gaussian_sum$gp <- factor(grid_Gaussian_sum$gp)


# plot, all panels, 
dev.new()
plot_grid(plvl[[1]],
          plvl[[2]],
          plvl[[3]],
          plvl[[4]],
          plvl[[5]],
          plvl[[6]],
          plvl[[7]],
          plvl[[8]],
          plvl[[9]],
          plvl[[10]])


# Figure 5F
# color contours
xy_bd_chull_df <- as.data.frame(xy_bd_chull)
gg_cont <- ggplot(grid_Gaussian_sum, aes(X, Y, z = Z)) +
  coord_fixed(ratio = 1) +
  scale_y_reverse()
for (j in 1:(length(glu_div)-1)) {
  gg_cont <- gg_cont + 
    geom_contour(data = grid_Gaussian_ls[[j]], aes(X,Y,z=Z), breaks = seq(gprange[[j]][1], gprange[[j]][2], length.out = 3), color = dolphin_col[j], alpha = 0.9, lwd = 1)
}
gg_cont <- gg_cont + 
  geom_polygon(data = xy_bd_chull_df, aes(x=phi_deg, y=theta_deg, z=0),colour = "black", alpha = 0) +
  theme_void() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
  geom_segment(aes(x = -8, y = 90, xend = 160, yend = 90), size = 2) +
  geom_segment(aes(x = 10, y = 180, xend = 19, yend = 180), size=2, lineend = "round") +
  annotate("text", x = 14.5, y = 185, label = "9°") +
  geom_segment(aes(x = 10, y = 180, xend = 10, yend = 171), size=2, lineend = "round") +
  annotate("text", x = 5, y = 175.5, label = "9°") +
  annotate("text", x = 50, y = 0, label = paste(sprintf("%.1f %%", 100*cutoff_per)))

# make alpha hull
conn_LC6_hull <- as.data.frame(conn_LC6)[c(-5646,-5679,-5640,-5631,-5654,-5652),]
gp_x <- factor(conn_LC6_hull[,"gp_x"], labels = seq(1,N_gp-1, by = 1))
conn_LC6_hull[,"gp_x"] <- gp_x
ii_ahull <- ashape(unique(conn_LC6_hull[, c("x","y")]), alpha = 7000)$edges[, 3:6]
ii_ahull <- as.data.frame(ii_ahull)
# make separate alpha hulls
ii_ahull_ls <- list()
for (j in 1:(length(glu_div)-1)) {
  ii_ahull_ls[[j]] <- ashape(unique(conn_LC6_hull[conn_LC6_hull$gp_x == j, c("x","y")]), alpha = 7000)$edges[, 3:6]
  ii_ahull_ls[[j]] <- as.data.frame(ii_ahull_ls[[j]])
}

# com of each division
range_syn <- matrix(nrow = length(glu_div)-1, ncol = 4)
for (j in 1:(length(glu_div)-1)) {
  tmp <- conn_LC6 %>%
    as_tibble() %>%
    filter(gp_x == j) %>%
    transmute(x, y) %>%
    as.data.frame() 
  range_syn[j,1:2] <- range(tmp$x)
  range_syn[j,3:4] <- range(tmp$y)
}

# num of syn in each div
N_syn_LC6 <- c()
for (j in 1:(length(glu_div)-1)) {
  N_syn_LC6[j] <- sum(conn_LC6LC6$gp_x == j)
}

# dolphin in colors
conn_LC6LC6 <- as.data.frame(conn_LC6LC6)
gp_x <- factor(conn_LC6LC6[,"gp_x"], labels = seq(1,N_gp-1, by = 1))
conn_LC6LC6[,"gp_x"] <- gp_x

pglu <- ggplot(conn_LC6LC6) + 
  geom_point(aes(x = x, y = y, colour = gp_x), shape = 16)
for (j in 1:(length(glu_div)-1)) { #add ahulls
  pglu <- pglu + 
    geom_segment(data = ii_ahull_ls[[j]], aes(x = x1, y = y1, xend = x2, yend = y2)) +
    annotate("text", x = mean(range_syn[j,c(1,2)]), y = mean(range_syn[j,c(3,4)]), label = paste(sprintf("%.1f %%", 100*N_syn_LC6[j]/sum(N_syn_LC6))))
}
pglu <- pglu +
  scale_colour_manual(values = dolphin_col) +
  guides(col = guide_legend(title = "groupings")) +
  theme(legend.position = "bottom") +
  theme_bw() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  scale_y_reverse() +
  ylim(205000, 180000)+
  coord_fixed(ratio = 1) +
  geom_segment(aes(x = 380000, y = 205000, xend = 390000, yend = 205000), size=2, lineend = "round") +
  annotate("text", x = 385000, y = 204000, label = "10 µm") + 
  labs(title = "LC6 Synapses distribution in glomerulus")

windows(record = F, width = 16, height = 10)
ggdraw() +
  draw_plot(pglu + theme(legend.justification = "bottom"), 0, 0.51, 1, 0.5) +
  draw_plot(gg_cont + theme(axis.title = element_text()), 0, 0, 0.5, 0.5)



# Figure S6B, 10 individual compartments
ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])
plvl <- list()
simdata <- list()
simdata_df <- list()
mat_names <- paste("compartment", seq(1,10),sep = " ")

for (j in 1:(n_lvl-1)) {
  getPalette <- colorRampPalette(c("white", dolphin_col[j]))
  pal_tar_red <- getPalette(n_lvl - 1)
  
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(conn_pre)) {
    x0 <- (xy_com[[k]]["phi_deg"])
    y0 <- (xy_com[[k]]["theta_deg"])
    A <- ID_wt[[j]][k,2]
    grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
  }
  simdata[[j]] <- grid_Gaussian
  simdata_df[[j]] <- simdata[[j]]
  colnames(simdata_df[[j]]) <- c("x","y","z")
  simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
  simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  plvl[[j]] <- ggplot(simdata_df[[j]], aes(x, y, z = z)) +
    geom_raster(aes(fill = equalSpace), interpolate = F) +
    scale_fill_manual(values = pal_tar_red, guide=FALSE) +
    theme_void() +
    geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE, lwd = 2) +
    # geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
    # annotate("text", x = -5, y = 185, label = "9°") +
    # geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
    # annotate("text", x = -15, y = 175.5, label = "9°") +
    geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
    geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
    geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
    scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
    scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
    # labs(title = mat_names[j]) +
    geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = seq(0,1,by = 0.1), color = dolphin_col[j], alpha = 1, lwd = 1)+
    coord_fixed(ratio = 1)
}

dev.new()
plvl[[1]]
dev.new()
plvl[[2]]
dev.new()
plvl[[3]]
dev.new()
plvl[[4]]
dev.new()
plvl[[5]]
dev.new()
plvl[[6]]
dev.new()
plvl[[7]]
dev.new()
plvl[[8]]
dev.new()
plvl[[9]]
dev.new()
plvl[[10]]



# target neuron RF wt by synapses ---------------------------------------------------------------------------------

a1 = -0.6; b1 = 1; c1 = 1.3; d1 = -200000 #separate LO and glomerulus
conn_target <- list()
for (j in 1:length(neu_target)) {
  tb_conn <- matrix(ncol = 7)
  for (k in 1:length(neu)) {
    conn_fromto <- catmaid_get_connectors_between(pre_skids = neu_target[[j]]$skid, post_skids = neu[[k]]$skid)
    if (is.null(conn_fromto)) {
      fromto_LO <- 0
      fromto_glu <- 0
    } else {
      mat_conn <- data.matrix(conn_fromto[, c("connector_x", "connector_y", "connector_z")])
      fromto_LO <- sum(mat_conn %*% c(a1, b1, c1) + d1 > 0)
      fromto_glu <- dim(conn_fromto)[1] - fromto_LO
    }
    conn_tofrom <- catmaid_get_connectors_between(pre_skids = neu[[k]]$skid, post_skids = neu_target[[j]]$skid)
    if (is.null(conn_tofrom)) {
      tofrom_LO <- 0
      tofrom_glu <- 0
    } else {
      mat_conn <- data.matrix(conn_tofrom[, c("connector_x", "connector_y", "connector_z")])
      tofrom_LO <- sum(mat_conn %*% c(a1, b1, c1) + d1 > 0)
      tofrom_glu <- dim(conn_tofrom)[1] - tofrom_LO
    }
    Nconn_glu <- fromto_glu + tofrom_glu
    tb_conn_tmp <- matrix(c(j, k, fromto_LO, fromto_glu, tofrom_LO, tofrom_glu, Nconn_glu), nrow = 1)
    tb_conn <- rbind(tb_conn, tb_conn_tmp)
  }
  conn_target[[j]] <- as.data.frame(tb_conn[-1,])
  colnames(conn_target[[j]]) <- c("target","LC6","fromto_LO","fromto_glu","tofrom_LO","tofrom_glu", "Nconn_glu")
}

conn_tt <- matrix(ncol = length(neu_target), nrow = length(neu_target)) #target to target
for (j in 1:length(neu_target)) {
  for (k in 1:length(neu_target)) {
    conn_fromto <- catmaid_get_connectors_between(pre_skids = neu_target[[j]]$skid, post_skids = neu_target[[k]]$skid)
    if (!is.null(conn_fromto)) {
      mat_conn <- data.matrix(conn_fromto[, c("connector_x", "connector_y", "connector_z")])
      fromto <- dim(conn_fromto)[1]
      conn_tt[j,k] <- fromto
    }
  }
}


# Figure S7
# PLOT,  Gaussian around each com with synap# as height,  cp data via binning
ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])
plvl <- list()
simdata <- list()
simdata_df <- list()
LC6_tar_median <- list()
for (j in 1:length(conn_target)) {
  LC6_tar_median[[j]] <- (quantile(conn_target[[j]]$tofrom_glu, c(0.0)))
}
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
    A <- conn_target[[j]][k,"tofrom_glu"]
    if (A >= LC6_tar_median[[j]]) { # selected neuron
      grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
    }
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


# Figure 6C
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
      A <- conn_target_agglo_sum[[j]][k,"tofrom_glu"]
      if (A > 0) { # selected neuron
        grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
      }
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
      geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
      geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
      geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
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




# # physiology data --------------------------------------------------------------------------------------------
# 
# library(R.matlab)
# expBi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss825_mean_RFmap_n=4.mat") #bi
# expBi <- as.matrix(expBi[[1]])
# expIpsi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss2036_mean_RFmap_n=7.mat") #ipsi
# expIpsi <- as.matrix(expIpsi[[1]])
# 
# slope_ex <- (90-83.28)/(117-16)
# (90-83.28)/(117-16)*(9+16)
# 
# n_lvl <- 11
# breaks_3 <- seq(0,1,length.out = n_lvl)
# getPalette <- colorRampPalette(brewer.pal(9, "Oranges"), space = "Lab", bias = 1.8)
# pal_tar_ex <- getPalette(n_lvl - 1)
# 
# x2 <- seq(0-4.25,117-4.25,by = 9) # -18 to 99, with sim data
# y2 <- seq(54, 108,by = 9)
# xygrid2 <- expand.grid(x2, y2)
# 
# expBi_df <- data.frame(xygrid2, as.vector(t(expBi)))
# colnames(expBi_df) <- c("x","y","z")
# expBi_df$equalSpace <- cut(expBi_df$z, breaks_3)
# expBi_df$z <- expBi_df$z*1
# dev.new()
# ggplot(expBi_df, aes(x, y, z = z)) +
#   geom_raster(data = expBi_df, aes(x,y,fill = equalSpace), interpolate = F) +
#   scale_fill_manual(values = pal_tar_ex, guide_legend("dF/F"), labels = paste(seq(0.1,1,length.out = 10))) +
#   geom_segment(aes(x = -9, y = 112, xend = 0, yend = 112), size=2, lineend = "round") +
#   annotate("text", x = 0, y = 116, label = "9°") +
#   geom_segment(aes(x = -9, y = 112, xend = -9, yend = 103), size=2, lineend = "round") +
#   annotate("text", x = -14.5, y = 102, label = "9°") +
#   # geom_segment(aes(x = 16, y = 49.5, xend = 16, yend = 112.5), size=2, lineend = "round") +
#   geom_segment(aes(x = -9, y = 90-lefty, xend = 16, yend = 90), size=2, lineend = "round") +
#   geom_segment(aes(x = 16, y = 90, xend = 117, yend = 90-righty), size=2, lineend = "round") +
#   # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
#   scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
#   scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
#   geom_segment(aes(x = 16, y = 49.5, xend = 16, yend = 112.5), size = 2, lineend = "round") +
#   coord_fixed(ratio = 1)
# 
# expIpsi_df <- data.frame(xygrid2, as.vector(t(expIpsi)))
# colnames(expIpsi_df) <- c("x","y","z")
# expIpsi_df[14, "z"] <- expIpsi_df[14, "z"] - 0.01
# expIpsi_df$equalSpace <- cut(expIpsi_df$z, breaks_3)
# expIpsi_df$z <- expIpsi_df$z*1
# dev.new()
# ggplot(expIpsi_df, aes(x, y, z = z)) +
#   # geom_raster(aes(fill = equalSpace), interpolate = F) +
#   # scale_fill_manual(values = pal_tar_ex, guide_legend("synp den"), labels = paste(seq(0.1,1,length.out = 10))) +
#   geom_raster(data = expIpsi_df, aes(x,y,fill = equalSpace), interpolate = F) +
#   scale_fill_manual(values = pal_tar_ex, guide_legend("dF/F"), labels = paste(seq(0.1,1,length.out = 10))) +
#   geom_segment(aes(x = -9, y = 112, xend = 0, yend = 112), size=2, lineend = "round") +
#   annotate("text", x = 0, y = 116, label = "9°") +
#   geom_segment(aes(x = -9, y = 112, xend = -9, yend = 103), size=2, lineend = "round") +
#   annotate("text", x = -14.5, y = 102, label = "9°") +
#   # geom_segment(aes(x = 16, y = 49.5, xend = 16, yend = 112.5), size=2, lineend = "round") +
#   geom_segment(aes(x = -9, y = 90-lefty, xend = 16, yend = 90), size=2, lineend = "round") +
#   geom_segment(aes(x = 16, y = 90, xend = 117, yend = 90-righty), size=2, lineend = "round") +
#   # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
#   scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
#   scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
#   geom_segment(aes(x = 16, y = 49.5, xend = 16, yend = 112.5), size = 2, lineend = "round") +
#   labs(title = paste("Target ipsi exp")) +
#   coord_fixed(ratio = 1)
# 
# 
# 
# 
# # Figure 6D
# # -- add up the same type
# mat_names <- c("bi_all", "ipsi_all")
# conn_target_agglo_sum <- list()
# conn_target_agglo_sum[[1]] <- conn_target[[1]] + conn_target[[2]] + conn_target[[3]] + conn_target[[4]]
# conn_target_agglo_sum[[2]] <- conn_target[[5]]+conn_target[[6]]+conn_target[[7]]+conn_target[[8]]+conn_target[[9]]
# 
# 
# ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,1], xy_bd_chull[,2])
# plvl <- list()
# simdata <- list()
# simdata_df <- list()
# LC6_tar_median <- list()
# 
# for (j in 1:length(conn_target_agglo_sum)) {
#   grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
#   ii <- grid_Gaussian[,1] > 16
#   grid_Gaussian[ii,2] <- slope_ex*(grid_Gaussian[ii,1] - 16) + grid_Gaussian[ii,2]
#   grid_Gaussian$Z = 0
#   colnames(grid_Gaussian) <- c("X","Y","Z")
#   for (k in 1:length(neu)) {
#     x0 <- (xy_com[[k]]["phi_deg"])
#     y0 <- (xy_com[[k]]["theta_deg"])
#     A <- conn_target_agglo_sum[[j]][k,"tofrom_glu"]
#     # if (A >= LC6_tar_median[[j]][j2]) { # selected neuron
#     grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
#     # }
#   }
#   grid_Gaussian[ii,2] <- round(- slope_ex*(grid_Gaussian[ii,1] - 16) + grid_Gaussian[ii,2])
#   simdata[[j]] <- grid_Gaussian
#   simdata_df[[j]] <- simdata[[j]]
#   colnames(simdata_df[[j]]) <- c("x","y","z")
#   simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
#   simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
#   simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
#   simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
# 
#   if (j <= 1) {
#     plvl[[j]] <- ggplot(expBi_df, aes(x, y, z = z))+
#       geom_raster(data = expBi_df, aes(x,y,fill = equalSpace), interpolate = F) +
#       scale_fill_manual(values = pal_tar_ex, guide_legend("dF/F"), labels = paste(seq(0.1,1,length.out = 10))) +
#       labs(title = paste("Target bi exp vs ", mat_names[j], ", 70%", "N=",sum(conn_target_agglo_sum[[j]][,"tofrom_glu"]),sep = " ")) +
#       geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2)
#   } else {
#     plvl[[j]] <- ggplot(expIpsi_df, aes(x, y, z = z))+
#       geom_raster(data = expIpsi_df, aes(x,y,fill = equalSpace), interpolate = F) +
#       scale_fill_manual(values = pal_tar_ex, guide_legend("dF/F"), labels = paste(seq(0.1,1,length.out = 10))) +
#       labs(title = paste("Target ipsi exp vs", mat_names[j], ", 70%", "N=",sum(conn_target_agglo_sum[[j]][,"tofrom_glu"]),sep = " ")) +
#       geom_contour(data = simdata_df[[j]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2)
#   }
# 
#   plvl[[j]] <- plvl[[j]] +
#     geom_segment(aes(x = -9, y = 112, xend = 0, yend = 112), size=2, lineend = "round") +
#     annotate("text", x = 0, y = 116, label = "9°") +
#     geom_segment(aes(x = -9, y = 112, xend = -9, yend = 103), size=2, lineend = "round") +
#     annotate("text", x = -14.5, y = 102, label = "9°") +
#     # geom_segment(aes(x = 16, y = 49.5, xend = 16, yend = 112.5), size=2, lineend = "round") +
#     geom_segment(aes(x = -9, y = 90-lefty, xend = 16, yend = 90), size=2, lineend = "round") +
#     geom_segment(aes(x = 16, y = 90, xend = 117, yend = 90-righty), size=2, lineend = "round") +
#     # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
#     scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
#     scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
#     geom_segment(aes(x = 16, y = 20, xend = 16, yend = 160), size = 2, lineend = "round") +
#     # geom_segment(aes(x = 0, y = 90, xend = 180, yend = 90), size = 2) +
#     coord_fixed(ratio = 1)
# }
# 
# windows(record = F, width = 8, height = 8)
# plvl[[1]]
# windows(record = F, width = 8, height = 8)
# plvl[[2]]



# target neurons in dolphin glomerulus ----------------------------------------------------------------------------------------------

N_gp <- 11

conn_LC6_tar <- list()
LC6_wt <- list() #neuron skid with synapse weight in each division for each target
for (j in 1:length(neu_target)) {
  tmp <- matrix(ncol = 4)
  for (k in 1:length(neu)) {
    tmp2 <- matrix(ncol = 3)
    conn_tofrom <- catmaid_get_connectors_between(pre_skids = neu[[k]]$skid, post_skids = neu_target[[j]]$skid)
    if (!is.null(conn_tofrom)) {
      tmp2 <- rbind(tmp2, data.matrix(conn_tofrom[, c("connector_x", "connector_y", "connector_z")]))
    }
    if (dim(tmp2)[1] == 2) {
      tmp <- rbind(tmp, c(rep(k,dim(tmp2)[1]-1), tmp2[-1,]))
    } else {
      tmp <- rbind(tmp, cbind(rep(k,dim(tmp2)[1]-1), tmp2[-1,]))
    }
  }
  tmp <- tmp[-1,] %>%
    as_tibble() %>%
    mutate(gp_x = 0)
  colnames(tmp) <- c("ID", "x","y","z","gp_x")
  
  for (k in 1:(length(glu_div)-1)) {
    tmp %<>% as_tibble() %>%
      mutate(gp_x = gp_x + (x>glu_div[k]))
  }
  
  ID_wt <- list() #neuron skid with synapse weight in each division
  for (k in 1:length(glu_div)) {
    df_tmp <- tmp %>%
      filter(gp_x == k) %>%
      transmute(ID) %>%
      data.frame()
    mat_tmp <- matrix(ncol = 2, nrow = length(neu))
    for (m in 1:length(neu)) {
      mat_tmp[m,] <- c(m, sum(df_tmp$ID %in% m))  
    }
    ID_wt[[k]] <- mat_tmp
  }
  
  conn_LC6_tar[[j]] <- as.data.frame(tmp)
  LC6_wt[[j]] <- ID_wt
}
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


# Figure S8B, bar plot
num_tot <- colSums(LC6_ipsi_dol)
for (j in 1:5) {
  LC6_ipsi_dol[,j] <- LC6_ipsi_dol[,j]/num_tot[j]
}
dev.new()
barplot(LC6_ipsi_dol[,c(3,2,4,5,1)], main = 'LC6 to ipsi in dolphin compartments', ylim = c(-0.1, 1.2), col = dolphin_col )
# legend(x = 5, y = 800, legend = c(paste(seq(1,10))), cex = 1.3, fill = dolphin_col, horiz = F)
text(x = seq(0.7,5.5, length.out = 5), y = rep(1.1,5), labels = num_tot[c(3,2,4,5,1)])
text(x = seq(0.7,5.5, length.out = 5), y = rep(-0.05,5), labels = ipsi_skid[c(3,2,4,5,1)])


# separate targets
LC6_wt_biL <- list() # first 2 are biL, 2 biR, last 5 are ipsi
LC6_wt_biR <- list()
LC6_wt_bi <- list()
LC6_wt_ipsi <- list()
for (j in 1:(length(glu_div)-1)) {
  LC6_wt_biL[[j]] <- LC6_wt[[1]][[j]] + LC6_wt[[2]][[j]]
  LC6_wt_biR[[j]] <- LC6_wt[[3]][[j]] + LC6_wt[[4]][[j]]
  LC6_wt_bi[[j]] <- LC6_wt[[1]][[j]] + LC6_wt[[2]][[j]] + LC6_wt[[3]][[j]] + LC6_wt[[4]][[j]]
  LC6_wt_ipsi[[j]] <- LC6_wt[[5]][[j]] + LC6_wt[[6]][[j]] + LC6_wt[[7]][[j]] + LC6_wt[[8]][[j]] + LC6_wt[[9]][[j]]
}

# num of synapse per division
N_syn <- matrix(nrow = length(glu_div)-1, ncol = 3)
for (j in 1:(length(glu_div)-1)) {
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
for (j in 1:(length(glu_div)-1)) {
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
pglu <- ggplot(conn_LC6_tar_bi) +
  geom_point(aes(x = x, y = y, colour = gp_x), shape = 16) 
  
for (j in 1:(length(glu_div)-1)) { #add ahulls
  pglu <- pglu + 
    geom_segment(data = ii_ahull_ls[[j]], aes(x = x1, y = y1, xend = x2, yend = y2))
}
pglu <- pglu +
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
  annotate("text", x = 385000, y = 204000, label = "10 µm")+
  labs(title = "Synapses distribution in glumerulus")
for (j in 1:(length(glu_div)-2)) {
  glu_bd <- data.frame(x1 = range_syn[j,2], x2 = range_syn[j,2], y1 = range_syn[j,3], y2 = range_syn[j,4])
  # pglu <- pglu + geom_segment(data = glu_bd, aes(x = x1, y = y1, xend = x2, yend = y2))
}
for (j in 1:(length(glu_div)-1)) {
  pglu <- pglu + 
    # annotate("text", x = com_syn[j,1], y = com_syn[j,2], label = paste(sprintf("%.1f %%", 100*N_syn[j,2]/sum(N_syn[,2]))))
    annotate("text", x = com_syn[j,1], y = com_syn[j,2], label = paste(sprintf("%.1f ", 100*N_syn[j,2]/sum(N_syn[,2]))))
}

# PLOT, glu with 3 groups
windows(record = F, width = 16, height = 10)
ggdraw() +
  draw_plot(pglu + theme(legend.justification = "bottom"), 0, 0.51, 1, 0.5) +
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
for (j in 1:(length(glu_div)-1)) {
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
# for (j in 1:(length(glu_div)-1)) {
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
pglu <- ggplot(conn_LC6_tar_ipsi) +
  geom_point(aes(x = x, y = y, colour = gp_x), shape = 16) 
  # geom_segment(data = ii_ahull, aes(x = x1, y = y1, xend = x2, yend = y2)) +
for (j in 1:(length(glu_div)-1)) { #add ahulls
  pglu <- pglu + 
    geom_segment(data = ii_ahull_ls[[j]], aes(x = x1, y = y1, xend = x2, yend = y2))
}
pglu <- pglu +
  scale_colour_manual(values = dolphin_col) +
  guides(col = guide_legend(title = "groupings", override.aes = list(shape = rep(".", 10)))) +
  theme_void() +
  theme(legend.position = "bottom") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  # scale_y_reverse() +
  ylim(205000, 180000)+
  coord_fixed(ratio = 1) +
  geom_segment(aes(x = 380000, y = 205000, xend = 390000, yend = 205000), size=2, lineend = "round") +
  annotate("text", x = 385000, y = 204000, label = "10 µm")+
  labs(title = "Synapses distribution in glumerulus")
for (j in 1:(length(glu_div)-2)) {
  glu_bd <- data.frame(x1 = range_syn[j,2], x2 = range_syn[j,2], y1 = range_syn[j,3], y2 = range_syn[j,4])
}
for (j in 1:(length(glu_div)-1)) {
  pglu <- pglu + 
    annotate("text", x = com_syn[j,1], y = com_syn[j,2], label = paste(sprintf("%.1f ", 100*N_syn[j,3]/sum(N_syn[,3]))))
}

# PLOT, glo with 3 groups
windows(record = F, width = 16, height = 10)
ggdraw() +
  draw_plot(pglu + theme(legend.justification = "bottom"), 0, 0.51, 1, 0.5) +
  draw_plot(gg_cont + theme(axis.title = element_text()), 0, 0, 0.5, 0.5) 



# Figure 6A, all targets with volume ref

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


# -- cp LM
nopen3d()
par3d('windowRect' = c(100,100,1100,1100))
# shade3d(v14, alpha=0.05)
shade3d(glo_vol, col = 'grey90', alpha = 0.3)
ipsi_pal <- brewer.pal(6, "RdYlBu")[c(1,2,3,5,6)]
tar <- neu_ipsi[[2]]
plot3d(tar, col = ipsi_pal[1])
points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = ipsi_pal[1], size = 20)
tar <- neu_ipsi[[3]]
plot3d(tar, col = ipsi_pal[5])
points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = ipsi_pal[5], size = 20)
rgl.viewpoint(userMatrix = rotationMatrix(1*90/180*pi+pi/2,1,0,0) %*% rotationMatrix(0,0,1,0), zoom = 1)

# Figure S7, individual targets RF
for (j in 1:2) {
  nopen3d()
  par3d('windowRect' = c(100,100,1100,1100))
  rgl.viewpoint(userMatrix = rotationMatrix(1*90/180*pi+pi/2,1,0,0) %*% rotationMatrix(0,0,1,0), zoom = 1)
  # shade3d(v14, alpha=0.05)
  shade3d(glo_vol, col = 'grey90', alpha = 0.3)
  tar <- neu_biR[[j]]
  # tar <- neu_biL[[j]]
  # tar <- neu_ipsi[[j]]
  plot3d(tar, col = ipsi_pal[2*j], lwd = 3)
  points3d(xyzmatrix(tar$d[match(tar$tags$soma, tar$d$PointNo), ]), col = ipsi_pal[2*j], size = 20)
  title3d(paste("biR_", biR_skid[j]))
  rgl.snapshot(filename = paste("biR_",  j, ".png", sep = ''),fmt = "png")
}

# rgl.snapshot(filename = "ipsi_x2.png",fmt = "png")


# data from Mai ---------------------------------------------------------------------------------------------------

library(cowplot)
library(R.matlab)
expBi3 <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss825_bi_RF_indiv_tseries.mat") #bi
expBi3 <- expBi3[[1]]
expBi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss825_mean_RFmap_n=4.mat") #bi
# expBi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss825_mean_AUC_RFmap_n=4.mat") #bi
expBi <- as.matrix(expBi[[1]])


expIpsi3 <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss2036_ipsi_RF_indiv_tseries.mat") #ipsi
expIpsi3 <- expIpsi3[[1]]
expIpsi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss2036_mean_RFmap_n=7.mat") #ipsi
# expIpsi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss2036_mean_AUC_RFmap_n=7.mat") #ipsi
expIpsi <- as.matrix(expIpsi[[1]])

# indexing
ind_mai <- c(t(matrix(seq(1:98), byrow = F, ncol = 14)))
tar_pal <- brewer.pal(4,"RdYlBu")[c(1,3,4,2)]


# x2 <- seq(-18 - 2.25, 99 - 2.25, by = 9) # -18 to 99, with sim data
# y2 <- seq(54, 108,by = 9)
# xygrid2 <- expand.grid(x2, y2) # looming center position 

# loom position from Matt
loom_theta_mat <- read.csv('loom_center_theta.txt', sep = ',', header = F) / pi * 180
loom_theta_mat <- 180 - loom_theta_mat
loom_theta_mat <- loom_theta_mat[seq(7,1),]
loom_phi_mat <- read.csv('loom_center_phi.txt', sep = ',', header = F) / pi * 180
loom_phi_mat <- loom_phi_mat[seq(7,1),]

loom_theta <- melt(t(as.matrix(loom_theta_mat)))$value 
loom_phi <- melt(t(as.matrix(loom_phi_mat)))$value 

xygrid2 <- cbind(loom_phi, loom_theta) 

# shim
shim_xy <- read.csv('shim.csv', sep = ',', header = F)
shim_xy <- as.matrix(shim_xy)
shim_xy <- t(shim_xy) / pi * 180
shim_xy[,2] <- - shim_xy[,2] + 90
shim_xy <- shim_xy[shim_xy[,1] > min(loom_phi)-4.5 & shim_xy[,1] < max(loom_phi)+4.5 & shim_xy[,2] > min(loom_theta)-4.5 & shim_xy[,2] < max(loom_theta)+4.5, ]
# head(shim_xy,1); tail(shim_xy,1); max(shim_xy[,1]); min(shim_xy[,2])
ptadd <- c(max(shim_xy[,1]), min(shim_xy[,2]))
shim_xy <- rbind(shim_xy, ptadd)
shim_xy <- as.data.frame(shim_xy)
colnames(shim_xy) <- c('x','y')


# construct data frame
expBi_df <- data.frame(xygrid2, as.vector(t(expBi)))
colnames(expBi_df) <- c("x","y","z")
expIpsi_df <- data.frame(xygrid2, as.vector(t(expIpsi)))
colnames(expIpsi_df) <- c("x","y","z")

# # slope of LED panel, cf. S5B
# slope_ex <- (90-83.28)/(117-16)
# lefty <- slope_ex*abs(x2[1]) 
# righty <- slope_ex*abs(tail(x2,1))


# not used --------------------------------------------------------------------------------------------------------

# # plot indiv traces
# mai <- expBi3
# plt_exp <- list()
# for (j in 1:dim(mai)[3]) {
#   amp_max <- quantile(na.omit(c(mai[,,j])), probs = c(0.98)) # normalize to 98%
#   for (k in 1:dim(mai)[1]) {
#     # amp_max <- quantile(na.omit(c(exp_raw[[j]][,,k])), probs = c(0.98)) # normalize to 98% 
#     # exp_raw[[j]][,,k] <- exp_raw[[j]][,,k] / amp_max
#     df <- data.frame(x = seq(1,dim(mai)[2]), y = mai[ind_mai[k],,j]/amp_max)
#     plt_exp[[(j-1)*dim(mai)[1] + k]] <- ggplot(df, aes(x,y)) + 
#       geom_point() 
#       # ylim(0,1)
#       # theme_void() 
#   }
#   
# }
# dev.new()
# plot_grid(plotlist = plt_exp, ncol = 14, labels = seq(1,98))
# 
# # plot mean
# mai_mean <- rowMeans(mai, dims = 2, na.rm = T)
# amp_max <- quantile(na.omit(c(mai_mean)), probs = c(0.98)) # normalize to 98% 
# for (k in 1:dim(mai_mean)[1]) {
#   df <- data.frame(x = seq(1,dim(mai_mean)[2]), y = mai_mean[ind_mai[k],]/amp_max)
#   plt_exp[[k]] <- ggplot(df, aes(x,y)) + 
#     geom_point() +
#     ylim(0,1)
#   # theme_void() 
# }
# dev.new()
# plot_grid(plotlist = plt_exp, ncol = 14, labels = seq(1,98))

# stim_end <- 200 # loom during the first 150 points, 4 sec
# 
# # normalize, mean aross flies, then normalize again
# exp_raw <- list(expBi3, expIpsi3)
# pval <- list()
# j <- 1
#   # normalize within fly
#   for (k in 1:dim(exp_raw[[j]])[3]) {
#     amp_max <- quantile(na.omit(c(exp_raw[[j]][,,k])), probs = c(0.98)) # normalize to 98% 
#     exp_raw[[j]][,,k] <- exp_raw[[j]][,,k] / amp_max
#   }
#   
# # t-test
# mu_test <- sum(exp_raw[[j]][ind_mai[c(12,13,14,27,28,42,56)], 1:stim_end,]) / 7 / dim(exp_raw[[j]])[3]
# pval_tmp <- c()
# for (k in 1:dim(exp_raw[[j]])[1]) {
#   pval_tmp[k] <- t.test(colSums(exp_raw[[j]][k, 1:stim_end,]), mu = mu_test)$p.value
# }
# pval[[j]] <- pval_tmp
# 
#   # mean and normalize across flies
#   mai_mean <- rowMeans(mai, dims = 2, na.rm = T)
#   amp_max <- quantile(na.omit(c(mai_mean)), probs = c(0.98)) # normalize to 98% 
#   
#   max(na.omit(c(mai_mean))) / amp_max
#   min(na.omit(c(mai_mean))) / amp_max


# # -- coutour from data
# # x2 <- seq(0-4.25,117-4.25,by = 9) # -18 to 99, with sim data
# # y2 <- seq(54, 108,by = 9)
# x3 <- seq(-18,99,by = 9) # -18 to 99, with sim data
# y3 <- seq(54, 108,by = 9)
# xygrid3 <- expand.grid(x3, y3)
# 
# exp_df <- data.frame(xygrid3, as.vector(t(expBi)))
# exp_df <- data.frame(xygrid3, as.vector(t(expIpsi)))
# colnames(exp_df) <- c("x","y","z")
# 
# exp_loess <- loess (z ~ x * y, exp_df, degree = 2, span = 0.3)
# x4 <- seq(-18, 99, by = 1)
# y4 <- seq(54, 108, by = 1)
# xygrid4 <- expand.grid(x = x4, y = y4)
# loess_fit <- predict (exp_loess, xygrid4, se = T)
# exp_interp <- cbind(xygrid4, melt(loess_fit$fit)$value)
# colnames(exp_interp) <- c('x','y','z')
# 
# dev.new()
# image(x= x4, y= y4, z = loess_fit$fit,  asp = 1)
# points(exp_df)
# 
# cut <- quantile(exp_interp$z, probs = c(0.7))
# # cut <- 0.7 * max(exp_interp$z)
# exp_interp_cut <- exp_interp[exp_interp$z > cut, ]
# 
# dim(exp_interp_cut)[1]



#  Figure 6E -- time series mean + EM contour  --------------------------------------------------------------------

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
LC6_tar_median <- list()
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
    A <- conn_target_agglo_sum[[j]][k,"tofrom_glu"]
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
  
  plvl[[j]] <- plvl[[j]] + labs(title = paste("target", mat_names[j], ", 70%", "N=",sum(conn_target_agglo_sum[[j]][,"tofrom_glu"]),sep = " ")) +
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


# #-- mean exp contour with EM bi - ipsi heat map
#   ggplot() +
#   geom_polygon(data = shim_xy, aes(x,y), fill = 'black', alpha = 0.3) +
#   geom_raster(data = simdata_df_diff, aes(x, y, z = z, fill = equalSpace), interpolate = F) +
#   scale_fill_manual(values = pal_tar) +
#   geom_path(data = xy_bd_chull_df, aes(x = phi_deg, y = theta_deg), colour = 'black',  inherit.aes = FALSE) +
#   geom_contour(data = exp_interp, aes(x,y,z=z), breaks = c(0.61), color = "red", alpha = 1, lwd = 2) +
#   geom_segment(aes(x = -12, y = 135, xend = -3, yend = 135), size=2, lineend = "round") +
#   annotate("text", x = -3, y = 139, label = "9°") +
#   geom_segment(aes(x = -12, y = 135, xend = -12, yend = 124), size=2, lineend = "round") +
#   annotate("text", x = -17.5, y = 124, label = "9°") +
#   # geom_segment(aes(x = -9, y = 90-lefty, xend = 16, yend = 90), size=2, lineend = "round") +
#   # geom_segment(aes(x = 16, y = 90, xend = 117, yend = 90-righty), size=2, lineend = "round") +
#   geom_segment(aes(x = 0, y = 50, xend = 0, yend = 120), size = 2) + # range(loom_phi)
#   geom_segment(aes(x = 90, y = 50, xend = 90, yend = 120), size = 2) +
#   geom_segment(aes(x = -25, y = 90, xend = 115, yend = 90), size = 2) +
#   # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
#   scale_x_continuous(breaks = x_tick, labels = paste(c(0, 45, 90, 135, 180))) +
#   scale_y_reverse(breaks = y_tick, labels = paste(c(90, 45, 0, -45 ,-90))) +
#   theme_void() +
#   coord_fixed(ratio = 1)

# # areas of RF
# dim(ct_cut_pt[[1]])[1] / sum(ii_inpoly)
# 
# dim(ct_cut_pt[[2]])[1] / sum(ii_inpoly)


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
  plvl_ol[[j]] <- plvl_ol[[j]] + labs(title = paste("Target ipsi exp vs", mat_names[j], ", 70%", "N=",sum(conn_target_agglo_sum[[j]][,"tofrom_glu"]),sep = " ")) +
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
  
  
  
# Figure 3B and S5C
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

# Figure 3D alt
windows(record = F, width = 8, height = 8)
plvl_ol[[1]] 
windows(record = F, width = 8, height = 8)
plvl_ol[[2]] 

# Figure SX, difference EM, bi - ipsi, indiv exp
windows(record = F, width = 8, height = 8)
plvl_ol_em[[1]] 
windows(record = F, width = 8, height = 8)
plvl_ol_em[[2]] 

# Figure 6D
# contour map + exp contour
windows(record = F, width = 8, height = 8)
plvl_ol[[1]] + geom_contour(data = simdata_df[[1]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2)

windows(record = F, width = 8, height = 8)
plvl_ol[[2]] + geom_contour(data = simdata_df[[2]], aes(x,y,z=z), breaks = c(0.71), color = "blue", alpha = 1, lwd = 2)



# -- exp and EM overlap
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

library(Hmisc)
windows(record = F, width = 8, height = 8)
ggplot(ol_ratio, aes(x = pair, group = pair, y = ratio)) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize = 1) +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5) +
  scale_x_continuous(breaks=c(1,2,3,4), labels=c("bi-biEM", "ipsi-biEM", "bi-ipsiEM", "ipsi-ipsiEM"))+
  theme_bw()

# ggplot(ol_ratio, aes(x = pair, group = pair, y = ratio)) + geom_boxplot()+
#   scale_x_continuous(breaks=c(1,2,3,4), labels=c("bi-biEM", "ipsi-biEM", "bi-ipsiEM", "ipsi-ipsiEM")) +
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) 
  


# -- difference in EM contour
# simdata_df_diff <- simdata_df[[1]][,1:3]
# simdata_df_diff$z <- simdata_df_diff$z - simdata_df[[2]]$z
# simdata_df_diff$z <- simdata_df_diff$z / max(abs(simdata_df_diff$z))
# # simdata_df_diff$equalSpace <- cut(simdata_df_diff$z, seq(0,max(abs(simdata_df_diff$z)),length.out = n_lvl))
# # simdata_df_diff$equalSpace <- cut(simdata_df_diff$z, seq(min(simdata_df_diff$z),max(simdata_df_diff$z),length.out = n_lvl))
# simdata_df_diff$equalSpace <- cut(simdata_df_diff$z, seq(-1,1,length.out = n_lvl))

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
  geom_contour(data = bi_ct[[1]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[1], alpha = 1, lwd = 2) +
  geom_contour(data = bi_ct[[2]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[3], alpha = 1, lwd = 2) +
  geom_contour(data = bi_ct[[3]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[5], alpha = 1, lwd = 2) +
  geom_contour(data = bi_ct[[4]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[7], alpha = 1, lwd = 2) +
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
  geom_contour(data = ipsi_ct[[1]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[1], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[2]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[2], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[3]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[3], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[4]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[4], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[5]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[5], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[6]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[6], alpha = 1, lwd = 2) +
  geom_contour(data = ipsi_ct[[7]], aes(x,y,z=z), breaks = c(0.61 * max(exp_interp[,3])), color = ct_col[7], alpha = 1, lwd = 2) +
  coord_fixed(ratio = 1)

dev.off()


# user ------------------------------------------------------------------------------------------------------------

user_stat <- catmaid_get_contributor_stats(c(LC6_skid, ipsi_skid, biL_skid, biR_skid))

user_stat$node_contributors[order(user_stat$node_contributors$n),]
user_stat$pre_contributors[order(user_stat$pre_contributors$n), ]
user_stat$post_contributors[order(user_stat$post_contributors$n), ]
user_stat$review_contributors[order(user_stat$review_contributors$n), ]

udf <- catmaid_get_user_list()

# make a df of contributions 
jd <- full_join(user_stat$post_contributors, user_stat$pre_contributors, by = "id")
jd <- full_join(jd, user_stat$review_contributors, by = 'id')
jd <- full_join(jd, user_stat$node_contributors, by = 'id')
jd[is.na(jd)] <- 0
jd <- left_join(jd, udf) %>%
  # mutate(n = n.x + n.y + n.x.x + n.y.y) %>%
  mutate(n = n.x + n.y ) %>%
  select(id, n, full_name) %>%
  arrange(desc(n)) 

9


# Review stuff below------------------------------------------------------------------------------------------------


# approx retinotopy along axon-----------------------------------------------------------------------------------------------

cross3D <- function(a, b){
  prod <- matrix(c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1]), ncol = 1)
  return(prod)
}

Nn <- length(xyz_com)

xyz_com_m <- matrix(unlist(xyz_com), ncol = 3, byrow = T)
ind_1 <- seq(1,Nn)[xyz_com_m[,2] < 250000]
ind_2 <- seq(1,Nn)[xyz_com_m[,2] >= 250000]

nopen3d()
plot3d(LC6[ind_1], col = 'green')
plot3d(LC6[ind_2], col = 'blue')


which.min(xyz_com_m[,2]) #top = 30
which.max(xyz_com_m[,2]) #bottom = 61
nt <- 30
nb <- 12
# nt <- 34
# nb <- 3

# LC6 <- neu

nref <- LC6[[nb]] #origin
# o_xyz_m <- xyzmatrix(nref$d)
nrefrs <- resample(nref, stepsize = 8000)
o_xyz_m_rs <- xyzmatrix(nrefrs$d)

nref_z <- LC6[[nt]] # +z axis
z_xyz_m <- xyzmatrix(nref_z$d)

# nopen3d()
# plot3d(nref)
# points3d(xyzmatrix(nrefrs$d), size = 10)
# identify3d(xyzmatrix(nrefrs$d))

# which(sapply(nrefrs$SegList, function(x) 35 %in% x))
# nrefrs$SegList[[2]]

# ref_ind <- seq(35,55, length.out = 3) #use these nodes as ref points
# N_ind <- length(ref_ind)-1

# plot ref position along neu
nopen3d()
plot3d(nrefrs, WithNodes = F)
points3d(xyzmatrix(nrefrs$d), size =5)

# identify3d(xyzmatrix(nrefrs$d))

ref_ind <- c(32, 42, 84, 85) #only use first 3, last is to defind +x
N_ind <- length(ref_ind) - 1

nopen3d()
plot3d(nrefrs, WithNodes = F)
points3d(xyzmatrix(nrefrs$d[ref_ind,]), col = 'red',size = 10)


yz_ls <- list()
ref_d <- c()
for (j in 1:N_ind ) {
  xp <- o_xyz_m_rs[ref_ind[j+1],] - o_xyz_m_rs[ref_ind[j],] # +x axis along cable
  o_xyz <- o_xyz_m_rs[ref_ind[j],]
  z_ind <- which.min(rowSums(sweep(z_xyz_m, MARGIN = 2, STATS = o_xyz, '-')^2))
  z_xyz <- z_xyz_m[z_ind,] # +z, no.30
  zp <- z_xyz - o_xyz
  ref_d <- c(ref_d, sqrt(sum(zp^2)))
  
  
  pt_yz <- matrix(ncol = 2, nrow = length(LC6))
  for (k in 1:nrow(pt_yz)) {
    xyz_m <- xyzmatrix(LC6[[k]]$d)
    ind <- which.min(rowSums(sweep(xyz_m, MARGIN = 2, STATS = o_xyz, '-')^2))
    xyz <- xyz_m[ind,]
    v <- xyz - o_xyz
    ang <- acos(v %*% zp / sqrt(sum(v^2)) / sqrt(sum(zp^2)))
    ang <- ang * sign(xp %*% cross3D(v, zp) )
    mag <- sqrt(sum(v^2)) / sqrt(sum(zp^2)) #normalize to +z
    pt_yz[k,] <- c(sin(ang)*mag, cos(ang)*mag)
  }
  pt_yz[nt,] <- c(0,1)
  pt_yz[nb,] <- c(0,0)
  yz_ls[[j]] <- pt_yz
}

ref_d <- c(sqrt(sum((xyz_com[[nt]]-xyz_com[[nb]])^2)),ref_d)


# com on eye
yz_eye <- matrix(unlist(xy_com), ncol = 4, byrow = T)[,3:4]
zp <- yz_eye[nt,] - yz_eye[nb,]
angR <- pi - atan(zp[1]/zp[2])
yz_eye <- t( matrix(c(cos(angR), -sin(angR), sin(angR), cos(angR)), ncol = 2 ) %*%
  t(sweep(yz_eye, 2, yz_eye[nb,], '-')) )
yz_eye <- yz_eye / yz_eye[nt,2] #normalize
yz_eye <- data.frame(yz_eye)
colnames(yz_eye) <- c('y','z')

# make df for plotting
df <- matrix(ncol = 3, nrow = 0)
for (j in 1:length(LC6) ) {
  df <- rbind(df, as.numeric(c(yz_eye[j,], j)) )
  for (k in 1:(N_ind)) {
    df <- rbind(df, c(yz_ls[[k]][j,], j) )
  }
}
df <- as.data.frame(df)
colnames(df) <- c('y','z','n')
df$n <- factor(df$n)

df_cable <- df

# color palette
getPalette <- colorRampPalette(brewer.pal(9, "RdYlBu"))
cable_col <- getPalette(length(LC6))

# make plot
plt <- ggplot()
# for (j in 1:length(LC6)) {
for (j in c(nt,nb) ){
  df_plt <- df_cable[((N_ind+1)*(j-1)+1):((N_ind+1)*j),]
  df_plt <- cbind(df_plt, seq(1,N_ind+1)*0.1+1)
  colnames(df_plt)[4] <- 'N_ind'
  plt <- plt + 
    # geom_point(data = df_cable, aes(y, z,colour = n, size = N_ind) ) +
    # ggplot() +
    geom_point(data = df_plt, aes(y, z) ) +
    geom_path(data = df_plt, aes(x=y,y=z), colour = cable_col[j], lwd =2) +
    geom_point(data = yz_eye[j,], aes(x=y,y=z), size = 6, shape = 18)
}

dev.new()
plt + xlim(-1.5,1.5) +
  ylim(-1.5,1.5)



# -- repeat for glo spine but with glu_div

# first get the backbone within glo, and resample

start_xyz <- xyzmatrix(nrefrs$d[ref_ind[2],]) #starting node for all neu, the 2rd position

LC6_sp_glo_rs <- LC6 # initiate

for (ni in 1:length(LC6)) {
  tar <- LC6[[ni]]

  tar_xyz <- xyzmatrix(tar$d)
  ind <- which.min(rowSums(sweep(tar_xyz, MARGIN = 2, STATS = start_xyz, '-')^2))

  targ <- as.ngraph(tar)
  distal_points <- igraph::graph.dfs(targ, root=ind, unreachable=FALSE, neimode='out')$order
  tar_distal <- subset(tar, distal_points)

  sp <- spine(tar_distal, UseStartPoint = T)
  sp_glo <- subset(sp, pointsinside(sp, surf=glo.msh, rval='distance') > -5000)
  LC6_sp_glo_rs[[ni]] <- resample(sp_glo, stepsize = 400)
}

# second pick 2 ref
i_AP <- c(4, 10) #long but not i_excl, 4 ante, 10 post
# i_AP <- c(12, 30) # as for tract

nref_sp <- LC6_sp_glo_rs[[i_AP[1]]] #origin
# o_xyz_m <- xyzmatrix(nref_sp$d)
# nref_sp <- subset(nref_sp, pointsinside(nref_sp, surf=glo.msh, rval='distance') > -5000)
nrefrs_sp <- resample(nref_sp, stepsize = 1000)
o_xyz_m_rs_sp <- xyzmatrix(nrefrs_sp$d)

nref_z_sp <- LC6_sp_glo_rs[[i_AP[2]]] # +z axis
z_xyz_m_sp <- xyzmatrix(nref_z_sp$d)

# nopen3d()
# plot3d(nrefrs_sp)
# points3d(xyzmatrix(nrefrs_sp$d), size = 10)
# identify3d(xyzmatrix(nrefrs_sp$d))

# choose ref pt, mid point of glu_div
ref_ind_sp <- c() 
for (j in 1:10) {
  mid_x <- (glu_div[j] + glu_div[j+1])/2
  ii <- which.min((o_xyz_m_rs_sp[,1] - mid_x)^2)
  ref_ind_sp <- c(ref_ind_sp, ii)
}

N_ind_sp <- length(ref_ind_sp)

# # DEBUG
# nopen3d()
# plot3d(nrefrs_sp, WithNodes = F)
# # points3d(xyzmatrix(nrefrs_sp$d[ref_ind_sp,]), col = 'red',size = 10)
# points3d(o_xyz_m_rs_sp, col = 'grey',size = 5)
# points3d(o_xyz_m_rs_sp[ref_ind_sp,], col = 'red',size = 10)
# arrow3d(o_xyz_m_rs_sp[ref_ind_sp[j],], o_xyz_m_rs_sp[ref_ind_sp[j],] + xp*5e3, theta = pi / 9, n = 8, col = "red", type = "rotation")
# points3d((z_xyz_m_sp[z_ind,]), size =5, col='blue')
# arrow3d(o_xyz, o_xyz + zp, theta = pi / 9, n = 8, col = "red", type = "rotation")
# planes3d(xp[1],xp[2],xp[3],-d1, alpha=0.3)
# planes3d(xp[1],xp[2],xp[3],-d2, alpha=0.3)
# points3d(xyz_m[xyz_m_ind,], size =10, col = 'cyan')


yz_ls <- list()
ref_d <- c()
dL <- 500 # half thickness of cross section
for (j in 1:N_ind_sp) {
  # use pc1 to def x-axis
  xp_1 <- o_xyz_m_rs_sp[ref_ind_sp[j],] - o_xyz_m_rs_sp[ref_ind_sp[j]-1,] # +x axis along cable
  xp_pc <- o_xyz_m_rs_sp[(ref_ind_sp[j]-3):(ref_ind_sp[j]+3),] %>%
    prcomp() %>% 
    .[["rotation"]] %>% 
    .[,"PC1"]
  if (xp_1 %*% xp_pc > 0) {
    xp <- xp_pc
  } else {
    xp <- - xp_pc
  }
  
  o_xyz <- o_xyz_m_rs_sp[ref_ind_sp[j],] #origin
  d1 <- (o_xyz + xp*dL) %*% xp #distance to plane
  d2 <- (o_xyz - xp*dL) %*% xp #distance to plane
  
  # z_ind <- which.min(rowSums(sweep(z_xyz_m_sp, MARGIN = 2, STATS = o_xyz, '-')^2))
  # z_xyz <- z_xyz_m_sp[z_ind,] # +z
  
  # use 2 planes
  z_ind <- apply(z_xyz_m_sp, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
  if (sum(z_ind) == 1 ) {
    z_xyz <- z_xyz_m_sp[z_ind,] 
  } else {
    z_xyz <- colMeans(z_xyz_m_sp[z_ind,]) 
  }
  zp <- z_xyz - o_xyz
  zp <- zp - c(zp %*% xp)*xp #+z
  
  ref_d <- c(ref_d, sqrt(sum(zp^2))) #unit length
  
  
  pt_yz <- matrix(NA, ncol = 2, nrow = length(LC6))
  for (k in 1:nrow(pt_yz)) {
    xyz_m <- xyzmatrix(LC6_sp_glo_rs[[k]]$d)
    # xyz_m <- xyz_m[xyz_m[,1] > glu_div[j] & xyz_m[,1] < glu_div[j+1],]
    xyz_m_ind <- apply(xyz_m, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
    if (sum(xyz_m_ind) > 0) {
      if (sum(xyz_m_ind) == 1) { #if single pt
        xyz <- xyz_m[xyz_m_ind,] 
      } else {
        xyz <- colMeans(xyz_m[xyz_m_ind,]) 
      }
      # dd <- sqrt(rowSums(sweep(xyz_m, MARGIN = 2, STATS = o_xyz, '-')^2))
      # ind <- which.min(dd)
      # xyz <- xyz_m[ind,]
      v <- xyz - o_xyz
      v <- v - c(v %*% xp)*xp 
      ang <- acos(v %*% zp / sqrt(sum(v^2)) / sqrt(sum(zp^2)))
      # ang <- ang * sign(v %*% zp)
      ang <- ang * sign(xp %*% cross3D(zp, v) )
      mag <- sqrt(sum(v^2)) / sqrt(sum(zp^2)) #normalize to +z
      pt_yz[k,] <- c(sin(ang)*mag, cos(ang)*mag)
    }
  }
  pt_yz[i_AP[2],] <- c(0,1)
  pt_yz[i_AP[1],] <- c(0,0)
  yz_ls[[j]] <- pt_yz
}

 # add eye
ref_d <- c(sqrt(sum((xyz_com[[i_AP[2]]]-xyz_com[[i_AP[1]]])^2)),ref_d) 

# com on eye
yz_eye <- matrix(unlist(xy_com), ncol = 4, byrow = T)[,3:4]
zp <- yz_eye[i_AP[2],] - yz_eye[i_AP[1],]
angR <- pi - atan(zp[1]/zp[2])
yz_eye <- t( matrix(c(cos(angR), -sin(angR), sin(angR), cos(angR)), ncol = 2 ) %*%
               t(sweep(yz_eye, 2, yz_eye[i_AP[1],], '-')) )
yz_eye <- yz_eye / yz_eye[i_AP[2],2] #normalize
yz_eye <- data.frame(yz_eye)
colnames(yz_eye) <- c('y','z')

# make df for plotting
df <- matrix(ncol = 3, nrow = 0)
for (j in 1:length(LC6) ) {
  df <- rbind(df, as.numeric(c(yz_eye[j,], j)) )
  for (k in 1:N_ind_sp) {
    df <- rbind(df, c(yz_ls[[k]][j,], j) )
  }
}
df <- as.data.frame(df)
colnames(df) <- c('y','z','n')
df$n <- factor(df$n)

df_cable_sp <- df



#  approx retinotopy along axon NEW -------------------------------------------------------------------------------

# combine tract and glo

i_AP <- c(4, 10) #long but not i_excl, 4 ante, 10 post, 4 -> 10 as y-axis

# - first get the backbone, and resample

# nopen3d()
# plot3d(LC6[[i_AP[1]]], WithNodes = F)
# points3d(xyzmatrix(LC6[[i_AP[1]]]$d), size = 10)
# plot3d(LC6)
# identify3d(xyzmatrix(LC6[[i_AP[1]]]$d))

start_xyz <- xyzmatrix(LC6[[i_AP[1]]]$d[4727,]) #starting node for all neu

LC6_sp_rs <- LC6 # initiate

for (ni in 1:length(LC6)) {
  tar <- LC6[[ni]]
  
  tar_xyz <- xyzmatrix(tar$d)
  ind <- which.min(rowSums(sweep(tar_xyz, MARGIN = 2, STATS = start_xyz, '-')^2))
  
  targ <- as.ngraph(tar)
  distal_points <- igraph::graph.dfs(targ, root=ind, unreachable=FALSE, neimode='out')$order
  tar_distal <- subset(tar, distal_points)
  
  sp <- spine(tar_distal, UseStartPoint = T)
  # sp_glo <- subset(sp, pointsinside(sp, surf=glo.msh, rval='distance') > -5000)
  LC6_sp_rs[[ni]] <- resample(sp, stepsize = 400)
}


# - 2 ref set up ref frame
nref_sp <- LC6_sp_rs[[i_AP[1]]] #origin
# o_xyz_m <- xyzmatrix(nref_sp$d)
# nref_sp <- subset(nref_sp, pointsinside(nref_sp, surf=glo.msh, rval='distance') > -5000)
nrefrs_sp <- resample(nref_sp, stepsize = 1000)
o_xyz_m_rs_sp <- xyzmatrix(nrefrs_sp$d)

nref_y_sp <- LC6_sp_rs[[i_AP[2]]] # +y axis, NEW
y_xyz_m_sp <- xyzmatrix(nref_y_sp$d)

nopen3d()
plot3d(nrefrs_sp)
points3d(xyzmatrix(nrefrs_sp$d), size = 10)
plot3d(LC6_sp_rs)
shade3d(glo.msh, col='grey',alpha = 0.2)
# identify3d(xyzmatrix(nrefrs_sp$d))

# 7, 99 (before split), 175(before glo), 246 (afte split)

# choose ref pt, mid point of glu_div
ref_ind_sp <- c() 
for (j in 1:10) {
  mid_x <- (glu_div[j] + glu_div[j+1])/2
  ii <- which.min((o_xyz_m_rs_sp[175:323,1] - mid_x)^2) + 174 # 175:323 use pts within glo
  ref_ind_sp <- c(ref_ind_sp, ii)
}
# ref_ind_sp <- c(7, 99, ref_ind_sp[c(4,6)]) #all 65
ref_ind_sp <- c(7, 175, ref_ind_sp[c(9, 5)]) #better positions

N_ind_sp <- length(ref_ind_sp)

# DEBUG
nopen3d()
plot3d(nrefrs_sp, WithNodes = F)
# points3d(xyzmatrix(nrefrs_sp$d[ref_ind_sp,]), col = 'red',size = 10)
points3d(o_xyz_m_rs_sp, col = 'grey',size = 5)
points3d(o_xyz_m_rs_sp[ref_ind_sp,], col = 'red',size = 10)
arrow3d(o_xyz_m_rs_sp[ref_ind_sp[j],], o_xyz_m_rs_sp[ref_ind_sp[j],] + xp*1e3, theta = pi / 9, n = 8, col = "red", type = "rotation")
points3d((y_xyz_m_sp[y_ind,]), size =5, col='blue')
arrow3d(o_xyz, o_xyz + yp, theta = pi / 18, n = 8, col = "red", type = "rotation")
planes3d(xp[1],xp[2],xp[3],-d1, alpha=0.3)
planes3d(xp[1],xp[2],xp[3],-d2, alpha=0.3)
points3d(xyz_m[xyz_m_ind,], size =10, col = 'cyan')


yz_ls <- list()
ref_d <- c()
dL <- 500 # half thickness of cross section
for (j in 1:N_ind_sp) {
  # use pc1 to def x-axis
  xp_1 <- o_xyz_m_rs_sp[ref_ind_sp[j],] - o_xyz_m_rs_sp[ref_ind_sp[j]-1,] # +x axis along cable
  xp_pc <- o_xyz_m_rs_sp[(ref_ind_sp[j]-3):(ref_ind_sp[j]+3),] %>%
    prcomp() %>% 
    .[["rotation"]] %>% 
    .[,"PC1"]
  if (xp_1 %*% xp_pc > 0) {
    xp <- xp_pc
  } else {
    xp <- - xp_pc
  }
  
  o_xyz <- o_xyz_m_rs_sp[ref_ind_sp[j],] #origin
  d1 <- (o_xyz + xp*dL) %*% xp #distance to plane
  d2 <- (o_xyz - xp*dL) %*% xp #distance to plane
  
  # z_ind <- which.min(rowSums(sweep(y_xyz_m_sp, MARGIN = 2, STATS = o_xyz, '-')^2))
  # z_xyz <- y_xyz_m_sp[z_ind,] # +z
  
  # use 2 planes
  y_ind_1 <- apply(y_xyz_m_sp, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
  y_ind_2 <- sqrt(rowSums(sweep(y_xyz_m_sp, MARGIN = 2, STATS = o_xyz, '-')^2)) < 20000
  y_ind <- y_ind_1 & y_ind_2
  if (sum(y_ind) == 1 ) {
    y_xyz <- y_xyz_m_sp[y_ind,] 
  } else {
    y_xyz <- colMeans(y_xyz_m_sp[y_ind,]) 
  }
  yp <- y_xyz - o_xyz
  yp <- yp - c(yp %*% xp)*xp # +y
  
  ref_d <- c(ref_d, sqrt(sum(yp^2))) #unit length
  
  pt_yz <- matrix(NA, ncol = 2, nrow = length(LC6))
  for (k in 1:nrow(pt_yz)) {
    xyz_m <- xyzmatrix(LC6_sp_rs[[k]]$d)
    # xyz_m <- xyz_m[xyz_m[,1] > glu_div[j] & xyz_m[,1] < glu_div[j+1],]
    xyz_m_ind_1 <- apply(xyz_m, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
    xyz_m_ind_2 <- sqrt(rowSums(sweep(xyz_m, MARGIN = 2, STATS = o_xyz, '-')^2)) < 10000
    xyz_m_ind <- xyz_m_ind_1 & xyz_m_ind_2 # within 2 planes and 20 um within ref pt
    if (sum(xyz_m_ind) > 0) {
      if (sum(xyz_m_ind) == 1) { #if single pt
        xyz <- xyz_m[xyz_m_ind,] 
      } else {
        xyz <- colMeans(xyz_m[xyz_m_ind,]) 
      }
      v <- xyz - o_xyz
      v <- v - c(v %*% xp)*xp 
      ang <- acos(v %*% yp / sqrt(sum(v^2)) / sqrt(sum(yp^2)))
      # ang <- ang * sign(v %*% yp)
      ang <- ang * sign(xp %*% cross3D(yp, v) )
      mag <- sqrt(sum(v^2)) / sqrt(sum(yp^2)) #normalize to 1
      pt_yz[k,] <- c(cos(ang)*mag, sin(ang)*mag)
    }
  }
  pt_yz[i_AP[2],] <- c(1,0)
  pt_yz[i_AP[1],] <- c(0,0)
  yz_ls[[j]] <- pt_yz
}

# add eye
ref_d <- c(sqrt(sum((xyz_com[[i_AP[2]]]-xyz_com[[i_AP[1]]])^2)),ref_d) 

# com on eye
yz_eye <- matrix(unlist(xy_com), ncol = 4, byrow = T)[,c(4,3)]
yp <- yz_eye[i_AP[2],] - yz_eye[i_AP[1],]
angR <- atan(yp[2]/yp[1])
yz_eye <- t( matrix(c(cos(angR), sin(angR), -sin(angR), cos(angR)), ncol = 2 ) %*%
               t(sweep(yz_eye, 2, yz_eye[i_AP[1],], '-')) )
# yz_eye <- sweep(yz_eye, 2, yz_eye[i_AP[1],], '-') 
yz_eye <- yz_eye / yz_eye[i_AP[2],1] #normalize
yz_eye <- data.frame(yz_eye)
colnames(yz_eye) <- c('y','z')

# make df for plotting
df <- matrix(ncol = 3, nrow = 0)
for (j in 1:length(LC6) ) {
  df <- rbind(df, as.numeric(c(yz_eye[j,], j)) )
  for (k in 1:N_ind_sp) {
    df <- rbind(df, c(yz_ls[[k]][j,], j) )
  }
}
df <- as.data.frame(df)
colnames(df) <- c('y','z','n')
df$n <- factor(df$n)

df_cable_ax <- df


# pick 4 areas - NEW
i1 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[nt,], '-')^2), decreasing = F)[1:5] #dorsal
i2 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[31,], '-')^2), decreasing = F)[1:(5+0)] #back
i3 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[nb,], '-')^2), decreasing = F)[1:5] #ventral
i4 <- c(1,4,48,46,45)
i5 <- seq(1,65)[-c(i1,i2,i3,i4)]
i_excl <- c(39,14,19,8) #alt tract

col4_a1 <- c('blue', 'brown', 'cyan','pink')
col4 <- adjustcolor(col4_a1, alpha.f = 0.6)

# - Figure 5F
nopen3d()
par3d('windowRect' = c(100,100,2100,2100))
plot3d(LC6[i1], col = col4[1], lwd = 2, soma = T)
# plot3d(LC6[[i1[1]]], col = col4[1], lwd = 2)
plot3d(LC6[i2], col = col4[2], lwd = 2, soma = T)
plot3d(LC6[i3], col = col4[3], lwd = 2, soma = T)
# plot3d(LC6[[i3[1]]], col = col4[3], lwd = 2)
plot3d(LC6[i4], col = col4[4], lwd = 2, soma = T)
# points3d(xyzmatrix(nrefrs$d[ref_ind,]), size = 10)

# draw cross sections
cs4 <- matrix(c(0,1,1,0,-1,1,0,-1,-1,0,1,-1)*1e4, ncol = 3, byrow = T)

for (j in 1:4) {
  cspc <- prcomp(o_xyz_m_rs_sp[(ref_ind_sp[j]-3):(ref_ind_sp[j]+3),])
  # rotm <- cbind(cspc$rotation, cross3D(cspc$rotation[,1], cspc$rotation[,2]))
  rotm <- cspc$rotation
  quad_xyz <- sweep(cs4 %*% t(rotm), MARGIN = 2, STATS = xyzmatrix(nrefrs_sp$d[ref_ind_sp[j],]), '+') #rotate and center
  quads3d(quad_xyz, col = 'black', alpha = 0.9, lit=F)
}

rgl.viewpoint(userMatrix = rotationMatrix(30/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0))
rgl.pop(type = "light")
rgl.light(theta = 30, phi = 30)

# rgl.snapshot(filename = "cross section 3d.png",fmt = "png")

# approx retinotopy LO -> glo -------------------------------------------------------------------------------------

xy_com_m <- matrix(unlist(xy_com), ncol = 4, byrow = T)[,c(4,3)]

getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
synpc_order_pal <-  getPalette(65)


# xyz com in 3D
nopen3d()
for (j in 1:65) {
  text3d(xyz_com[[j]], texts = j)
}


# Figure 5D
windows(record = F, width = 8, height = 8)
bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
bd_theta <- seq(1, 180, by = 1)
xy_bd <- matrix(ncol = 2)
bd_grid <- expand.grid(bd_phi, bd_theta)
plot(bd_grid, ylim = rev(range(bd_grid$Var2)), type = "n", axes = FALSE, ann = F)
for (j in 1:length(xy_poly)) {
  xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
  if (j == 10) {
    polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "#d7191c", density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 3, pch = 3)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
  }
  else if (j == 51) {
    polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "#2c7bb6", density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 3, pch = 3)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
  }
  else {
    # polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "cyan", border = 'black', density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
    text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = j, pos = 1, offset = 0.3)
  }
  # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = synpc_order[j], pos = 1, offset = 0.3)
  
  # points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=synpc_order_pal[synpc_order[j]], cex = 1.5, pch = 20)
  # points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=synpc_order_pal[cable_med_order[j]], cex = 1.5, pch = 20)
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

# pick 4 areas
i1 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[nt,], '-')^2), decreasing = F)[1:5] #dorsal
i2 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[31,], '-')^2), decreasing = F)[1:(5+0)] #back
i3 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[nb,], '-')^2), decreasing = F)[1:5] #ventral
# i4 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[39,], '-')^2), decreasing = F)[1:(9+0)] #front
# i4 <- i4[!(i4 %in% altTract )]
i4 <- c(1,4,48,46,45)
i5 <- seq(1,65)[-c(i1,i2,i3,i4)]
i_excl <- c(39,14,19,8) #alt tract

# i1 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[15,], '-')^2), decreasing = F)[1:2] 
# i2 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[55,], '-')^2), decreasing = F)[1:2] 
# i3 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[38,], '-')^2), decreasing = F)[1:2] 
# i4 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[51,], '-')^2), decreasing = F)[1:2] 

col4_a1 <- c('blue', 'brown', 'cyan','pink')
col4 <- adjustcolor(col4_a1, alpha.f = 0.6)

# - Figure 5F
nopen3d()
plot3d(LC6[i1], col = col4[1], lwd = 2, soma = T)
# plot3d(LC6[[i1[1]]], col = col4[1], lwd = 2)
plot3d(LC6[i2], col = col4[2], lwd = 2, soma = T)
plot3d(LC6[i3], col = col4[3], lwd = 2, soma = T)
# plot3d(LC6[[i3[1]]], col = col4[3], lwd = 2)
plot3d(LC6[i4], col = col4[4], lwd = 2, soma = T)
# points3d(xyzmatrix(nrefrs$d[ref_ind,]), size = 10)

# draw cross sections
cs4 <- matrix(c(0,1,1,0,-1,1,0,-1,-1,0,1,-1)*2e4, ncol = 3, byrow = T)

j <- 1
cspc <- prcomp(xyzmatrix(nrefrs$d[c(ref_ind[j],ref_ind[j]+1),]))
rotm <- cbind(cspc$rotation, cross3D(cspc$rotation[,1], cspc$rotation[,2]))
quad_xyz <- sweep(cs4 %*% t(rotm), MARGIN = 2, STATS = xyzmatrix(nrefrs$d[ref_ind[j],]), '+') #rotate and center
quads3d(quad_xyz, col = 'grey', alpha = 1)

j <- 2
cspc <- prcomp(xyzmatrix(nrefrs$d[c(ref_ind[j],ref_ind[j]+1),]))
rotm <- cbind(cspc$rotation, cross3D(cspc$rotation[,1], cspc$rotation[,2]))
quad_xyz <- sweep(cs4 %*% t(rotm), MARGIN = 2, STATS = xyzmatrix(nrefrs$d[ref_ind[j],]), '+') #rotate and center
quads3d(quad_xyz, col = 'grey', alpha = 0.6)

j <- 3
cspc <- prcomp(xyzmatrix(nrefrs$d[c(ref_ind[j],ref_ind[j]+1),]))
rotm <- cbind(cspc$rotation, cross3D(cspc$rotation[,1], cspc$rotation[,2]))
quad_xyz <- sweep(cs4 %*% t(rotm), MARGIN = 2, STATS = xyzmatrix(nrefrs$d[ref_ind[j],]), '+') #rotate and center
quads3d(quad_xyz, col = 'grey', alpha = 0.6)

rgl.viewpoint(userMatrix = rotationMatrix(30/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0))
rgl.pop(type = "light")
rgl.light(theta = 30, phi = 30)

# rgl.postscript("crossSecton.pdf", "pdf")
# rgl.snapshot(filename = "corssSection.png",fmt = "png")

# - Figure 5C, chosen neurons i1, i2 etc
windows(record = F, width = 8, height = 8)

# pdf(file = "4 group LO.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
bd_theta <- seq(1, 180, by = 1)
xy_bd <- matrix(ncol = 2)
bd_grid <- expand.grid(bd_phi, bd_theta)
plot(bd_grid, ylim = rev(range(bd_grid$Var2)), type = "n", axes = FALSE, ann = F)
i1i <- 1; i2i <- 1; i3i <- 1; i4i <- 1
for (j in 1:length(xy_poly)) {
  xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
  if (j %in% i1) {
    # geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 18 )
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[1], cex = 3, pch = 16)
    # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = i1i, pos = 1, offset = 0.3)
    # i1i <- i1i + 1
  }
  else if (j %in% i2) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[2], cex = 3, pch = 16)
    # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = i2i, pos = 1, offset = 0.3)
    # i2i <- i2i + 1
  }
  else if (j %in% i3) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[3], cex = 3, pch = 16)
    # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = i3i, pos = 1, offset = 0.3)
    # i3i <- i3i + 1
  }
  else if (j %in% i4) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[4], cex = 3, pch = 16)
    # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = i4i, pos = 1, offset = 0.3)
    # i4i <- i4i + 1
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


# # where they are in glo
# # M <- par3d("userMatrix")
# nopen3d(userMatrix = M)
# points3d(syn_xyz_avg[i1,], col=col4[1],size=15)
# points3d(syn_xyz_avg[i2,], col=col4[2],size=15)
# points3d(syn_xyz_avg[i3,], col=col4[3],size=15)
# points3d(syn_xyz_avg[i4,], col=col4[4],size=15)
# for (j in i1) {
#   xyz <- xyzmatrix(LC6_glo_rs[[j]]$d)
#   # points3d(xyz, col = col4[1], size = 5)
#   points3d(matrix(colMeans(xyz), ncol = 3), col = col4[1], size = 15)
#   # plot3d(LC6[[j]], col = col4[1])
# }
# for (j in i2) {
#   xyz <- xyzmatrix(LC6_glo_rs[[j]]$d)
#   # points3d(xyz, col = col4[2], size = 5)
#   points3d(matrix(colMeans(xyz), ncol = 3), col = col4[2], size = 15)
#   # plot3d(LC6[[j]], col = col4[2])
# }
# for (j in i3) {
#   xyz <- xyzmatrix(LC6_glo_rs[[j]]$d)
#   # points3d(xyz, col = col4[3], size = 5)
#   points3d(matrix(colMeans(xyz), ncol = 3), col = col4[3], size = 15)
#   # plot3d(LC6[[j]], col = col4[3])
# }
# for (j in i4) {
#   xyz <- xyzmatrix(LC6_glo_rs[[j]]$d)
#   # points3d(xyz, col = col4[4], size = 5)
#   points3d(matrix(colMeans(xyz), ncol = 3), col = col4[4], size = 15)
#   # plot3d(LC6[[j]], col = col4[4])
# }


nopen3d()
for (j in i3) {
  conn <- connectors(LC6[[j]])
  conn_pre <- conn[conn$prepost== 0, c('x','y','z')]
  inglo <- rowSums(sweep(conn_pre,2,c(a,b,c), '*')) + d < 0
  points3d(conn_pre[inglo,], col = col4[2], size = 5)
}

# mean and x-median syn position
syn_xyz_avg <- matrix(ncol = 3, nrow = length(LC6))
syn_xyz_med_x <- matrix(ncol = 1, nrow = length(LC6))
for (j in 1:length(LC6)) {
  conn <- connectors(LC6[[j]])
  conn_pre <- conn[conn$prepost==1, c('x','y','z')]
  inglo <- rowSums(sweep(conn_pre,2,c(a,b,c), '*')) + d < 0
  syn_xyz_avg[j,] <- colMeans(conn_pre[inglo,])
  syn_xyz_med_x[j,] <- median(conn_pre[inglo,1])
}

# nopen3d()
# points3d(syn_xyz_avg)
synpc <- prcomp(syn_xyz_avg)
syn_xyz_avg_rot <- t( t(synpc$rotation) %*% t(syn_xyz_avg) )
nopen3d()
points3d(syn_xyz_avg_rot)

synpc_order <- order(syn_xyz_avg_rot[,1])

# dev.new()
# plot(syn_xyz_avg[,1])
# plot(syn_xyz_avg[,2])
# plot(syn_xyz_avg[,3])


# # - cable plot, 4 groups
# 
# plt <- ggplot()
# for (j in i1) {
#   df_plt <- df_cable[((N_ind+1)*(j-1)+1):((N_ind+1)*j),]
#   df_plt <- cbind(df_plt, seq(1,N_ind+1)*0.1+1)
#   colnames(df_plt)[4] <- 'N_ind'
#   plt <- plt + 
#     # geom_point(data = df_plt, aes(y, z) ) +
#     geom_path(data = df_plt, aes(x=y,y=z), colour = col4[1], lwd =2) +
#     geom_point(data = yz_eye[j,], aes(x=y,y=z), size = 6, shape = 18)
# }
# for (j in i2) {
#   df_plt <- df_cable[((N_ind+1)*(j-1)+1):((N_ind+1)*j),]
#   df_plt <- cbind(df_plt, seq(1,N_ind+1)*0.1+1)
#   colnames(df_plt)[4] <- 'N_ind'
#   plt <- plt + 
#     # geom_point(data = df_plt, aes(y, z) ) +
#     geom_path(data = df_plt, aes(x=y,y=z), colour = col4[2], lwd =2) +
#     geom_point(data = yz_eye[j,], aes(x=y,y=z), size = 6, shape = 18)
# }
# for (j in i3) {
#   df_plt <- df_cable[((N_ind+1)*(j-1)+1):((N_ind+1)*j),]
#   df_plt <- cbind(df_plt, seq(1,N_ind+1)*0.1+1)
#   colnames(df_plt)[4] <- 'N_ind'
#   plt <- plt + 
#     # geom_point(data = df_plt, aes(y, z) ) +
#     geom_path(data = df_plt, aes(x=y,y=z), colour = col4[3], lwd =2) +
#     geom_point(data = yz_eye[j,], aes(x=y,y=z), size = 6, shape = 18)
# }
# for (j in i4) {
#   df_plt <- df_cable[((N_ind+1)*(j-1)+1):((N_ind+1)*j),]
#   df_plt <- cbind(df_plt, seq(1,N_ind+1)*0.1+1)
#   colnames(df_plt)[4] <- 'N_ind'
#   plt <- plt + 
#     # geom_point(data = df_plt, aes(y, z) ) +
#     geom_path(data = df_plt, aes(x=y,y=z), colour = col4[4], lwd =2) +
#     geom_point(data = yz_eye[j,], aes(x=y,y=z), size = 6, shape = 18)
# }
# 
# dev.new()
# plt 

# cable plot axon -----------------------------------------------------------------------------------------------

# Figure 5F
windows(record = F, width = 8, height = 8)

# pdf(file = "4 group LO.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
bd_theta <- seq(1, 180, by = 1)
xy_bd <- matrix(ncol = 2)
bd_grid <- expand.grid(bd_phi, bd_theta)
plot(bd_grid, ylim = rev(range(bd_grid$Var2)), type = "n", axes = FALSE, ann = F)
i1i <- 1; i2i <- 1; i3i <- 1; i4i <- 1
for (j in 1:length(xy_poly)) {
  xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
  if (j %in% i1) {
    # geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 18 )
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[1], cex = 3, pch = 16)
    # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = i1i, pos = 1, offset = 0.3)
    # i1i <- i1i + 1
  }
  else if (j %in% i2) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[2], cex = 3, pch = 16)
    # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = i2i, pos = 1, offset = 0.3)
    # i2i <- i2i + 1
  }
  else if (j %in% i3) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[3], cex = 3, pch = 16)
    # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = i3i, pos = 1, offset = 0.3)
    # i3i <- i3i + 1
  }
  else if (j %in% i4) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[4], cex = 3, pch = 16)
    # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = i4i, pos = 1, offset = 0.3)
    # i4i <- i4i + 1
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



# m <- 1
# df_plt <- df_cable_ax[((N_ind_sp+1)*(c(i1,i2,i3,i4)-1)+m),]
# df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
#                                  rep(2,length(i2)),
#                                  rep(3,length(i3)),
#                                  rep(4,length(i4)))) )
# colnames(df_plt)[4] <- 'colcol'
# gpl <- ggplot() +
#   geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 16 ) +
#   scale_colour_manual(values = col4,
#                       breaks = c("1", "2", "3","4"),
#                       labels = c("top", "back", "bottom","front") ) +
#   guides(colour = guide_legend("groups") )+
#   # xlim(xrange) +
#   # ylim(yrange) +
#   scale_y_reverse() +
#   coord_fixed(ratio = 1) +
#   # theme_minimal() +
#   theme_void() +
#   geom_segment(aes(x = 0, y = -0.5, xend = 1/ref_d[m]*1e5/10, yend = -0.5), size=2, lineend = "round") +
#   annotate("text", x = 0.05, y = -0.6, label = '10 um') +
#   # annotate("text", x = -.5, y = -0.6, label = paste('retino-index = ', mean(rt_ind[,m]))) +
#   labs(title = '1st postion (in LO)')
# for (j in 1:4) {
#   for (k in 1:5) {
#     gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k]+.03, y=df_plt$z[5*(j-1)+k]+.03, label = k)
#     # gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k], y=df_plt$z[5*(j-1)+k], label = k, colour='white') 
#   }
# }
# windows(record = F, width = 8, height = 8)
# # pdf(file = "crossSecton 1.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
# gpl
# dev.off()



m <- 4 #5, 6
df_plt <- df_cable_ax[((N_ind_sp+1)*(c(i1,i2,i3,i4)-1)+m),]
df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
                                 rep(2,length(i2)),
                                 rep(3,length(i3)),
                                 rep(4,length(i4)))) )
colnames(df_plt)[4] <- 'colcol'
gpl <- ggplot() +
  geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 20 ) +
  scale_colour_manual(values = col4,
                      breaks = c("1", "2", "3","4"),
                      labels = c("top", "back", "bottom","front") ) +
  guides(colour = guide_legend("groups") )+
  scale_y_reverse() +
  # xlim(xrange) +
  # ylim(yrange) +
  coord_fixed(ratio = 1) +
  # theme_minimal() +
  theme_void() +
  geom_segment(aes(x = 0, y = -0.5, xend = 1/ref_d[m]*1e3, yend = -0.5), size=2, lineend = "round") +
  annotate("text", x = 0.05, y = -0.8, label = '1 um') +
  # annotate("text", x = -.5, y = -2.8, label = paste('retino-index = ', mean(rt_ind[,m]))) +
  # labs(title = paste('retino-index = ', round(mean(rt_ind[,m]),digits = 2)))
  labs(title = m-1)

# windows(record = F, width = 8, height = 8)
pdf(file = paste(m, ".pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl

dev.off()



# all 65
m <- 5
df_plt <- df_cable_ax[((N_ind_sp+1)*(c(i1,i2,i3,i4,i5)-1)+m),]
df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
                                 rep(2,length(i2)),
                                 rep(3,length(i3)),
                                 rep(4,length(i4)),
                                 rep(5,length(i5)))), c(rep(2, 20), 
                                                        rep(1, 45)) )
colnames(df_plt)[4] <- 'colcol'
colnames(df_plt)[5] <- 'dotsize'
# legend_size <- c(3, 4, 5, 6, 7)
gpl <- ggplot(data = df_plt, aes(x=y, y=z)) +
  geom_point(aes(colour = colcol, size = dotsize), shape = 20 ) +
  scale_colour_manual(values = c(col4, 'grey'),
                      breaks = c("1", "2", "3","4","5"),
                      labels = c("top", "back", "bottom","front","other") ) +
  guides(colour = guide_legend("groups") )+
  scale_y_reverse() +
  scale_size(range = c(6,10)) +
  # geom_text(hjust=0, vjust=0, label=c(i1,i2,i3,i4,i5) ) +
  coord_fixed(ratio = 1) +
  # theme_minimal() +
  theme_void() +
  geom_segment(aes(x = 0, y = -0.5, xend = 1/ref_d[m]*1e3, yend = -0.5), size=2, lineend = "round") +
  annotate("text", x = 0.05, y = -0.8, label = '1 um') +
  #guides(color= guide_legend(), size = guide_legend(override.aes = list(size = legend_size)) ) +
  # annotate("text", x = -.5, y = -2.8, label = paste('retino-index = ', mean(rt_ind[,m]))) +
  # labs(title = paste('retino-index = ', round(mean(rt_ind[,m]),digits = 2)))
  labs(title = m-1)

# windows(record = F, width = 8, height = 8)
pdf(file = paste("65 - ", m, ".pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl

dev.off()

# cable plot in glo -----------------------------------------------------------------------------------------------

# where they are in glo
# M <- par3d("userMatrix")
# nopen3d(userMatrix = M)
# nopen3d()
# plot3d(LC6_glo_rs[[i1[1]]], col = col4[1], WithNodes = F, lwd =2)
# plot3d(LC6_glo_rs[[i2[1]]], col = col4[2], WithNodes = F, lwd =2)
# plot3d(LC6_glo_rs[[i3[1]]], col = col4[3], WithNodes = F, lwd =2)
# plot3d(LC6_glo_rs[[i4[1]]], col = col4[4], WithNodes = F, lwd =2)


m <- 1
df_plt <- df_cable_sp[((N_ind_sp+1)*(c(i1,i2,i3,i4)-1)+m),]
df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
                                 rep(2,length(i2)),
                                 rep(3,length(i3)),
                                 rep(4,length(i4)))) )
colnames(df_plt)[4] <- 'colcol'
gpl <- ggplot() +
  geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 16 ) +
  scale_colour_manual(values = col4,
                      breaks = c("1", "2", "3","4"),
                      labels = c("top", "back", "bottom","front") ) +
  guides(colour = guide_legend("groups") )+
  # xlim(xrange) +
  # ylim(yrange) +
  coord_fixed(ratio = 1) +
  # theme_minimal() +
  theme_void() +
  geom_segment(aes(x = 0, y = -0.5, xend = 1/ref_d[m]*1e5/10, yend = -0.5), size=2, lineend = "round") +
  annotate("text", x = 0.05, y = -0.6, label = '10 um') +
  # annotate("text", x = -.5, y = -0.6, label = paste('retino-index = ', mean(rt_ind[,m]))) +
  labs(title = '1st postion (in LO)')
for (j in 1:4) {
  for (k in 1:5) {
    gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k]+.03, y=df_plt$z[5*(j-1)+k]+.03, label = k)
    # gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k], y=df_plt$z[5*(j-1)+k], label = k, colour='white') 
  }
}
windows(record = F, width = 8, height = 8)
# pdf(file = "crossSecton 1.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl
dev.off()


m <- 2 #5, 6
df_plt <- df_cable_sp[((N_ind_sp+1)*(c(i1,i2,i3,i4)-1)+m),]
df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
                                 rep(2,length(i2)),
                                 rep(3,length(i3)),
                                 rep(4,length(i4)))) )
colnames(df_plt)[4] <- 'colcol'
gpl <- ggplot() +
  geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 20 ) +
  scale_colour_manual(values = col4,
                      breaks = c("1", "2", "3","4"),
                      labels = c("top", "back", "bottom","front") ) +
  guides(colour = guide_legend("groups") )+
  # xlim(xrange) +
  # ylim(yrange) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  # theme_void() +
  # geom_segment(aes(x = 0, y = -2.5, xend = 1/ref_d[m]*1e3, yend = -2.5), size=2, lineend = "round") +
  # annotate("text", x = 0.05, y = -2.8, label = '1 um') +
  # annotate("text", x = -.5, y = -2.8, label = paste('retino-index = ', mean(rt_ind[,m]))) +
  # labs(title = paste('retino-index = ', round(mean(rt_ind[,m]),digits = 2)))
  labs(title = m-1)

# windows(record = F, width = 8, height = 8)
pdf(file = paste(m, ".pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl

dev.off()



# all 65
m <- 3 #5, 6
df_plt <- df_cable_sp[((N_ind_sp+1)*(c(i1,i2,i3,i4,i5)-1)+m),]
df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
                                 rep(2,length(i2)),
                                 rep(3,length(i3)),
                                 rep(4,length(i4)),
                                 rep(5,length(i5)))) )
colnames(df_plt)[4] <- 'colcol'
gpl <- ggplot(data = df_plt, aes(x=y, y=z)) +
  geom_point(aes(colour = colcol), size = 5, shape = 20 ) +
  scale_colour_manual(values = c(col4, 'grey'),
                      breaks = c("1", "2", "3","4","5"),
                      labels = c("top", "back", "bottom","front","other") ) +
  guides(colour = guide_legend("groups") )+
  geom_text(hjust=0, vjust=0, label=c(i1,i2,i3,i4,i5) ) +
  # xlim(xrange) +
  # ylim(yrange) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  # theme_void() +
  # geom_segment(aes(x = 0, y = -2.5, xend = 1/ref_d[m]*1e3, yend = -2.5), size=2, lineend = "round") +
  # annotate("text", x = 0.05, y = -2.8, label = '1 um') +
  # annotate("text", x = -.5, y = -2.8, label = paste('retino-index = ', mean(rt_ind[,m]))) +
  # labs(title = paste('retino-index = ', round(mean(rt_ind[,m]),digits = 2)))
  labs(title = m-1)

# windows(record = F, width = 8, height = 8)
pdf(file = paste("65 - ", m, ".pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl

dev.off()

#  - cable plot, 4 groups, 4 separate locations ---------------------------------------------------------------------

# xrange <- c(-1, 2)
# yrange <- c(-1, 3)

# Figure 5F

m <- 1
df_plt <- df_cable[((N_ind+1)*(c(i1,i2,i3,i4)-1)+m),]
df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
                                 rep(2,length(i2)),
                                 rep(3,length(i3)),
                                 rep(4,length(i4)))) )
colnames(df_plt)[4] <- 'colcol'
gpl <- ggplot() +
  geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 16 ) +
  scale_colour_manual(values = col4,
                      breaks = c("1", "2", "3","4"),
                      labels = c("top", "back", "bottom","front") ) +
  guides(colour = guide_legend("groups") )+
  # xlim(xrange) +
  # ylim(yrange) +
  coord_fixed(ratio = 1) +
  # theme_minimal() +
  theme_void() +
  geom_segment(aes(x = 0, y = -0.5, xend = 1/ref_d[m]*1e5/10, yend = -0.5), size=2, lineend = "round") +
  annotate("text", x = 0.05, y = -0.6, label = '10 um') +
  annotate("text", x = -.5, y = -0.6, label = paste('retino-index = ', mean(rt_ind[,m]))) +
  labs(title = '1st postion (in LO)')
for (j in 1:4) {
  for (k in 1:5) {
    gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k]+.03, y=df_plt$z[5*(j-1)+k]+.03, label = k)
    # gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k], y=df_plt$z[5*(j-1)+k], label = k, colour='white') 
  }
}
windows(record = F, width = 8, height = 8)

# pdf(file = "crossSecton 1.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl
dev.off()

m <- 2
df_plt <- df_cable[((N_ind+1)*(c(i1,i2,i3,i4)-1)+m),]
df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
                                 rep(2,length(i2)),
                                 rep(3,length(i3)),
                                 rep(4,length(i4)))) )
colnames(df_plt)[4] <- 'colcol'
gpl <- ggplot() +
  geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 20 ) +
  scale_colour_manual(values = col4,
                      breaks = c("1", "2", "3","4"),
                      labels = c("top", "back", "bottom","front") ) +
  guides(colour = guide_legend("groups") )+
  # xlim(xrange) +
  # ylim(yrange) +
  coord_fixed(ratio = 1) +
  # theme_minimal() +
  theme_void() +
  geom_segment(aes(x = 0, y = -2.5, xend = 1/ref_d[m]*1e3, yend = -2.5), size=2, lineend = "round") +
  annotate("text", x = 0.05, y = -2.8, label = '1 um') +
  # annotate("text", x = -.5, y = -2.8, label = paste('retino-index = ', mean(rt_ind[,m]))) +
  labs(title = paste('2nd postion (close to LO), retino-index = ', round(mean(rt_ind[,m]),digits = 2)))
# for (j in 1:4) {
#   for (k in 1:5) {
#     gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k]+.2, y=df_plt$z[5*(j-1)+k]+.2, label = k)
#     # gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k], y=df_plt$z[5*(j-1)+k], label = k, colour='white') 
#   }
# }
# windows(record = F, width = 8, height = 8)

pdf(file = "crossSecton 2.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl
dev.off()

m <- 3
df_plt <- df_cable[((N_ind+1)*(c(i1,i2,i3,i4)-1)+m),]
df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
                                 rep(2,length(i2)),
                                 rep(3,length(i3)),
                                 rep(4,length(i4)))) )
colnames(df_plt)[4] <- 'colcol'
gpl <- ggplot() +
  geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 20 ) +
  scale_colour_manual(values = col4,
                      breaks = c("1", "2", "3","4"),
                      labels = c("top", "back", "bottom","front") ) +
  guides(colour = guide_legend("groups") )+
  # xlim(xrange) +
  # ylim(yrange) +
  coord_fixed(ratio = 1) +
  # theme_minimal() +
  theme_void() +
  geom_segment(aes(x = 0, y = -2.5, xend = 1/ref_d[m]*1e3, yend = -2.5), size=2, lineend = "round") +
  annotate("text", x = 0.5, y = -2.8, label = '1 um') +
  labs(title = paste('3rd, retino-index = ', round(mean(rt_ind[,m]),digits = 2)))
# for (j in 1:4) {
#   for (k in 1:5) {
#     gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k]+.4, y=df_plt$z[5*(j-1)+k]+.4, label = k)
#     # gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k], y=df_plt$z[5*(j-1)+k], label = k, colour='white') 
#   }
# }
# windows(record = F, width = 8, height = 8)

pdf(file = "crossSecton 3.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl
dev.off()


m <- 4
df_plt <- df_cable[((N_ind+1)*(c(i1,i2,i3,i4)-1)+m),]
df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
                                 rep(2,length(i2)),
                                 rep(3,length(i3)),
                                 rep(4,length(i4)))) )
colnames(df_plt)[4] <- 'colcol'
gpl <- ggplot() +
  geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 20 ) +
  # geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), shape = '.' ) +
  scale_colour_manual(values = col4,
                      breaks = c("1", "2", "3","4"),
                      labels = c("top", "back", "bottom","front") ) +
  guides(colour = guide_legend("groups") )+
  # xlim(xrange) +
  # ylim(yrange) +
  coord_fixed(ratio = 1) +
  # theme_minimal() +
  theme_void() +
  geom_segment(aes(x = 0, y = -2.5, xend = 1/ref_d[m]*1e3, yend = -2.5), size=2, lineend = "round") +
  annotate("text", x = 0.5, y = -2.8, label = '1 um') +
  labs(title = paste('4th, retino-index = ', round(mean(rt_ind[,m]),digits = 2)))
# for (j in 1:4) {
#   for (k in 1:5) {
#     gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k]+.1, y=df_plt$z[5*(j-1)+k]+.1, label = k)
#     # gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k], y=df_plt$z[5*(j-1)+k], label = k, colour='white') 
#   }
# }
# windows(record = F, width = 8, height = 8)
# gpl

pdf(file = "crossSecton 4.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl
dev.off()

# j <- 3
# df_plt <- df_cable[((N_ind+1)*(c(i1,i2,i3,i4,i5)-1)+j),]
# df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
#                                  rep(2,length(i2)),
#                                  rep(3,length(i3)),
#                                  rep(4,length(i4)),
#                                  rep(5,length(i5)))) )
# colnames(df_plt)[4] <- 'colcol'
# gpl <- ggplot() +
#   geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 18 ) +
#   scale_colour_manual(values = c(col4,'grey'),
#                       breaks = c("1", "2", "3","4","5"),
#                       labels = c("top", "back", "bottom","front","rest") ) +
#   guides(colour = guide_legend("groups") )+
#   # xlim(xrange) +
#   # ylim(yrange) +
#   coord_fixed(ratio = 1) +
#   theme_minimal() +
#   labs(title = '3rd')
# for (j in 1:4) {
#   for (k in 1:5) {
#     gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k]+.05, y=df_plt$z[5*(j-1)+k]+.05, label = k)
#     # gpl <- gpl + annotate("text", x=df_plt$y[5*(j-1)+k], y=df_plt$z[5*(j-1)+k], label = k, colour='white') 
#   }
# }
# dev.new()
# gpl


# j <- 4
# df_plt <- df_cable[((N_ind+1)*(c(i1,i2,i3,i4)-1)+j),]
# df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
#                                  rep(2,length(i2)),
#                                  rep(3,length(i3)),
#                                  rep(4,length(i4)))) )
# colnames(df_plt)[4] <- 'colcol'
# dev.new()
# ggplot() +
#   geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 6, shape = 18 ) +
#   scale_colour_manual(values = col4,
#                       breaks = c("1", "2", "3","4"),
#                       labels = c("top", "back", "bottom","front") ) +
#   guides(colour = guide_legend("groups") )+
#   # xlim(xrange) +
#   # ylim(yrange) +
#   labs(title = '4th position')



# all 65 LC at 1st location
j <- 1
# df_plt <- df_cable[((N_ind+1)*(c(i1,i2,i3,i4)-1)+j),]
df_plt <- df_cable[((N_ind+1)*(seq(1,65)-1)+j),]
# df_plt <- cbind(df_plt, factor(rep(c(1,2,3,4),each = N_ind)) )
# df_plt <- cbind(df_plt, factor(seq(1,65))) 
# colnames(df_plt)[4] <- 'colcol'

gg <- ggplot() +
  geom_point(data = df_plt, aes(x=y, y=z), size = 0, shape = 18 ) +
  labs(title = '1st postion (close to LO)')
for (j in 1:65) {
  gg <- gg + annotate("text",x=df_plt$y[j], y=df_plt$z[j], label = j) 
}
dev.new()  
gg

  

# permutation index -----------------------------------------------------------------------------------------------

# fun arclength
arcLength <- function(p1, p2) {
  p1_xyz <- c(sin(p1[1])*cos(p1[2]), sin(p1[1])*sin(p1[2]), cos(p1[1]))
  p2_xyz <- c(sin(p2[1])*cos(p2[2]), sin(p2[1])*sin(p2[2]), cos(p2[1]))
  c2 <- sum((p1_xyz - p2_xyz)^2)
  arc_ang <- acos((2-c2)/2)
  return(arc_ang*1)
}

# library(permutations) #conflict with sth

# bubble sort function, x is member of Sn, symmetric group of n letters from an ordered set
bubble_sort <- function(x){ 
  n <- length(x)
  N_swap <- 0
  repeat{
    swapped = FALSE
    for (j in 1:(n-1)) {
      if (x[j] > x[j+1]) {
        x[c(j,j+1)] <- x[c(j+1,j)]
        swapped = TRUE
        N_swap <- N_swap + 1
      }
    }
    n <- n - 1
    if (swapped == FALSE | n == 1) {
      break()
    }
  }
  return(list(x, N_swap))
}

# alt, find the nth letter and insert at nth position, works for set(seq(1,n))
bubble_sort_swapCount <- function(x){
  n <- length(x)
  N_swap <- 0
  for (j in 1:(n-1)) {
    ind <- which(x %in% j)
    x <- x[-ind]
    N_swap <- N_swap + ind - 1
  }
  return(N_swap)
}


# # -- test time
# j <- 8
# mat <- gtools::permutations(j,j)
# count <- c()
# start_time <- Sys.time()
# for (k in 1:dim(mat)[1]) {
#   count <- c(count, bubble_sort_swapCount(mat[k,]) )
#   # count <- c(count, bubble_sort(mat[k,])[2] )
# }
# end_time <- Sys.time()
# end_time - start_time






# 
# o34 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[34,], '-')^2), decreasing = F) #dorsal
# 
# # -- order avg syn position
# o34_syn <- order(rowSums(sweep(syn_xyz_avg, 2, syn_xyz_avg[34,], '-')^2)[o34]) #order first by o34, then by distance
# bubble_sort_swapCount(o34_syn)
# 
# # -- order x median syn position
# o34_syn_med <- order(((syn_xyz_med_x - syn_xyz_med_x[34])^2)[o34]) #order first by o34, then by distance
# bubble_sort_swapCount(o34_syn_med)




#  skeleton med, mean, tip along glomerulus ----------------------------------------------------------------------------------

nopen3d()
points3d(conn_LC6[,c('x','y','z')])
# identify3d(conn_LC6[,c('x','y','z')])

# make a glo mesh
conn_LC6_gloHull <- conn_LC6[-c(4604, 1498, 1343, 1342,  122, 2345, 2782, 5828, 5182,
                                4018,  575, 5646, 5679, 5631, 5640, 5654, 5652),
                             c('x','y','z')]
glo.as <- ashape3d(as.matrix(conn_LC6_gloHull), alpha = 10000)
glo.msh <- as.mesh3d(glo.as)

# filter for neurites within the mesh
LC6_glo <- nlapply(LC6, subset, function(x) pointsinside(x, surf=glo.msh, rval='distance') > -5000)


nopen3d()
plot3d(LC6_glo)
shade3d(glo.msh, alpha = 1)
# points3d(conn_LC6_gloHull)

# resample
LC6_glo_rs <- resample(LC6_glo, stepsize = 400)


# # - neurites com within compartments
# xdiv <- quantile(range(conn_LC6_gloHull$x), seq(0,1,length.out = 4))
# com_cable <- list()
# for (j in 1:3) {
#   tmp <- matrix(ncol = 3, nrow = 65)
#   for (k in 1:65) {
#     xyz <- xyzmatrix(LC6_glo_rs[[k]]$d)
#     xyz_div <- xyz[xyz[,1] > xdiv[j] & xyz[,1] < xdiv[j+1],]
#     if (length(xyz_div) > 3) {
#       tmp[k,] <- colMeans(xyz_div)
#     }
#     # else {
#     #   tmp[k,] <- xyz_div
#     # }
#   }
#   com_cable[[j]] <- tmp
# }


# - neurite median in x
cable_glo_avg <- matrix(ncol = 3, nrow = 65)
cable_glo_med_x <- c()
cable_glo_med <- matrix(ncol = 3, nrow = 65)
cable_glo_mean_x <- c()
cable_glo_mid_x <- c()
cable_glo_min10_x <- c()
cable_glo_min10 <- matrix(ncol = 3, nrow = 65)
cable_glo_min1_x <- c()
cable_glo_min1 <- matrix(ncol = 3, nrow = 65)
for (j in 1:65) {
  xyz <- xyzmatrix(LC6_glo_rs[[j]]$d)
  cable_glo_avg[j,] <- colMeans(xyz)
  cable_glo_med_x <- c(cable_glo_med_x, median(xyz[,1]) )
  cable_glo_med[j,] <- xyz[which.min(abs(xyz[,1]-median(xyz[,1]))), ]
  cable_glo_mean_x <- c(cable_glo_mean_x, mean(xyz[,1]) )
  cable_glo_mid_x <- c(cable_glo_mid_x, mean(range(xyz[,1])) )
  xyz10 <- xyz[xyz[,1] < quantile(xyz[,1], 0.1),]
  cable_glo_min10_x <- c(cable_glo_min10_x, mean(xyz10[,1]) )
  cable_glo_min10[j,] <- colMeans(xyz10)
  xyz1 <- xyz[xyz[,1] < quantile(xyz[,1], 0.01),]
  cable_glo_min1_x <- c(cable_glo_min1_x, mean(xyz1[,1]) )
  cable_glo_min1[j,] <- colMeans(xyz1)
}


# - PLOT, order by median x in LO
cable_glo_x_o <- order(cable_glo_med_x)
# cable_glo_x_o <- order(cable_glo_mean_x)
# cable_glo_x_o <- order(cable_glo_mid_x)
# cable_glo_x_o <- order(cable_glo_min10_x)
# Figure 5D
windows(record = F, width = 8, height = 8)
bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
bd_theta <- seq(1, 180, by = 1)
xy_bd <- matrix(ncol = 2)
bd_grid <- expand.grid(bd_phi, bd_theta)
plot(bd_grid, ylim = rev(range(bd_grid$Var2)), type = "n", axes = FALSE, ann = F)
for (j in 1:length(xy_poly)) {
  xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
  # points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=synpc_order_pal[cable_glo_x_o[j]], cex = 1.5, pch = 20)
  points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=synpc_order_pal[match(j,cable_glo_x_o)], cex = 1.5, pch = 20)
  # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = cable_glo_x_o[j], pos = 1, offset = 0.3)
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
title('x median')


# Figure 6E 
# use 3 colors
# synpc_order_pal_3 <- brewer.pal(3, "Blues")
synpc_order_pal_3 <- c(brewer.pal(5, "RdBu")[c(1,2)], "grey40")

# cable_glo_x_o <- order(cable_glo_avg)
# cable_glo_x_o <- order(cable_glo_min10_x)
cable_glo_x_o <- order(cable_glo_med_x)
# cable_glo_x_o <- order(cable_glo_min10_x)

cable_glo_x_o <- data.frame(order = cable_glo_x_o)
cable_glo_x_o$gp <- 0
cable_glo_x_o$gp[cable_glo_x_o$order[1:22]] <- 1
cable_glo_x_o$gp[cable_glo_x_o$order[23:44]] <- 2
cable_glo_x_o$gp[cable_glo_x_o$order[43:65]] <- 3

windows(record = F, width = 8, height = 8)

# pdf(file = "x min 10.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
# pdf(file = "x median.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
bd_theta <- seq(1, 180, by = 1)
xy_bd <- matrix(ncol = 2)
bd_grid <- expand.grid(bd_phi, bd_theta)
plot(bd_grid, ylim = rev(range(bd_grid$Var2)), type = "n", axes = FALSE, ann = F)
for (j in 1:length(xy_poly)) {
  xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
  points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=synpc_order_pal_3[cable_glo_x_o[j,'gp']], cex = 3, pch = 16)
  # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = cable_glo_x_o[j], pos = 1, offset = 0.3)
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
# title('x min 10')
title('x median')
# title('x min 1')

dev.off()

# Figure 6F

# 3 color glo
nopen3d()
par3d('windowRect' = c(100,100,3100,1100))
rgl.viewpoint(userMatrix = rotationMatrix(30/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0)
              , zoom = 0.3)

# points3d(cable_glo_avg, col = synpc_order_pal_3[cable_glo_x_o[,'gp']], size = 30)
# points3d(cable_glo_med, col = synpc_order_pal_3[cable_glo_x_o[,'gp']], size = 30)
points3d(cable_glo_min10, col = synpc_order_pal_3[cable_glo_x_o[,'gp']], size = 30)
# points3d(cable_glo_min1, col = synpc_order_pal_3[cable_glo_x_o[,'gp']], size = 10)
# shade3d(glo.msh, alpha = 0.1)
shade3d(glo.msh, col='cyan',alpha = 0.05)

rgl.pop(type = "light")
rgl.light(theta = 30, phi = 30)

# title3d('x median')
title3d('cable min 10')

# rgl.snapshot(filename = "3gp min 10.png",fmt = "png")

# dev.off()

# df <- cable_glo_x_o
# df$x <- cable_glo_min10_x
# dev.new()
# ggplot(df) + geom_point(aes(x, y=0, colour = synpc_order_pal_3[gp]), size =10)
# ggplot(df) + geom_point(aes(x, y=0, colour = gp))

# x min 10, rot scatter plot ------------------------------------------------------------------------------------------

# add a rotation

k_fit <- c()
pval <- c()
rsq <- c()
xy_com_m_tp <- xy_com_m[,c(2,1)] / 180 * pi
cable_glo_med_x_norm <- (cable_glo_med_x - glu_div[1])/ diff(range(glu_div)) 
# rot axis
vr <- c(cos(sum(range(xy_com_m[,1]))/2/180*pi), sin(sum(range(xy_com_m[,1]))/2/180*pi),0)
for (j in 1:36) {
  # ang <- j*10/180*pi
  # M_rot <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol = 2, byrow = T)
  # xy_com_m_rot <- xy_com_m %*% M_rot
  
  ang <- j*10
  R_mat <- quaternion3D(vr, ang)
  
  xy_com_m_rot <- sph2cartZ(cbind(1, xy_com_m_tp)) %*% t(R_mat) #ccw on S2
  xy_com_m_rot <- cart2sphZ(xy_com_m_rot)[, c(3,2)]
  xy_com_m_rot[,1] <- if_else(xy_com_m_rot[,1] > pi, xy_com_m_rot[,1] - 2*pi, xy_com_m_rot[,1])
  xy_com_m_rot <- xy_com_m_rot/pi*180
  
  # normalize
  # com_rot_norm <- xy_com_m_rot[,1] / 360
  x_com_m_rot_norm <- (xy_com_m_rot[,1] - min(xy_com_m_rot[,1])) / diff(range(xy_com_m_rot[,1]))
  
  
  # dev.new()
  # plot(xy_com_m_rot)
  # points(xy_com_m_rot[c(1,2),], cex = 2, col ='red', pch = 16)
  # points(xy_com_m_rot[c(3,4),], cex = 2, col ='blue', pch = 16)
  
  # df <- cbind.data.frame(xy_com_m_rot[,1], cable_glo_med_x)
  
  # normalize to [0, 1]
  # cable_glo_min10_x_norm <- cable_glo_min10_x - min(cable_glo_min10_x)
  # cable_glo_min10_x_norm <- cable_glo_min10_x_norm / max(cable_glo_min10_x_norm)
  
  
  df <- cbind.data.frame(x_com_m_rot_norm, cable_glo_med_x_norm)
  colnames(df) <- c('x','y')
  linReg <- lm(y ~ x, data = df)
  # df$fitted <- linReg$fitted.values
  k_fit <- c(k_fit, linReg$coefficients['x'])
  pval <- c(pval, summary(linReg)$coefficients['x', 'Pr(>|t|)'] ) #p-value
  rsq <- c(rsq, summary(linReg)$r.squared ) # R^2
  
}
which.max(k_fit)
which.min(k_fit)

dev.new()
# pdf(file = "reg slope vs ang.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
plot(1:36, k_fit, type='p', pch = 16
     ,xaxt='none', yaxt='n'
     ,xlab = 'rotation angle', ylab = 'reg slope'
     ,xlim = c(0,36), ylim = c(-.5,0.5)
     )
axis(1, at = seq(0, 36, by = 9), labels = seq(0, 360, by = 90), las=1) 
axis(2, at = seq(-0.5, 0.5, by = 0.5), labels = seq(-0.5, 0.5, by = 0.5), las=2) 

dev.off()

# best choice
j <- 3
# j <- 21 #min
# j <- 31
# j <- 4

# ang <- j*10
# R_mat <- quaternion3D(vr, ang)
# 
# xy_com_m_rot <- sph2cartZ(cbind(1, xy_com_m_tp)) %*% t(R_mat) #ccw on S2
# xy_com_m_rot <- cart2sphZ(xy_com_m_rot)[, c(3,2)]
# xy_com_m_rot[,1] <- if_else(xy_com_m_rot[,1] > pi, xy_com_m_rot[,1] - 2*pi, xy_com_m_rot[,1])
# xy_com_m_rot <- xy_com_m_rot/pi*180
# 
# # ang <- j*10/180*pi
# # M_rot <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol = 2, byrow = T)
# # 
# # xy_com_m_rot <- xy_com_m %*% M_rot
# 
# # normalize
# com_rot_norm <- xy_com_m_rot[,1] / 360
# cable_glo_med_x_norm <- (cable_glo_med_x - glu_div[1])/ diff(range(glu_div)) 

# -- RI for all rot ang
N_swap_rot <- matrix(ncol = 18, nrow = 65)
for (k in 1:18) {
  # ang <- (k-1)*10/180*pi
  # M_rot <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol = 2, byrow = T)
  # xy_com_m_rot <- xy_com_m %*% M_rot

  ang <- (k-1)*10
  R_mat <- quaternion3D(vr, ang)
  
  xy_com_m_rot <- sph2cartZ(cbind(1, xy_com_m_tp)) %*% t(R_mat) #ccw on S2
  xy_com_m_rot <- cart2sphZ(xy_com_m_rot)[, c(3,2)]
  xy_com_m_rot[,1] <- if_else(xy_com_m_rot[,1] > pi, xy_com_m_rot[,1] - 2*pi, xy_com_m_rot[,1])
  xy_com_m_rot <- xy_com_m_rot/pi*180
  
  # normalize
  x_com_m_rot_norm <- (xy_com_m_rot[,1] - min(xy_com_m_rot[,1])) / diff(range(xy_com_m_rot[,1]))
  
  for (j in 1:65) {
    # arc length
    oj <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[j,])), decreasing = F)

    o_rot <- order(((x_com_m_rot_norm - x_com_m_rot_norm[j])^2)[oj], decreasing = F)

    N_swap_rot[j,k] <- bubble_sort_swapCount(o_rot)
  }
}

swapMax <- 64*(64-1)/2
N_swap_rot_index <- 1 - 2*N_swap_rot/swapMax


df <- data.frame(cbind(seq(1,65), N_swap_rot_index)) 
colnames(df) <- c('neu', seq(0,170,by = 10))
dfm <- melt(df, id = 'neu')

dev.new()
# pdf(file = "index rotation.pdf", width = 10, height = 4,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  # geom_boxplot() +
  # geom_point() +
  # geom_jitter(colour = 'black', size = 4, height = 0, width = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', size = 5,position = position_nudge(x = 0, y = 0) ) +
  coord_cartesian(ylim = c(-.5, 1)) +
  theme_minimal_grid() +
  labs(title = "index rotation LO", x='Angle', y='RI') 

dev.off()



# plot 2d
# dev.new()
# plot(xy_com_m_rot)
# points(xy_com_m_rot[c(1,2),], cex = 2, col ='red', pch = 16)
# points(xy_com_m_rot[c(3,4),], cex = 2, col ='blue', pch = 16)

# PLOT 3d eg. neuron
which.min(cable_glo_min10_x_norm)

# # -- example neuron
# nopen3d()
# par3d('windowRect' = c(100,100,1100,1100))
# j <- 48
# plot3d(LC6_glo_rs[[j]], col = col4[4], lwd = 3, WithNodes = F)
# xyz <- xyzmatrix(LC6_glo_rs[[j]]$d)
# xyz10 <- xyz[xyz[,1] < quantile(xyz[,1], 0.1),]
# points3d(matrix(colMeans(xyz10),ncol=3), size = 30, col = "pink")
# j <- 3
# plot3d(LC6_glo_rs[[j]], col = col4[3], lwd = 3, WithNodes = F)
# xyz <- xyzmatrix(LC6_glo_rs[[j]]$d)
# xyz10 <- xyz[xyz[,1] < quantile(xyz[,1], 0.1),]
# points3d(matrix(colMeans(xyz10),ncol=3), size = 30, col = "cyan")
# 
# rgl.viewpoint(userMatrix = rotationMatrix(30/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0))
# rgl.pop(type = "light")
# rgl.light(theta = 30, phi = 30)
# 
# shade3d(glo.msh, alpha = 0.1)
# title3d('min 10 eg')


# -- plot regression
dev.new()
# par(mfrow=c(4,2))

# pdf(file = "median reg 2.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

df <- cbind.data.frame(com_rot_norm, cable_glo_med_x_norm)
colnames(df) <- c('x','y')
linReg <- lm(y ~ x, data = df)
df$fitted <- linReg$fitted.values
pval <- summary(linReg)$coefficients['x', 'Pr(>|t|)'] #p-value
rsq <- summary(linReg)$r.squared # R^2
plot(df[,1:2], pch=16, cex=2, col='black', xaxt='none', yaxt='n',
     # xlim = c(0,.6), ylim = c(0.3,0.8), 
     xlim = c(-0.5,0.5), ylim = c(0.3,0.8), 
     # xlim = c(-1,1), ylim = c(0,1), 
     cex.axis=2, cex.lab=2, cex.main=2, 
     main = "lin reg", xlab = "position in LO", ylab = "position in glo")
# axis(1, at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5), las=1) 
# axis(2, at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5), las=2) 
# plot(df[,1:2], pch=16, cex=2, col=4,
#      xlim = c(20,250), ylim = c(320000,430000), main = j,
#      cex.axis=2, cex.lab=2, cex.main=2, xaxt='n',yaxt='n', ann=F)
lines(df[,c(1,3)], lwd=2, col='black')
text(0, 0.8, labels = bquote(paste("slope = ", .(signif(linReg$coefficients[2],2)))),
     adj = 0, srt = 0, cex = 1)
# text(100, 420000, labels = bquote(paste("p-value = ", .(signif(pval,2)))),
#      adj = 0, srt = 0, cex = 1)
# text(100, 420000, labels = bquote(paste("p-value = ", .(signif(pval,2)))),
#      adj = 0, srt = 0, cex = 1)

dev.off()

df <- cbind.data.frame(xy_com_m_rot[,1], cable_glo_min10_x)
colnames(df) <- c('x','y')
linReg <- lm(y ~ x, data = df)
df$fitted <- linReg$fitted.values
pval <- summary(linReg)$coefficients['x', 'Pr(>|t|)'] #p-value
rsq <- summary(linReg)$r.squared # R^2
points(df[,1:2], pch=16, cex=2, col=5,
     xlim = c(20,250), ylim = c(320000,430000), main = j,
     cex.axis=2, cex.lab=2, cex.main=2,xaxt='n',yaxt='n', ann=F)
lines(df[,c(1,3)], lwd=2, col='black')
text(100, 350000, labels = bquote(paste("p-value = ", .(signif(pval,2)))),
     adj = 0, srt = 0, cex = 1)
title('median and min 10 reg')

# Figure 6C, 2d LO line, together with projection onto x-/y-axis
windows(record = F, width = 10.5, height = 10.5)

# pdf(file = "proj and best fit.pdf", width = 10.5, height = 10.5,pointsize=12,family="Helvetica", useDingbats = F)

bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
bd_theta <- seq(1, 180, by = 1)
xy_bd <- matrix(ncol = 2)
bd_grid <- expand.grid(bd_phi, bd_theta)
# range(bd_grid$Var1)
plot(bd_grid, xlim = c(-20,215), ylim = c(190,-50), type = "n", axes = FALSE, ann = F)
i1i <- 1; i2i <- 1; i3i <- 1; i4i <- 1
for (j in 1:length(xy_poly)) {
  xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
  if (j %in% i1) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[1], cex = 3, pch = 16)
  }
  else if (j %in% i2) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[2], cex = 3, pch = 16)
  }
  else if (j %in% i3) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[3], cex = 3, pch = 16)
  }
  else if (j %in% i4) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[4], cex = 3, pch = 16)
  }
  else {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="grey", cex = 2, pch = 16)
    # points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="grey", type = '.') 
  }
}
xy_bd <- xy_bd[-1,]
hpts_2 <- chull(xy_bd)
hpts_2 <- c(hpts_2, hpts_2[1])
xy_bd_chull <- xy_bd[hpts_2,] # hull edge points
polygon(xy_bd_chull)
# lines(rbind(c(-11,180), c(-2,180)), lwd = 1)
# text(-5, 180, labels = expression(paste("9",degree)), pos = 1, offset = 0.3)
# lines(rbind(c(-11,180), c(-11,171)), lwd = 1)
# text(-11, 175, labels = expression(paste("9",degree)), pos = 2, offset = 0.2)
# lines(rbind(c(0,0), c(0,180)), lwd = 1) 
# text(0, -5, labels = "front", pos = 1, offset = 0)
# lines(rbind(c(90,0), c(90,180)), lwd = 1)
# text(90, -5, labels = expression(paste("side 90",degree)), pos = 1, offset = 0)
# lines(rbind(c(-12,90), c(162,90)), lwd = 1)
# text(-17, 90, labels = "equator", pos = 1, offset = 0, srt = 90)

# add projection on x and y
lines(rbind(c(-10,190), c(160,190)), lwd = 1)
points(cbind(xy_com_m[i1,1], 190), pch = 16, col = col4[1], cex = 2)
points(cbind(xy_com_m[i2,1], 190), pch = 16, col = col4[2], cex = 2 )
points(cbind(xy_com_m[i3,1], 190), pch = 16, col = col4[3], cex = 2 )
points(cbind(xy_com_m[i4,1], 190), pch = 16, col = col4[4], cex = 2 )
points(cbind(xy_com_m[-c(i1,i2,i3,i4),1], 190), pch = 16, col = 'grey', cex = 1 )
lines(rbind(c(-20,0), c(-20,180)), lwd = 1)
points(cbind(-20, xy_com_m[i1,2]), pch = 16, col = col4[1], cex = 2)
points(cbind(-20, xy_com_m[i2,2]), pch = 16, col = col4[2], cex = 2 )
points(cbind(-20, xy_com_m[i3,2]), pch = 16, col = col4[3], cex = 2 )
points(cbind(-20, xy_com_m[i4,2]), pch = 16, col = col4[4], cex = 2 )
points(cbind(-20, xy_com_m[-c(i1,i2,i3,i4),2]), pch = 16, col = 'grey', cex = 1 )


# add rotation line for min 10 projection
j <- 31
ang <- j*10/180*pi
M_rot <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol = 2, byrow = T)

xy_com_m_rot <- xy_com_m %*% M_rot
xy_com_m_rot_inv <- cbind(xy_com_m_rot[,1],0) %*% t(M_rot)
line_rot <- cbind(c(20,220),0) %*% t(M_rot)
# points(seq(10,140),  - seq(10,140)*tan(ang), type = 'l', lwd = 2, col='darkgreen')

points(sweep(xy_com_m_rot_inv[i1,],2,c(70,-70),'+'),  pch = 16, col = col4[1], cex = 2 )
points(sweep(xy_com_m_rot_inv[i2,],2,c(70,-70),'+'),  pch = 16, col = col4[2], cex = 2 )
points(sweep(xy_com_m_rot_inv[i3,],2,c(70,-70),'+'),  pch = 16, col = col4[3], cex = 2 )
points(sweep(xy_com_m_rot_inv[i4,],2,c(70,-70),'+'),  pch = 16, col = col4[4], cex = 2 )
points(sweep(xy_com_m_rot_inv[-c(i1,i2,i3,i4),],2,c(70,-70),'+'), pch = 16, col = 'grey', cex = 1 )

points(sweep(line_rot,2,c(70,-70),'+'), type = 'l', lwd = 2, col='black')

title( "proj and best fit.pdf", )

dev.off()

# # alpha palette
# mycols <- adjustcolor(palette(), alpha.f = 0.3)
# opal <- palette(mycols)



# reti_stats ------------------------------------------------------------------------------------------------------


xy_com_m_tp <- xy_com_m[,c(2,1)] / 180 * pi
reti_stats <- matrix(ncol = 4, nrow = 65)
for (j in 1:nrow(reti_stats)) {
  oj <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[j,])))
  # oj <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[j,], '-')^2), decreasing = F)
  cable_avg <- order(rowSums(sweep(cable_glo_avg, 2, cable_glo_avg[j,], '-')^2)[oj]) 
  cable_med <- order(((cable_glo_med_x - cable_glo_med_x[j])^2)[oj]) 
  syn_avg <- order(rowSums(sweep(syn_xyz_avg, 2, syn_xyz_avg[j,], '-')^2)[oj]) 
  syn_med <- order(((syn_xyz_med_x - syn_xyz_med_x[j])^2)[oj]) 
  
  # cable_avg_1 <- order(rowSums(sweep(com_cable[[2]], 2, com_cable[[2]][j,], '-')^2)[oj]) 
  # cable_avg_2 <- order(rowSums(sweep(com_cable[[3]], 2, com_cable[[3]][j,], '-')^2)[oj]) 
  
  reti_stats[j,] <- c(
    bubble_sort_swapCount(cable_avg),
    bubble_sort_swapCount(cable_med),
    bubble_sort_swapCount(syn_avg),
    bubble_sort_swapCount(syn_med)
    #, bubble_sort_swapCount(cable_avg_1),
    # bubble_sort_swapCount(cable_avg_2)
    )
}

# PLOT
swapMax <- 64*(64-1)/2
reti_stats_index <- 1 - 2*reti_stats/swapMax

df <- data.frame(cbind(seq(1,65), reti_stats_index))
colnames(df) <- c('neu', 'cable_avg', 'cable_med','syn_avg', 'syn_med')
# colnames(df) <- c('neu', 'cable_avg', 'cable_med','syn_avg', 'syn_med','sec1', 'sec2')
# colnames(df) <- c('neu','sec1', 'sec2')
dfm <- melt(df, id = 'neu')
dev.new()
ggplot(dfm, aes(x=neu, y=value, shape=variable) ) + 
  geom_point(size=3)+
  # scale_color_manual(labels = c("T999", "T888"), values = c("blue", "red")) +
  # coord_cartesian(xlim = c(1, 64), ylim = c(0,2080)) +
  coord_cartesian(xlim = c(1, 64), ylim = c(-1,1)) +
  labs(title = "num of adjacent transpositions, arc length",
       x='neuron', y='num of swaps',
       shape = 'quantity') 
  

# re-ordering of nodes/ syn distances in glo sectors-----------------------------------------------------------------


# -- pre-synapse in glo
LC6_pre_glo <- list() #pre-synapse xyz 
LC6_post_glo <- list() #post-synapse xyz 
for (j in 1:length(LC6)) {
  conn <- connectors(LC6[[j]])
  conn_pre <- conn[conn$prepost== 0, c('x','y','z')]
  inglo <- rowSums(sweep(conn_pre,2,c(a,b,c), '*')) + d < 0
  LC6_pre_glo[[j]] <- conn_pre[inglo,]
  
  # ii_pre <-  tar$d$PointNo %in% conn[conn$prepost==0,]$treenode_id
  # conn_pre <- xyzmatrix(tar$d[ii_pre,])
  
  conn_post <- conn[conn$prepost== 1, c('x','y','z')]
  inglo <- rowSums(sweep(conn_post,2,c(a,b,c), '*')) + d < 0
  LC6_post_glo[[j]] <- conn_post[inglo,]
}

LC6_pre_glo_all <- do.call(rbind.data.frame, LC6_pre_glo)

# check numbers within each compartment
glo_div_pre <- quantile(LC6_pre_glo_all[,1], probs = seq(0,1,length.out = 11)) 
N_pre_comp <- matrix(ncol = 10, nrow = 65)
for (j in 1:65) {
  xx <- LC6_pre_glo[[j]][,1]
  for (k in 1:10) {
    # N_pre_comp[j,k] <- sum(xx > glo_div_pre[k] & xx < glo_div_pre[k+1])
    N_pre_comp[j,k] <- sum(xx > glu_div[k] & xx < glu_div[k+1])
  }
}

colSums(N_pre_comp == 0)

# try 4, 5, 6, 8 and combined
#manual add each
N_swap_pre_ls <- list()
N_swap_pre_3_ls <- list()

for (s in c(4,5,6,8)) {
  # s <- 8
  LC6_pre_glo_comp <- list() # selected compartment
  jxyz_3 <- list()
  for (j in 1:65) {
    xx <- LC6_pre_glo[[j]][,1]
    # ii <- xx > glo_div_pre[s] & xx < glo_div_pre[s+1]
    ii <- xx > glu_div[s] & xx < glu_div[s+1] # USE old
    LC6_pre_glo_comp[[j]] <- LC6_pre_glo[[j]][ii,]
    # i3 <- xx > glu_div[s-1] & xx < glu_div[s+1+1]
    # jxyz_3[[j]] <- LC6_pre_glo[[j]][ii,] #check nb compartment
  }
  
  N_swap_pre <- matrix(ncol = 2, nrow = 65)
  for (j in 1:65) {
    # oj <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[j,], '-')^2), decreasing = F) 
    oj <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[j,])), decreasing = F)
    jxyz <- LC6_pre_glo_comp[[j]] 
    jswap <- c()
    for (k in 1:nrow(jxyz)) {
      dd <- c() #distance
      for (m in 1:65) {
        # mdd <- min(rowSums(sweep(jxyz_3[[m]], 2, as.numeric(jxyz[k,]), '-')^2)) #shortes to a neuron
        mdd <- min(rowSums(sweep(LC6_pre_glo_comp[[m]], 2, as.numeric(jxyz[k,]), '-')^2)) #shortes to a neuron
        dd <- c(dd, mdd)
      }
      
      order_k <- order(dd[oj])  # re-order
      # order_k <- order(dd[sample(65)])  # random
      
      jswap <- c(jswap, bubble_sort_swapCount(order_k)) #num of swaps
    }
    N_swap_pre[j,] <- c(mean(jswap), sd(jswap))
  }
  
  
  N_swap_pre_ls[[s]] <- N_swap_pre  
  # N_swap_pre_3_ls[[s]] <- N_swap_pre
}

# combine 4,5,6
LC6_pre_glo_comp <- list() # selected compartment
for (j in 1:65) {
  xx <- LC6_pre_glo[[j]][,1]
  ii <- xx > glu_div[4] & xx < glu_div[7]
  LC6_pre_glo_comp[[j]] <- LC6_pre_glo[[j]][ii,]
}
N_swap_pre_456 <- matrix(ncol = 2, nrow = 65)
N_swap_pre_ran <- matrix(ncol = 2, nrow = 65)
for (j in 1:65) {
  oj <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[j,])), decreasing = F)
  jxyz <- LC6_pre_glo_comp[[j]] 
  jswap <- c()
  jswap_ran <- c()
  for (k in 1:nrow(jxyz)) {
    dd <- c() #distance
    for (m in 1:65) {
      mdd <- min(rowSums(sweep(LC6_pre_glo_comp[[m]], 2, as.numeric(jxyz[k,]), '-')^2)) #shortes to a neuron
      dd <- c(dd, mdd)
    }
    
    order_k <- order(dd[oj])  # re-order
    jswap <- c(jswap, bubble_sort_swapCount(order_k)) #num of swaps
    
    order_k_ran <- order(dd[sample(65)])  # random
    jswap_ran <- c(jswap_ran, bubble_sort_swapCount(order_k_ran)) #num of swaps
    
  }
  N_swap_pre_456[j,] <- c(mean(jswap), sd(jswap))
  N_swap_pre_ran[j,] <- c(mean(jswap_ran), sd(jswap_ran))
}



# # a random one
# s <- 8
# LC6_pre_glo_comp <- list() # selected compartment
# jxyz_3 <- list()
# for (j in 1:65) {
#   xx <- LC6_pre_glo[[j]][,1]
#   ii <- xx > glo_div_pre[s] & xx < glo_div_pre[s+1]
#   LC6_pre_glo_comp[[j]] <- LC6_pre_glo[[j]][ii,]
#   i3 <- xx > glo_div_pre[s-1] & xx < glo_div_pre[s+1+1]
#   jxyz_3[[j]] <- LC6_pre_glo[[j]][ii,] #check nb compartment
# }
# 
# N_swap_pre <- matrix(ncol = 2, nrow = 65)
# for (j in 1:65) {
#   # oj <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[j,], '-')^2), decreasing = F) 
#   oj <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[j,])), decreasing = F)
#   jxyz <- LC6_pre_glo_comp[[j]] 
#   jswap <- c()
#   for (k in 1:nrow(jxyz)) {
#     dd <- c() #distance
#     for (m in 1:65) {
#       mdd <- min(rowSums(sweep(jxyz_3[[m]], 2, as.numeric(jxyz[k,]), '-')^2)) #shortes to a neuron
#       # mdd <- min(rowSums(sweep(LC6_pre_glo_comp[[m]], 2, as.numeric(jxyz[k,]), '-')^2)) #shortes to a neuron
#       dd <- c(dd, mdd)
#     }
#     
#     # order_k <- order(dd[oj])  # re-order
#     order_k <- order(dd[sample(65)])  # random
#     
#     jswap <- c(jswap, bubble_sort_swapCount(order_k)) #num of swaps
#   }
#   N_swap_pre[j,] <- c(mean(jswap), sd(jswap))
# }
# 
# N_swap_pre_ran <- N_swap_pre  

# df <- data.frame(cbind(seq(1,65),
#                        1 - 2*N_swap_pre_ls[[1]][,1]/swapMax,
#                        1 - 2*N_swap_pre_ls[[2]][,1]/swapMax,
#                        1 - 2*N_swap_pre_ls[[3]][,1]/swapMax))
# colnames(df) <- c('neu', 'sec4', 'sec5', 'sec6')
# dfm <- melt(df, id = 'neu')
# dev.new()
# ggplot(dfm, aes(x=neu, y=value, shape=variable) ) + 
#   geom_point(size=3)+
#   # scale_color_manual(labels = c("T999", "T888"), values = c("blue", "red")) +
#   # coord_cartesian(xlim = c(1, 64), ylim = c(0,2080)) +
#   coord_cartesian(xlim = c(1, 64), ylim = c(-1,1)) +
#   labs(title = "num of adjacent transpositions for 3 sections, pre synp",
#        x='neuron', y='num of swaps',
#        shape = 'quantity') 



# arti mapping --------------------------------------------------------------------------------------------

# position after rotation:  pt1 <- pt0 %*% t(R_mat)
quaternion3D <- function(vr, ang){
  ang <- ang / 180 * pi
  vr <- vr / sqrt(sum(vr^2))
  
  qr <- cos(ang/2)
  qi <- vr[1]*sin(ang/2)
  qj <- vr[2]*sin(ang/2)
  qk <- vr[3]*sin(ang/2)
  
  R <- matrix(c(
    1-2*(qj^2+qk^2), 2*(qi*qj-qk*qr), 2*(qi*qk+qj*qr),
    2*(qi*qj+qk*qr), 1-2*(qi^2+qk^2), 2*(qj*qk-qi*qr),
    2*(qi*qk-qj*qr), 2*(qj*qk+qi*qr), 1-2*(qi^2+qj^2)),
    ncol = 3, byrow = T)
  
  return(R)
}

# theta is from z-axis
cart2sphZ <- function(xyz){ 
  r <- sqrt(rowSums(xyz^2))
  x <- xyz[,1] / r
  y <- xyz[,2] / r
  z <- xyz[,3] / r
  theta <- acos(z)
  phi <- 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))
  
  return(cbind(r, theta, phi)) #theta [0,pi], phi [0, 2*pi]
}

sph2cartZ <- function(rtp){
  r <- rtp[,1]
  t <- rtp[,2]
  p <- rtp[,3]
  x <- r * cos(p) * sin(t)
  y <- r * sin(p) * sin(t)
  z <- r * cos(t)
  
  return(cbind(x,y,z))
}

# #  cylinder loop azim
# ae <- xy_com_m
# ae[,1] <- ae[,1] / diff(range(ae[,1])) * 360 #exppand to 360 deg
# tp <- ae[,c(2,1)] / 180 * pi #theta phi
# # xyz <- sph2cart(cbind(tp,1))
# xyz <- sph2cartZ(cbind(1,tp))
# xyz[,3] <- xyz[,3] * 4
# 
# # cylinder loop after ccw 90deg rotation
# R <- quaternion3D(c(1,1,0), -90)
# ae <- xy_com_m
# tp <- ae[,c(2,1)] / 180 * pi #theta phi
# xyz2 <- sph2cartZ(cbind(1,tp))
# xyz2 <- xyz2 %*% t(R)
# tp <- cart2sphZ(xyz2)[,2:3]
# tp[,2] <- if_else(tp[,2]<=pi, tp[,2], tp[,2]-2*pi )
# tp[,2] <- tp[,2] / diff(range(tp[,2])) * 2.3*pi #exppand to 360 deg
# xyz2 <- sph2cartZ(cbind(1,tp))
# xyz2[,3] <- xyz2[,3] * 4

# # best fit
# j <- 31 # ccw 310 deg
# ang <- j*10/180*pi
# M_rot <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol = 2, byrow = T)
# xy_com_m_rot <- xy_com_m %*% M_rot

# check projection, rotate 90 and use azim coord
# sum(range(xy_com_m[,1]))/2
# sum(range(xy_com_m[,2]))/2
vr <- c(cos(sum(range(xy_com_m[,1]))/2/180*pi), sin(sum(range(xy_com_m[,1]))/2/180*pi),0)
R_mat <- quaternion3D(vr, ang = 90)

xy_com_m_rot_90 <- sph2cartZ(cbind(1, xy_com_m_tp)) %*% t(R_mat)
xy_com_m_rot_90 <- cart2sphZ(xy_com_m_rot_90)[, c(3,2)]
xy_com_m_rot_90[,1] <- if_else(xy_com_m_rot_90[,1] > pi, xy_com_m_rot_90[,1] - 2*pi, xy_com_m_rot_90[,1])



j <- 3 #21
ang <- j*10
vr <- c(cos(sum(range(xy_com_m[,1]))/2/180*pi), sin(sum(range(xy_com_m[,1]))/2/180*pi),0)
R_mat <- quaternion3D(vr, ang)

xy_com_m_rot <- sph2cartZ(cbind(1, xy_com_m_tp)) %*% t(R_mat) #ccw on S2
xy_com_m_rot <- cart2sphZ(xy_com_m_rot)[, c(3,2)]
xy_com_m_rot[,1] <- if_else(xy_com_m_rot[,1] > pi, xy_com_m_rot[,1] - 2*pi, xy_com_m_rot[,1])
xy_com_m_rot <- xy_com_m_rot/pi*180


# - plot 4 groups after rotation
windows(record = F, width = 8, height = 8)

# pdf(file = "4 group LO.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

xy_com_tmp <- xy_com_m_rot
colnames(xy_com_tmp) <- c('phi_deg', 'theta_deg')

# bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
# bd_theta <- seq(1, 180, by = 1)
# xy_bd <- matrix(ncol = 2)
# bd_grid <- expand.grid(bd_phi, bd_theta)
plot(xy_com_tmp, ylim = rev(range(xy_com_tmp[,2])), type = "n")#, axes = FALSE, ann = F)
for (j in 1:length(xy_poly)) {
  # xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
  if (j %in% i1) {
    points(xy_com_tmp[j,'phi_deg'], xy_com_tmp[j,'theta_deg'], col=col4_a1[1], cex = 3, pch = 16)
  }
  else if (j %in% i2) {
    points(xy_com_tmp[j,'phi_deg'], xy_com_tmp[j,'theta_deg'], col=col4_a1[2], cex = 3, pch = 16)
  }
  else if (j %in% i3) {
    points(xy_com_tmp[j,'phi_deg'], xy_com_tmp[j,'theta_deg'], col=col4_a1[3], cex = 3, pch = 16)
  }
  else if (j %in% i4) {
    points(xy_com_tmp[j,'phi_deg'], xy_com_tmp[j,'theta_deg'], col=col4_a1[4], cex = 3, pch = 16)
  }
  else {
    points(xy_com_tmp[j,'phi_deg'], xy_com_tmp[j,'theta_deg'], col="grey", cex = 2, pch = 16) 
  }
  text(xy_com_tmp[j,'phi_deg'], xy_com_tmp[j,'theta_deg'], labels = j, pos = 1)
}
# xy_bd <- xy_bd[-1,]
# hpts_2 <- chull(xy_bd)
# hpts_2 <- c(hpts_2, hpts_2[1])
# xy_bd_chull <- xy_bd[hpts_2,] # hull edge points
# polygon(xy_bd_chull)
# lines(rbind(c(-11,180), c(-2,180)), lwd = 1)
# text(-5, 180, labels = expression(paste("9",degree)), pos = 1, offset = 0.3)
# lines(rbind(c(-11,180), c(-11,171)), lwd = 1)
# text(-11, 175, labels = expression(paste("9",degree)), pos = 2, offset = 0.2)
# lines(rbind(c(0,0), c(0,180)), lwd = 1) 
# text(0, -5, labels = "front", pos = 1, offset = 0)
# lines(rbind(c(90,0), c(90,180)), lwd = 1)
# text(90, -5, labels = expression(paste("side 90",degree)), pos = 1, offset = 0)
# lines(rbind(c(-12,90), c(162,90)), lwd = 1)
# text(-17, 90, labels = "equator", pos = 1, offset = 0, srt = 90)

dev.off()


# - arti index
N_swap_man <- matrix(ncol = 4, nrow = 65)
for (j in 1:65) {
  # arc length
  oj <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[j,])), decreasing = F)
  
  # - along elev or azim
  o_phi <- order(((xy_com_m[,1] - xy_com_m[j,1])^2)[oj]) 
  o_theta <- order(((xy_com_m[,2] - xy_com_m[j,2])^2)[oj]) 
  o_rot <- order(((xy_com_m_rot[,1] - xy_com_m_rot[j,1])^2)[oj]) 
  o_rot_90 <- order(((xy_com_m_rot_90[,1] - xy_com_m_rot_90[j,1])^2)[oj]) 
  
  # o_xyz <- order(rowSums(sweep(xyz, 2, xyz[j,], '-')^2)[oj]) 
  # o_xyz2 <- order(rowSums(sweep(xyz2, 2, xyz2[j,], '-')^2)[oj]) 
  
  N_swap_man[j,] <- c(
    bubble_sort_swapCount(o_phi)
    ,    bubble_sort_swapCount(o_theta)
    ,    bubble_sort_swapCount(o_rot)
    ,    bubble_sort_swapCount(o_rot_90)
    # ,    bubble_sort_swapCount(o_xyz)
    # ,    bubble_sort_swapCount(o_xyz2)
  )
}

swapMax <- 64*(64-1)/2
N_swap_man_index <- 1 - 2*N_swap_man/swapMax

df <- data.frame(cbind(seq(1,65), N_swap_man_index))
colnames(df) <- c('neu', 'phi', 'theta','cylin z', 'cylin 60')
# colnames(df) <- c('neu','sec1', 'sec2')
dfm <- melt(df, id = 'neu')
dev.new()
ggplot(dfm, aes(x=neu, y=value, shape=variable) ) + 
  geom_point(size=3)+
  # scale_color_manual(labels = c("T999", "T888"), values = c("blue", "red")) +
  # coord_cartesian(xlim = c(1, 64), ylim = c(0,2080)) +
  coord_cartesian(xlim = c(1, 64), ylim = c(-1,1)) +
  labs(title = "num of adjacent transpositions, arc length, manual map",
       x='neuron', y='num of swaps',
       shape = 'quantity') 


# - index
N_swap_1D_rot <- matrix(ncol = 18, nrow = 65)
for (k in 1:18) {
  ang <- (k-1)*10
  R_mat <- quaternion3D(vr, ang)
  
  xy_com_m_rot <- sph2cartZ(cbind(1, xy_com_m_tp)) %*% t(R_mat) #ccw on S2
  xy_com_m_rot <- cart2sphZ(xy_com_m_rot)[, c(3,2)]
  xy_com_m_rot[,1] <- if_else(xy_com_m_rot[,1] > pi, xy_com_m_rot[,1] - 2*pi, xy_com_m_rot[,1])
  xy_com_m_rot <- xy_com_m_rot/pi*180
  
  # # normalize
  # x_com_m_rot_norm <- (xy_com_m_rot[,1] - min(xy_com_m_rot[,1])) / diff(range(xy_com_m_rot[,1]))
  
  for (j in 1:65) {
    o_phi <- order((xy_com_m_rot[,1] - xy_com_m_rot[j,1])^2) 
    
    # o_min10_x <- order(((cable_glo_min10_x - cable_glo_min10_x[j])^2)[o_phi]) 
    o_med_x <- order(((cable_glo_med_x - cable_glo_med_x[j])^2)[o_phi]) 
    
    N_swap_1D_rot[j,k] <- bubble_sort_swapCount(o_med_x)
  }
}

swapMax <- 64*(64-1)/2
N_swap_1D_rot_index <- 1 - 2*N_swap_1D_rot/swapMax


df <- data.frame(cbind(seq(1,65), N_swap_1D_rot_index)) 
colnames(df) <- c('neu', seq(0,170,by = 10))
dfm <- melt(df, id = 'neu')

dev.new()
# pdf(file = "index rotation 1D vs glo med.pdf", width = 10, height = 4,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  # geom_boxplot() +
  # geom_point() +
  # geom_jitter(colour = 'black', size = 4, height = 0, width = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', size = 5,position = position_nudge(x = 0, y = 0) ) +
  coord_cartesian(ylim = c(-.5, 1)) +
  theme_minimal_grid() +
  labs(title = "index 1D vs glo median sweep", x='Angle', y='RI') 

dev.off()



# - index ordering along 1D subspace VS glo axis
N_swap_1D <- matrix(ncol = 6, nrow = 65)

j <- 3 #21
ang <- j*10
vr <- c(cos(sum(range(xy_com_m[,1]))/2/180*pi), sin(sum(range(xy_com_m[,1]))/2/180*pi),0)
R_mat <- quaternion3D(vr, ang)

xy_com_m_max <- sph2cartZ(cbind(1, xy_com_m_tp)) %*% t(R_mat) #ccw on S2
xy_com_m_max <- cart2sphZ(xy_com_m_max)[, c(3,2)]
xy_com_m_max[,1] <- if_else(xy_com_m_max[,1] > pi, xy_com_m_max[,1] - 2*pi, xy_com_m_max[,1])
xy_com_m_max <- xy_com_m_max/pi*180

# ojx <- order(xy_com_m[,1], decreasing = F)
# ojy <- order(xy_com_m[,2], decreasing = F)
# ojbf <- order(xy_com_m_rot[,1], decreasing = F)
for (j in 1:65) {
  
  o_phi <- order((xy_com_m[,1] - xy_com_m[j,1])^2) 
  o_rot <- order((xy_com_m_max[,1] - xy_com_m_max[j,1])^2) 
  o_rot_90 <- order((xy_com_m_rot_90[,1] - xy_com_m_rot_90[j,1])^2)
  
  o_med_x <- order(((cable_glo_med_x - cable_glo_med_x[j])^2)[o_phi]) 
  o_med_y <- order(((cable_glo_med_x - cable_glo_med_x[j])^2)[o_rot_90]) 
  o_med_bf <- order(((cable_glo_med_x - cable_glo_med_x[j])^2)[o_rot]) 
  
  o_min10_x <- order(((cable_glo_min10_x - cable_glo_min10_x[j])^2)[o_phi]) 
  o_min10_y <- order(((cable_glo_min10_x - cable_glo_min10_x[j])^2)[o_rot_90]) 
  o_min10_bf <- order(((cable_glo_min10_x - cable_glo_min10_x[j])^2)[o_rot]) 
  
  # o_med_x <- order(((cable_glo_med_x - cable_glo_med_x[j])^2)[ojx]) 
  # o_med_y <- order(((cable_glo_med_x - cable_glo_med_x[j])^2)[ojy]) 
  # o_med_bf <- order(((cable_glo_med_x - cable_glo_med_x[j])^2)[ojbf]) 
  # 
  # o_min10_x <- order(((cable_glo_min10_x - cable_glo_min10_x[j])^2)[ojx]) 
  # o_min10_y <- order(((cable_glo_min10_x - cable_glo_min10_x[j])^2)[ojy]) 
  # o_min10_bf <- order(((cable_glo_min10_x - cable_glo_min10_x[j])^2)[ojbf]) 
  
  N_swap_1D[j,] <- c(
    bubble_sort_swapCount(o_med_x),
    bubble_sort_swapCount(o_med_y),
    bubble_sort_swapCount(o_med_bf),
    bubble_sort_swapCount(o_min10_x),
    bubble_sort_swapCount(o_min10_y),
    bubble_sort_swapCount(o_min10_bf)
  )
}

swapMax <- 64*(64-1)/2
N_swap_1D_index <- 1 - 2*N_swap_1D/swapMax

# # - test rotation angle
# N_swap_cyn <- matrix(ncol = 6, nrow = 65)
# for (k in 1:6) {
#   R <- quaternion3D(c(1,1,0), -30*(k-1))
#   ae <- xy_com_m
#   tp <- ae[,c(2,1)] / 180 * pi #theta phi
#   xyz2 <- sph2cartZ(cbind(1,tp))
#   xyz2 <- xyz2 %*% t(R)
#   tp <- cart2sphZ(xyz2)[,2:3]
#   tp[,2] <- if_else(tp[,2]<=pi, tp[,2], tp[,2]-2*pi )
#   tp[,2] <- tp[,2] / diff(range(tp[,2])) * 2.3*pi #exppand to 360 deg
#   xyz2 <- sph2cartZ(cbind(1,tp))
#   xyz2[,3] <- xyz2[,3] * 2
#   for (j in 1:65) {
#     # arc length
#     oj <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[j,])), decreasing = F)
#     o_xyz2 <- order(rowSums(sweep(xyz2, 2, xyz2[j,], '-')^2)[oj]) 
#     
#     N_swap_cyn[j,k] <- bubble_sort_swapCount(o_xyz2)
#   }
# }
# # -violin
# swapMax <- 65*(65-1)/2
# N_swap_man_index <- 1 - 2*N_swap_cyn/swapMax
# df <- data.frame(cbind(seq(1,65), 
#                        N_swap_man_index)
# )
# colnames(df) <- c('neu', 'cylin z', 'cylin 30', 'cyl60', 'cyl90','cyl120','cyl150')
# 
# dfm <- melt(df, id = 'neu')
# dev.new()
# ggplot(dfm, aes(factor(variable), value) ) + 
#   geom_violin(scale = 'width') +
#   geom_jitter(height = 0, width = 0.05) +
#   ylim(-1,1) +
#   labs(title = "violin plot",
#        x='metrics', y='index') 


# -- cartoon transformation
nopen3d()
# points3d(xyz[-c(i1,i2,i3,i4),], size = 10, col = 'grey')
points3d(xyz[-c(i1,i2,i3,i4),], size = 9, col = 'black')
points3d(xyz[i1,], size = 15, col = col4[1])
points3d(xyz[i2,], size = 15, col = col4[2])
points3d(xyz[i3,], size = 15, col = col4[3])
points3d(xyz[i4,], size = 15, col = col4[4])
spheres3d(0,0,0,1, col='grey', alpha=1)

shade3d(ellipse3d(diag(3),centre=c(0,0,0),scale=c(1,1,4)), col='grey')

# use ref points
xyz_a <- cbind(cos(seq(0,(2*pi-0.2),length.out = 6)), sin(seq(0,(2*pi-0.2),length.out = 6)), 0)
# getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
a_col <- rev(brewer.pal(9,'Blues'))
xyz_e <- cbind(0, sin(seq(0,pi,length.out = 6)), cos(seq(0,pi,length.out = 6)))
xyz_e[,3] <- xyz_e[,3] * 4
# getPalette <- colorRampPalette(brewer.pal(6, "Reds"))
e_col <- rev(brewer.pal(9,'Reds'))

nopen3d(userMatrix = M)
for (j in 1:6) {
  points3d(matrix(xyz_a[j,],ncol = 3), size = 15, col = a_col[j])
  points3d(matrix(xyz_e[j,],ncol = 3), size = 15, col = e_col[j])
}

# spheres3d(0,0,0,0.98, col='grey', alpha=0.5, lit = F)

shade3d(ellipse3d(diag(3)*0.12,centre=c(0,0,0),scale=c(1,1,4)), col='grey', alpha =.5, lit=F)



# ordering of 4 cross sections --------------------------------------------------------------------------------------------

swapMax <- 19*(19-1)/2

yz_cs <- list()
for (j in 1:4) {
  yz_cs[[j]] <- as.matrix(df_cable[((N_ind+1)*(c(i1,i2,i3,i4)-1)+j), c(1,2)])
}
 
rt_ind <- matrix(ncol = 4, nrow = 20)
for (j in 1:20) {
  oj <- order( rowSums(sweep(yz_cs[[1]], 2, yz_cs[[1]][j,], '-')^2))
  for (k in 1:4) {
    o_cs <- order(rowSums(sweep(yz_cs[[k]], 2, yz_cs[[k]][j,], '-')^2)[oj]) 
    rt_ind[j,k] <- 1 - 2 * bubble_sort_swapCount(o_cs) / swapMax
  }
} 

# plot
df <- as.data.frame(cbind(c(i1,i2,i3,i4), rt_ind) )
colnames(df) <- c('neu', '1','2','3','4')
dfm <- melt(df, id = 'neu')
dev.new()
# pdf(file = "crossSection index.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
ggplot(dfm, aes(factor(variable), value) ) + 
  # geom_boxplot() +
  # geom_point(colour = 'pink') +
  geom_jitter(height = 0, width = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar",
               colour = 'red', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'red', size = 5 ) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal_grid() +
  labs(title = "4 cross sections", x='position', y='retino-index') 
# dev.off()


# all 65
swapMax <- 64*(64-1)/2

yz_cs <- list()
for (j in 1:4) {
  yz_cs[[j]] <- as.matrix(df_cable[((N_ind+1)*(seq(1,65)-1)+j), c(1,2)])
}

rt65_ind <- matrix(ncol = 4, nrow = 65)
for (j in 1:65) {
  oj <- order( rowSums(sweep(yz_cs[[1]], 2, yz_cs[[1]][j,], '-')^2))
  for (k in 1:4) {
    o_cs <- order(rowSums(sweep(yz_cs[[k]], 2, yz_cs[[k]][j,], '-')^2)[oj]) 
    rt65_ind[j,k] <- 1 - 2 * bubble_sort_swapCount(o_cs) / swapMax
  }
} 

# plot
col_4_65 <- rep('grey', 65)
col_4_65 <- adjustcolor(col_4_65, alpha.f = 0.6)
col_4_65[c(i1,i2,i3,i4)] <- c(rep(col4[1],5),rep(col4[2],5),rep(col4[3],5),rep(col4[4],5))

df <- as.data.frame(cbind(seq(1,65), rt65_ind) )
colnames(df) <- c('neu', '1','2','3','4')
dfm <- melt(df, id = 'neu')
dfm$neu <- factor(dfm$neu)

dev.new()
# pdf(file = "4 crossSection index.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  geom_jitter(aes(colour=neu),size = 4,height = 0, width = 0.1) +
  scale_colour_manual(values = col_4_65, guide=FALSE) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', size = 5, position = position_nudge(x = 0.2, y = 0)) +
  coord_cartesian(ylim = c(-0.5, 1)) +
  theme_minimal_grid() +
  labs(title = "4 cross sections", x='position', y='RI') 
dev.off()


# ordering of 10 cross sections in glo--------------------------------------------------------------------------------------------

# TODO
# NA is least when order,
# when target is NA, index = 1



yz_cs_glo <- list()
for (j in 1:11) {
  yz_cs_glo[[j]] <- as.matrix(df_cable_sp[((N_ind_sp+1)*(c(i1,i2,i3,i4)-1)+j), c(1,2)])
}

xy_com_m_tp_20 <- xy_com_m_tp[c(i1,i2,i3,i4),]
rt_ind_glo <- matrix(NA, ncol = 11, nrow = 20)
for (j in 1:11) {
  yz <- yz_cs_glo[[j]] #y z in cable
  ii <- which(!is.na(yz[,1]))
  yz <- yz[ii,]
  tp <- xy_com_m_tp_20[ii,] #theta phi in LO
  swapMax <- (nrow(tp)-1)*(nrow(tp)-2)/2
  for (k in 1:nrow(tp)) {
    oj <- order(apply(tp, 1, function(x) arcLength(x, tp[k,])) )
    o_cs <- order(rowSums(sweep(yz, 2, yz[k,], '-')^2)[oj]) 
    rt_ind_glo[ii[k],j] <- 1 - 2 * bubble_sort_swapCount(o_cs) / swapMax
  }
}

# for (j in 1:20) {
#   oj <- order( rowSums(sweep(yz_cs_glo[[1]], 2, yz_cs_glo[[1]][j,], '-')^2))
#   for (k in 1:11) {
#     o_cs <- order(rowSums(sweep(yz_cs_glo[[k]], 2, yz_cs_glo[[k]][j,], '-')^2)[oj]) 
#     swapMax <- 19*(19-1)/2
#     rt_ind_glo[j,k] <- 1 - 2 * bubble_sort_swapCount(o_cs) / swapMax
#   }
# } 

# plot
df <- as.data.frame(cbind(c(i1,i2,i3,i4), rt_ind_glo) )
colnames(df) <- c('neu', seq(1,11))
dfm <- melt(df, id = 'neu')
dfm <- dfm[!is.na(dfm$value),]
# tmp <- dfm[21:220,]  # rid of NA / 1
# tmp <- tmp[!(tmp$value==1),]
# dfm <- rbind(dfm[1:20,], tmp)
dev.new()
# pdf(file = "crossSection index.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
ggplot(dfm, aes(factor(variable), value) ) + 
  # geom_boxplot() +
  # geom_point(colour = 'pink') +
  geom_jitter(height = 0, width = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar",
               colour = 'red', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'red', size = 5 ) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal_grid() +
  labs(title = "11 cross sections", x='position', y='retino-index') 
dev.off()



# - all 65

yz_cs_glo <- list()
for (j in 1:11) {
  yz_cs_glo[[j]] <- as.matrix(df_cable_sp[((N_ind_sp+1)*(seq(1,65)-1)+j), c(1,2)])
}

xy_com_m_tp_65 <- xy_com_m_tp
rt_ind_glo <- matrix(NA, ncol = 11, nrow = 65)
for (j in 1:11) {
  yz <- yz_cs_glo[[j]] #y z in cable
  ii <- which(!is.na(yz[,1]))
  yz <- yz[ii,]
  tp <- xy_com_m_tp_65[ii,] #theta phi in LO
  swapMax <- (nrow(tp)-1)*(nrow(tp)-2)/2
  for (k in 1:nrow(tp)) {
    oj <- order(apply(tp, 1, function(x) arcLength(x, tp[k,])) )
    o_cs <- order(rowSums(sweep(yz, 2, yz[k,], '-')^2)[oj]) 
    rt_ind_glo[ii[k],j] <- 1 - 2 * bubble_sort_swapCount(o_cs) / swapMax
  }
}


# yz_cs_glo <- list()
# for (j in 1:11) {
#   yz_cs_glo[[j]] <- as.matrix(df_cable_sp[((N_ind+1)*(seq(1,65)-1)+j), c(1,2)])
# }
# 
# rt65_ind_glo <- matrix(ncol = 11, nrow = 65)
# for (j in 1:65) {
#   oj <- order( rowSums(sweep(yz_cs_glo[[1]], 2, yz_cs_glo[[1]][j,], '-')^2))
#   for (k in 1:11) {
#     o_cs <- order(rowSums(sweep(yz_cs_glo[[k]], 2, yz_cs_glo[[k]][j,], '-')^2)[oj])
#     rt65_ind_glo[j,k] <- 1 - 2 * bubble_sort_swapCount(o_cs) / swapMax
#   }
# }

# plot
col_4_65 <- rep('grey', 65)
col_4_65 <- adjustcolor(col_4_65, alpha.f = 0.6)
col_4_65[c(i1,i2,i3,i4)] <- c(rep(col4[1],5),rep(col4[2],5),rep(col4[3],5),rep(col4[4],5))

df <- as.data.frame(cbind(seq(1,65), rt_ind_glo) )
colnames(df) <- c('neu', seq(1,11))
dfm <- melt(df, id = 'neu')
dfm$neu <- factor(dfm$neu)
dfm <- dfm[!is.na(dfm$value),]
# tmp <- dfm[66:715,] 
# tmp <- tmp[!(tmp$value==1),]
# dfm <- rbind(dfm[1:65,], tmp)

dev.new()
# pdf(file = "4 crossSection index.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  geom_jitter(aes(colour=neu),size = 4,height = 0, width = 0.1) +
  scale_colour_manual(values = col_4_65, guide=FALSE) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', size = 5, position = position_nudge(x = 0.2, y = 0)) +
  coord_cartesian(ylim = c(-0.5, 1)) +
  theme_minimal_grid() +
  labs(title = "11 cross sections", x='position', y='RI') 
# dev.off()


# ordering 65 axon ------------------------------------------------------------------------------------------------

yz_cs_glo <- list()
for (j in 2:5) {
  yz_cs_glo[[j]] <- as.matrix(df_cable_ax[((N_ind_sp+1)*(seq(1,65)-1)+j), c(1,2)])
}

xy_com_m_tp_65 <- xy_com_m_tp
rt_ind_glo <- matrix(NA, ncol = 5, nrow = 65)

# in LO
tp <- xy_com_m_tp_65 #theta phi in LO
swapMax <- (nrow(tp)-1)*(nrow(tp)-2)/2
oj <- order(apply(tp, 1, function(x) arcLength(x, tp[k,])) )
o_cs <- order(apply(tp, 1, function(x) arcLength(x, tp[k,]))[oj]) 
rt_ind_glo[,1] <- 1 - 2 * bubble_sort_swapCount(o_cs) / swapMax

# in axon
for (j in 2:5) {
  yz <- yz_cs_glo[[j]] #y z in cable
  ii <- which(!is.na(yz[,1]))
  yz <- yz[ii,]
  tp <- xy_com_m_tp_65[ii,] #theta phi in LO
  swapMax <- (nrow(tp)-1)*(nrow(tp)-2)/2
  for (k in 1:nrow(tp)) {
    oj <- order(apply(tp, 1, function(x) arcLength(x, tp[k,])) )
    o_cs <- order(rowSums(sweep(yz, 2, yz[k,], '-')^2)[oj]) 
    rt_ind_glo[ii[k],j] <- 1 - 2 * bubble_sort_swapCount(o_cs) / swapMax
  }
}


# plot
col_4_65 <- rep('grey', 65)
col_4_65 <- adjustcolor(col_4_65, alpha.f = 0.6)
col_4_65[c(i1,i2,i3,i4)] <- c(rep(col4[1],5),rep(col4[2],5),rep(col4[3],5),rep(col4[4],5))

df <- as.data.frame(cbind(seq(1,65), rt_ind_glo) )
colnames(df) <- c('neu', seq(1,5))
dfm <- melt(df, id = 'neu')
dfm$neu <- factor(dfm$neu)
dfm <- dfm[!is.na(dfm$value),]
# tmp <- dfm[66:715,] 
# tmp <- tmp[!(tmp$value==1),]
# dfm <- rbind(dfm[1:65,], tmp)

dev.new()
# pdf(file = "4 crossSection index.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  geom_jitter(aes(colour=neu),size = 4,height = 0, width = 0.1) +
  scale_colour_manual(values = col_4_65, guide=FALSE) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', size = 5, position = position_nudge(x = 0.2, y = 0)) +
  coord_cartesian(ylim = c(-0.5, 1)) +
  theme_minimal_grid() +
  labs(title = "11 cross sections", x='position', y='RI') 
dev.off()

#  compile all into a violin plot ---------------------------------------------------------------------------------
# Figure 5Supp

ks.test(reti_stats_index[,2], N_swap_man_index[,2])
ks.test(df[,3], df[,4])

# swapMax <- 64*(64-1)/2
# df <- data.frame(cbind(seq(1,65), 
#                        reti_stats_index[,c(2,4)], #cable/syn median
#                        N_swap_man_index,
#                        1 - 2*N_swap_pre_ls[[4]][,1]/swapMax,
#                        1 - 2*N_swap_pre_3_ls[[4]][,1]/swapMax,
#                        1 - 2*N_swap_pre_ls[[5]][,1]/swapMax,
#                        1 - 2*N_swap_pre_3_ls[[5]][,1]/swapMax,
#                        1 - 2*N_swap_pre_ls[[6]][,1]/swapMax,
#                        1 - 2*N_swap_pre_3_ls[[6]][,1]/swapMax,
#                        1 - 2*N_swap_pre_ls[[8]][,1]/swapMax,
#                        1 - 2*N_swap_pre_3_ls[[8]][,1]/swapMax,
#                        1 - 2*N_swap_pre_ran[,1]/swapMax)
# )
# # colnames(df) <- c('neu', 'cable_avg', 'cable_med','syn_avg', 'syn_med','sec1', 'sec2',
# #                   'sec4','sec5-3', 'sec5', 'sec5-3', 'sec6', 'sec6-3',
# #                   'phi', 'theta','cylin z', 'cylin 60')
# colnames(df) <- c('neu', 'cable_med','syn_med',
#                   'phi', 'theta','cylin z', 'cylin 90',
#                   'sec4', 'sec4-3','sec5', 'sec5-3', 'sec6','sec6-3', 'sec8','sec8-3', 'sec8 random'
#                   )
# # colnames(df) <- c('neu','sec1', 'sec2')
# dfm <- melt(df, id = 'neu')
# dev.new()
# ggplot(dfm, aes(factor(variable), value) ) + 
#   geom_violin(scale = 'width') +
#   geom_jitter(height = 0, width = 0.05) +
#   ylim(-1,1) +
#   labs(title = "violin plot",
#        x='metrics', y='index') 
# 
# 
# # -- whisker plot
# dev.new()
# ggplot(dfm, aes(factor(variable), value) ) + 
#   geom_boxplot() +
#   geom_jitter(height = 0, width = 0.05) +
#   ylim(-1,1) +
#   theme_bw() +
#   labs(title = "box plot",
#        x='metrics', y='index') 


# -- median whisker plot
swapMax <- 64*(64-1)/2

# - Figure 5 supp A, median of cable and syn
df <- data.frame(cbind(seq(1,65), reti_stats_index[,c(2,4)]))
colnames(df) <- c('neu', 'cable median', 'syn median')
dfm <- melt(df, id = 'neu')
dev.new()
# pdf(file = "index median.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
ggplot(dfm, aes(factor(variable), value) ) + 
  # geom_boxplot() +
  # geom_point() +
  geom_jitter(aes(colour=variable), size = 4, height = 0, width = 0.1) +
  scale_colour_manual(values = c("chartreuse", "forestgreen"), guide=FALSE) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', position = position_nudge(x = 0.2, y = 0), size = 5 ) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal_grid() +
  labs(title = "median", x='quantity', y='RI') 

# dev.off()

# Figure 6A
# example neuron

nopen3d()
par3d('windowRect' = c(100,100,3100,1100))
rgl.viewpoint(userMatrix = rotationMatrix(30/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0)  
              , zoom = 0.3)

j <- 1
xyz <- xyzmatrix(LC6_glo_rs[[j]]$d)
points3d(matrix(c(median(xyz[,1]), median(xyz[,2]), median(xyz[,3])), ncol = 3),
         size = 50, col="#af8dc3")

conn <- connectors(LC6[[j]])
conn_pre <- conn[conn$prepost==1, c('x','y','z')]
inglo <- rowSums(sweep(conn_pre,2,c(a,b,c), '*')) + d < 0
points3d(conn_pre[inglo,], col = "lightseagreen", size =  15)
points3d(matrix(c(median(conn_pre[inglo,1]), median(conn_pre[inglo,2]), median(conn_pre[inglo,3])), ncol = 3),
         size = 50, col='forestgreen')

plot3d(LC6_glo_rs[[j]], col = 'black', lwd = 3, WithNodes = F, lit=F)

# rgl.pop(type = "light")
# rgl.light(theta = 30, phi = 30)

shade3d(glo.msh, alpha = 0.02, lit=F)
title3d('purple cable, dark green synap')

# rgl.snapshot(filename = "median.png",fmt = "png")


# - arti
df <- data.frame(cbind(seq(1,65), N_swap_man_index[,c(1,2,4,3)]))
colnames(df) <- c('neu', 'AP', 'DV old', 'DV', 'rot') #,'cylin z', 'cylin 90')
dfm <- melt(df, id = 'neu')

dev.new()
# pdf(file = "arti mapping.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  # geom_boxplot() +
  # geom_point() +
  geom_jitter(size = 4,height = 0, width = 0.1, colour='blue') +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', size = 5,position = position_nudge(x = 0.2, y = 0) ) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal_grid() +
  labs(title = "artificial mapping", x='quantity', y='RI') 

# dev.off()

# Figure 6B
# 4 scatters together
df <- data.frame(cbind(seq(1,65), reti_stats_index[,c(2,4)],  N_swap_man_index[,c(1,4)])) 
colnames(df) <- c('neu', 'cable median', 'syn median', 'AP', 'DV')
dfm <- melt(df, id = 'neu')

dev.new()
# pdf(file = "4 scatters.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  # geom_boxplot() +
  # geom_point() +
  geom_jitter(aes(colour=variable), size = 4, height = 0, width = 0.1, alpha = 0.5) +
  # scale_colour_manual(values = c(synpc_order_pal_3[1], "forestgreen", rep('black',3)), guide=FALSE) +
  scale_colour_manual(values = c("#af8dc3", "forestgreen", rep('black',3)), guide=FALSE) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', size = 5,position = position_nudge(x = 0.2, y = 0) ) +
  coord_cartesian(ylim = c(-.5, 1)) +
  theme_minimal_grid() +
  labs(title = "4 mapping", x='quantity', y='RI') 

dev.off()

# - 1D
df <- data.frame(cbind(seq(1,65), N_swap_1D_index))
colnames(df) <- c('neu', 'med x', 'med y', 'med bf', 'min10 x', 'min10 y', 'min10 bf')
dfm <- melt(df, id = 'neu')
dev.new()
ggplot(dfm, aes(factor(variable), value) ) + 
  # geom_boxplot() +
  # geom_point() +
  geom_jitter(size = 2,height = 0, width = 0.1) +
  # stat_summary(fun.min = function(z) { quantile(z,0.25) },
  #              fun.max = function(z) { quantile(z,0.75) },
  #              geom = "errorbar",
  #              colour = 'red', size = 1, width = .2 ) +
  # stat_summary(fun = median, geom = "point",  colour = 'red', size = 5 ) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal_grid() +
  labs(title = "glo VS 1D", x='quantity', y='RI') 



# - sectors
swapMax <- 64*(64-1)/2
df <- data.frame(cbind(seq(1,65), 1 - 2*N_swap_pre_ls[[4]][,1]/swapMax,
                       1 - 2*N_swap_pre_ls[[5]][,1]/swapMax,
                       1 - 2*N_swap_pre_ls[[6]][,1]/swapMax,
                       1 - 2*N_swap_pre_456[,1]/swapMax,
                       1 - 2*N_swap_pre_ran[,1]/(65*(65-1)/2)) )
colnames(df) <- c('neu', 'sec4', 'sec5','sec6', 'sec456', 'sec456 random')
dfm <- melt(df, id = 'neu')

# Figure 6H
dev.new()
# pdf(file = "sector RI.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  geom_jitter(aes(colour=variable),size = 4,height = 0, width = 0.1) +
  scale_colour_manual(values = c(dolphin_col_456, 'black', 'black'), guide=FALSE) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', position = position_nudge(x = 0.2, y = 0),size = 5 ) +
  coord_cartesian(ylim = c(-0.5, 1)) +
  theme_minimal_grid() +
  labs(title = "sector RI", x='quantity', y='RI') 

dev.off()


library(colorspace)   ## hsv colorspace manipulations

## Function for desaturating colors by specified proportion
sat <- 2
X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(dolphin_col[4:6]))
X[2,] <- 1
dolphin_col_456 <- hsv(X[1,], X[2,], X[3,])
  
# Figure 6G, 3 sectors
nopen3d()
par3d('windowRect' = c(100,100,1600,1600))
rgl.viewpoint(userMatrix = rotationMatrix(30/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0)
              , zoom = 0.3)

# points3d(LC6_pre_glo_all, col ='grey', size=1)
ii <- LC6_pre_glo_all[,1] > glu_div[4] & LC6_pre_glo_all[,1] < glu_div[5]
points3d(LC6_pre_glo_all[ii,], col = dolphin_col_456[1], size = 10)
ii <- LC6_pre_glo_all[,1] > glu_div[5] & LC6_pre_glo_all[,1] < glu_div[6]
points3d(LC6_pre_glo_all[ii,], col = dolphin_col_456[2], size = 10)
ii <- LC6_pre_glo_all[,1] > glu_div[6] & LC6_pre_glo_all[,1] < glu_div[7]
points3d(LC6_pre_glo_all[ii,], col = dolphin_col_456[3], size = 10)
# ii <- LC6_pre_glo_all[,1] > glu_div[7] & LC6_pre_glo_all[,1] < glu_div[8]
# points3d(LC6_pre_glo_all[ii,], col ='brown')

shade3d(glo.msh, alpha = 0.1,lit = F)
title3d('sectors in synapo glo')
rgl.snapshot(filename = "sector.png",fmt = "png")


# Figure 6G, eg neuron3
# sec 4 mesh
ii <- conn_LC6_gloHull[,1] > glu_div[4] & conn_LC6_gloHull[,1] < glu_div[5]
glo4.as <- ashape3d(as.matrix(conn_LC6_gloHull[ii,]), alpha = 10000)
glo4.msh <- as.mesh3d(glo4.as)


nopen3d()
par3d('windowRect' = c(100,100,2100,1100))
rgl.viewpoint(userMatrix = rotationMatrix(30/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0)
              , zoom = 0.3)

# shade3d(glo4.msh, alpha = 0.1, col = dolphin_col_456[1], lit = F)
# planes3d(1,0,0, -glu_div[4], col = dolphin_col_456[1], lit=F)
# planes3d(1,0,0, -glu_div[5], col = dolphin_col_456[1], lit=F)

j <- 5
conn <- connectors(LC6[[j]])
conn_pre <- conn[conn$prepost==0, c('x','y','z')]
inglo <- rowSums(sweep(conn_pre,2,c(a,b,c), '*')) + d < 0
conn_pre <- conn_pre[inglo,]
insec <- conn_pre[,1] > glu_div[4] & conn_pre[,1] < glu_div[5]
conn_pre <- conn_pre[insec,]
points3d(conn_pre[1,], size = 50, col = 'black', alpha =0.8)
plot3d(LC6_glo[[j]], col = 'black', lwd = 10, WithNodes = F, alpha =0.5)
# plot3d(prune_strahler(LC6_glo[[j]], 1), col = 'black', lwd = 2, WithNodes = F)
bdot <- as.matrix(conn_pre[1,])

j <- 4
conn <- connectors(LC6[[j]])
conn_pre <- conn[conn$prepost==0, c('x','y','z')]
inglo <- rowSums(sweep(conn_pre,2,c(a,b,c), '*')) + d < 0
conn_pre <- conn_pre[inglo,]
insec <- conn_pre[,1] > glu_div[4] & conn_pre[,1] < glu_div[5]
conn_pre <- conn_pre[insec,]
points3d(conn_pre, size = 40, col = 'red', alpha =0.8)
plot3d(LC6_glo[[j]], col = 'red', lwd = 10, WithNodes = F, alpha =0.5)
ii <- which.min(rowSums(sweep(as.matrix(conn_pre),2,bdot,'-')^2))
segments3d(rbind(bdot, conn_pre[ii,]), lwd = 5, col = 'purple')

# rgl.snapshot(filename = "sector eg.png",fmt = "png")

# pick target, pick top 50% LC, look at them in LO and maybe glo --------------------------------------------------

cut_perc <- 0.3

# -- #2
syn_n <- conn_target[[2]]$tofrom_glu #syn from target to LC
neu_order <- order(syn_n, decreasing = T) #order descending
cut_ind <- match(TRUE, (cumsum(syn_n[neu_order]) > sum(syn_n)*cut_perc)) # 50% cut
neu_ind <- neu_order[1:cut_ind] # top 50%
neu_ind_t2 <- neu_ind

pal_tmp <- synpc_order_pal
pal_tmp[neu_order[1:22]] <- synpc_order_pal_3[1]
pal_tmp[neu_order[23:44]] <- synpc_order_pal_3[2]
pal_tmp[neu_order[45:65]] <- synpc_order_pal_3[3]

# Figure 5D
windows(record = F, width = 8, height = 8)
plot(bd_grid, ylim = rev(range(bd_grid$Var2)), type = "n", axes = FALSE, ann = F)
for (j in 1:length(xy_poly)) {
  points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=pal_tmp[j], cex = 3, pch = 20) #pch=1 circle, 32+j ASCII
  # points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 1, pch = 20) #pch=1 circle, 32+j ASCII
}
# for (j in neu_ind) {
#   text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = syn_n[j], pos = 1, col = 'red', offset = 0.3)
# }
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
title('target #2')

# -- #5
syn_n <- conn_target[[5]]$tofrom_glu #syn from target to LC
neu_order <- order(syn_n, decreasing = T) #order descending
cut_ind <- match(TRUE, (cumsum(syn_n[neu_order]) > sum(syn_n)*cut_perc)) # 50% cut
neu_ind <- neu_order[1:cut_ind] # top 50%
neu_ind_t5 <- neu_ind

pal_tmp <- synpc_order_pal
pal_tmp[neu_order[1:22]] <- synpc_order_pal_3[1]
pal_tmp[neu_order[23:44]] <- synpc_order_pal_3[2]
pal_tmp[neu_order[45:65]] <- synpc_order_pal_3[3]

# Figure 5D
windows(record = F, width = 8, height = 8)
plot(bd_grid, ylim = rev(range(bd_grid$Var2)), type = "n", axes = FALSE, ann = F)
for (j in 1:length(xy_poly)) {
  # points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 1, pch = 20) #pch=1 circle, 32+j ASCII
  points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=pal_tmp[j], cex = 3, pch = 20) #pch=1 circle, 32+j ASCII
}
# for (j in neu_ind) {
#   text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = syn_n[j], pos = 1, col = 'red', offset = 0.3)
# }
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
title('target #5')


# -- #9
syn_n <- conn_target[[9]]$tofrom_glu #syn from target to LC
neu_order <- order(syn_n, decreasing = T) #order descending
cut_ind <- match(TRUE, (cumsum(syn_n[neu_order]) > sum(syn_n)*cut_perc)) # 50% cut
neu_ind <- neu_order[1:cut_ind] # top 50%
neu_ind_t9 <- neu_ind

pal_tmp <- synpc_order_pal
pal_tmp[neu_order[1:22]] <- synpc_order_pal_3[1]
pal_tmp[neu_order[23:44]] <- synpc_order_pal_3[2]
pal_tmp[neu_order[45:65]] <- synpc_order_pal_3[3]

# Figure 5D
windows(record = F, width = 8, height = 8)
plot(bd_grid, ylim = rev(range(bd_grid$Var2)), type = "n", axes = FALSE, ann = F)
for (j in 1:length(xy_poly)) {
  # points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 1, pch = 20) #pch=1 circle, 32+j ASCII
  points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=pal_tmp[j], cex = 3, pch = 20) #pch=1 circle, 32+j ASCII
}
# for (j in neu_ind) {
#   text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = syn_n[j], pos = 1, col = 'red', offset = 0.3)
# }
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
title('target #9')

# - syn localization ? --> not much
nopen3d()
neu_ind <- neu_ind_t9
for (j in neu_ind) {
  # points3d(matrix(colMeans(LC6_pre_glo[[j]]),ncol = 3), col ='blue', size =10)
  # points3d(LC6_pre_glo[[j]], col ='blue', size =10)
  # text3d(matrix(colMeans(xyz),ncol = 3), texts = j)
  conn <- catmaid_get_connectors_between(pre_skids = neu[[j]]$skid,
                                 post_skids = neu_target[[9]]$skid)
  points3d(conn[,c('connector_x','connector_y', 'connector_z')], col = 'blue', size = 10)
}
plot3d(LC6[neu_ind_t9], col ='blue', lwd = 0.2)
plot3d(neu_target[[9]], col ='blue', lwd = 1)

neu_ind <- neu_ind_t5
for (j in neu_ind) {
  # points3d(matrix(colMeans(LC6_pre_glo[[j]]),ncol = 3), col ='cyan', size =10)
  # points3d(LC6_pre_glo[[j]], col ='red', size =10)
  conn <- catmaid_get_connectors_between(pre_skids = neu[[j]]$skid,
                                         post_skids = neu_target[[5]]$skid)
  points3d(conn[,c('connector_x','connector_y', 'connector_z')], col = 'red', size = 10)
}
plot3d(LC6[neu_ind_t5], col ='red', lwd = 0.2)
plot3d(neu_target[[5]], col ='red', lwd = 1)

neu_ind <- neu_ind_t2
for (j in neu_ind) {
  # points3d(matrix(colMeans(LC6_pre_glo[[j]]),ncol = 3), col ='pink', size =10)
  # points3d(LC6_pre_glo[[j]], col ='pink', size =10)
  conn <- catmaid_get_connectors_between(pre_skids = neu[[j]]$skid,
                                         post_skids = neu_target[[2]]$skid)
  points3d(conn[,c('connector_x','connector_y', 'connector_z')], col = 'cyan', size = 5)
}

# histo of syn, bi VS ipsi ----------------------------------------------------------------------------------------
# x <- seq(5,25,by = 5) - 2.5
x <- seq(0,23,by = 1) 
syn_hist <- matrix(ncol = 9, nrow = length(x))
for (j in 1:9) {
  syn_n <- conn_target[[j]]$tofrom_glu
  # y <- hist(syn_n, breaks = seq(-1,23))$density
  # y <- hist(syn_n, breaks = seq(-1,23))$count
  y <- hist(syn_n, breaks = seq(-1,23,by = 1))$count
  syn_hist[,j] <- y
}

syn_hist <- as.data.frame(syn_hist)
df <- cbind(x, syn_hist)
colnames(df)[1] <- 'num'

dfm <- melt(df, id = 'num')
dev.new()
ggplot(dfm, aes(x=num, y=value, colour=variable) ) + 
  geom_line(lwd = 3) +
  scale_colour_manual(values = brewer.pal(9, 'RdBu')) +
                      # breaks = c("1", "2", "3","4"),
                      # labels = c("top", "back", "bottom","front") ) 
  theme_bw() +
  labs(title = "syn histogram 4 bi + 5 ipsi ",
       x='num of syn', y='density') 


# - cumsum
x <- seq(1,65) 
syn_cumsum <- matrix(ncol = 9, nrow = length(x))
for (j in 1:9) {
  syn_n <- conn_target[[j]]$tofrom_glu
  syn_cumsum[,j] <- cumsum(sort(syn_n,decreasing = T))
  syn_cumsum[,j] <- syn_cumsum[,j]/max(syn_cumsum[,j])
}
syn_cumsum <- as.data.frame(syn_cumsum[,3:9])
df <- cbind(x, syn_cumsum)
colnames(df)[1] <- 'num'

dfm <- melt(df, id = 'num')

getPalette <- colorRampPalette(brewer.pal(9, "RdYlBu"))
dolphin_col <- getPalette(11)
col9 <- dolphin_col[-c(5,6)]

dev.new()
ggplot(dfm, aes(x=num, y=value, colour=variable) ) + 
  geom_line(lwd = 2) +
  scale_colour_manual(values = col9, labels = names(sapply(neu_target, function(x) x$skid)))  +
  theme_bw() +
  guides(colour = guide_legend(title="target neurons") ) +
  labs(title = "syn histogram 4 bi + 5 ipsi ",
       x='num of syn', y='density') 

# 1d proj ---------------------------------------------------------------------------------------------------------

j <- 1

# EM
mm <- simdata_df[[j]]
# mm <- mm[mm$x >= -9 & mm$x <= 111, ]
mm <- mm[mm$x >= -9 & mm$x <= 80, ]
mm <- mm[mm$y >= 63 & mm$y <= 108, ]
pj1_sim <- tapply(mm$z, mm$x, sum) #along elev or azim
pj1_sim <- pj1_sim/max(pj1_sim)

# ephys
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
fdr[[j]] <- p.adjust(pval_tmp, method = 'BH')

exp_df <- data.frame(xygrid2, as.vector(t(exp_raw_mean[[j]])))
colnames(exp_df) <- c("x","y","z")
exp_loess <- loess (z ~ x * y, exp_df, degree = 2, span = 0.3)
loess_fit <- predict (exp_loess, xygrid4, se = T)
exp_interp <- cbind(xygrid4, melt(loess_fit$fit)$value)
colnames(exp_interp) <- c('x','y','z')
exp_interp <- exp_interp[!is.na(exp_interp[,'z']), ]

mm <- exp_interp
# mm <- mm[mm$x >= -9 & mm$x <= 111, ]
mm <- mm[mm$x >= -9 & mm$x <= 80, ]
mm <- mm[mm$y >= 63 & mm$y <= 108, ]
pj1_ephys <- tapply(mm$z, mm$x, sum)
pj1_ephys <- pj1_ephys/max(pj1_ephys)

dev.new()
plot(pj1_sim, col = 'pink')
points(pj1_ephys, col = 'blue')





# LC6LC6 RI -------------------------------------------------------------------------------------------------------
N_swap_synLC6LC6 <- c()
for (j in 1:65) {
  oj <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[j,])), decreasing = F)
  Nsyn_tmp <- c()
  for (k in 1:65) {
    if (k == j) {
      Nsyn_tmp <- c(Nsyn_tmp, 0)
    } else {
      i_bool <- (dist_Nconn$from == j & dist_Nconn$to == k)  | (dist_Nconn$from == k & dist_Nconn$to == j) 
      Nsyn_tmp <- c(Nsyn_tmp, dist_Nconn[i_bool, "Nconn_glu"])
    }
  }
  if (length(Nsyn_tmp) != 65 ) {
    print(j)
  }
  oNsyn <- order(Nsyn_tmp[oj], decreasing = F) 
  N_swap_synLC6LC6 <- c(N_swap_synLC6LC6, bubble_sort_swapCount(oNsyn))
}


swapMax <- 64*(64-1)/2
1 - 2*N_swap_synLC6LC6/swapMax

df <- data.frame(cbind(seq(1,65), 1 - 2*N_swap_synLC6LC6/swapMax) )
colnames(df) <- c('neu', 'Nsyn')
dfm <- melt(df, id = 'neu')

dev.new()
pdf(file = "LC6LC6 Nsyn RI.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  geom_jitter(colour="blue",size = 4,height = 0, width = 0.1) +
  # scale_colour_manual(values = c(dolphin_col_456, 'blue', 'blue'), guide=FALSE) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', position = position_nudge(x = 0.2, y = 0),size = 5 ) +
  coord_cartesian(ylim = c(-1, 1)) +
  theme_minimal_grid() +
  labs(title = "LC6LC6 Nsyn RI", x='quantity', y='RI') 

dev.off()


# target Nsyn vs 3 groups LC --------------------------------------------------------------------------------------

Nconn_target_3gp <- matrix(ncol = 9, nrow = 3)
for (j in 1:9) {
  for (k in 1:3) {
    # Nconn_target_3gp[k,j] <- sum(conn_target[[j]]$Nconn_glu[cable_glo_x_o$gp == k])
    Nconn_target_3gp[k,j] <- sum(conn_target[[j]]$fromto_glu[cable_glo_x_o$gp == k])
  }
}

df <- data.frame(cbind(1:3, Nconn_target_3gp))
colnames(df) <- c('gp', seq(1,9))
dfm <- melt(df, id='gp')
# dfm$variable <- factor(dfm$variable)
dev.new()
ggplot(dfm, aes(x=gp, y=value, colour = variable)) +
  geom_point(size = 10)


# RI cartoon and stats --------------------------------------------------------------------------------------------------

# func, calculate indiv RI 
RI <- function(pt0, pt1){
  pt0 <- as.matrix(pt0)
  pt1 <- as.matrix(pt1)
  
  if (nrow(pt0) != nrow(pt1)) {
    stop("pt num don't match")
  }
  
  nr <- nrow(pt0)
  swapMean <- (nr-1)*(nr-2)/4
  
  N_swap <- c()
  for (j in 1:nr) {
    o_pt0 <- sweep(pt0, 2, pt0[j,],'-')^2 %>% rowSums() %>% order()
    N <- sweep(pt1, 2, pt1[j,],'-')^2 %>% rowSums() %>% .[o_pt0] %>% order() %>% bubble_sort_swapCount()
    N_swap <- c(N_swap, N)
  }
  return(1 - N_swap/swapMean)
}


# func, calculate avg RI for 2D points config
RI_2D_avg <- function(pt0, pt1){
  nr <- nrow(pt0)
  swapMean <- (nr-1)*(nr-2)/4
  N_swap <- c()
  for (j in 1:nr) {
    o_pt0 <- sweep(pt0, 2, pt0[j,],'-')^2 %>% rowSums() %>% order()
    N <- sweep(pt1, 2, pt1[j,],'-')^2 %>% rowSums() %>% .[o_pt0] %>% order() %>% bubble_sort_swapCount()
    N_swap <- c(N_swap, N)
  }
  return(1 - sum(N_swap)/nr/swapMean)
}

# func, calculate avg RI for 2D points config proj on x and y axis
RI_2D_proj_avg <- function(pt0, pt1){
  x0 <- pt0[,1]
  y0 <- pt0[,2]
  
  x1 <- pt1[,1]
  y1 <- pt1[,2]
  
  nr <- nrow(pt0)
  swapMean <- (nr-1)*(nr-2)/4
  N_swap_x <- c()
  N_swap_y <- c()
  for (j in 1:nr) {
    o_x0 <- (x0 - x0[j])^2 %>% order()
    N <- (x1 - x1[j])^2 %>% .[o_x0] %>% order() %>% bubble_sort_swapCount()
    N_swap_x <- c(N_swap_x, N)
    
    o_y0 <- (y0 - y0[j])^2 %>% order()
    N <- (y1 - y1[j])^2 %>% .[o_y0] %>% order() %>% bubble_sort_swapCount()
    N_swap_y <- c(N_swap_y, N)
  }
  return(c(1 - sum(N_swap_x)/nr/swapMean, 1 - sum(N_swap_y)/nr/swapMean) )
}

# -- pts on grid
N <- 8
x <- seq(-N/2+0.5, N/2-0.5, length.out = 8)
y <- seq(-N/2+0.5, N/2-0.5, length.out = 8)
pt0 <- expand.grid(x, y) + cbind(runif(64,min = -1,max = 1)*0.1, runif(64,min = -1,max = 1)*0.1)
pt0 <- as.matrix(pt0)


# -- pts on grid
Nr <- 6
Nc <- 11
y <- seq(-Nr/2+0.5, Nr/2-0.5, length.out = Nr)
x <- seq(-Nc/2+0.5, Nc/2-0.5, length.out = Nc)
pt0 <- expand.grid(x, y) + cbind(runif(Nr*Nc,min = -1,max = 1)*0.1, runif(Nr*Nc,min = -1,max = 1)*0.1)
pt0 <- as.matrix(pt0)


# - mapping

# swap points
pt1 <- pt0
pt1[c(43,59),] <- pt1[c(59,43),]
# pt1[c(1,64),] <- pt1[c(64,1),]
# pos01 <- match(pt0[,1], pt1[,1])



# random
pt1 <- pt0[sample(64),]

# swap rows
pt1 <- pt0
pt1[(1:8)+2*8,] <- pt0[(1:8)+8*7,]
pt1[(1:8)+7*8,] <- pt0[(1:8)+2*8,]

# 1D proj
pt1 <- pt0
pt1[,2] <- pt1[,2] / 2
pt1 <- pt1 + cbind(runif(Nr*Nc,min = -1,max = 1)*0.1, runif(Nr*Nc,min = -1,max = 1)*0.1)

# cos
pt1 <- pt0
pt1[,2] <- cos(pt1[,2]/2)

# cylinder
pt1 <- cbind(cos(pt0[,1]/8*2*pi), sin(pt0[,1]/8*2*pi), pt0[,2]/7*2*pi)

nopen3d()
points3d(pt1, col = col_64Blue, size = 10)
axes3d(c('x','y','z')); title3d('','','x','y','z')
# rgl.bbox(color="gray90", alpha = 1, 
#          # emission="#333377", specular="#333377", shininess=5,
#          xlen = 0, ylen = 0, zlen = 0)

# Gaussian
z <- 10 * exp( - rowSums(pt0^2) / 2 / 4^2)
pt1 <- cbind(pt0,z)
colnames(pt1) <- c('x','y','z')

nopen3d()
points3d(pt1, col = col_64Blue, size = 10)
rgl.bbox(color="gray90", alpha = 1, 
         # emission="#333377", specular="#333377", shininess=5,
         xlen = 0, ylen = 0, zlen = 0)


# jittering
pt1 <- expand.grid(x, y) + cbind(runif(64,min = -1,max = 1)*0.5, runif(64,min = -1,max = 1)*0.5)
pt1 <- as.matrix(pt1)

# triangle
# pt1 <- expand.grid(x, y) + cbind(runif(64,min = -1,max = 1)*0.5, runif(64,min = -1,max = 1)*0.5)
pt1 <- pt0
for (j in 1:8) {
  pt1[seq(1,8)+(j-1)*8,1] <- pt1[seq(1,8)+(j-1)*8,1]*0.2*j
}

# local scramble
pt1 <- pt0
ii <- pt0[,1] > -2 & pt0[,1] < 2 & pt0[,2] > -2 & pt0[,2] < 2
pt1[ii,] <- pt0[ii,][sample(sum(ii)),] + cbind(runif(sum(ii),min = -1,max = 1)*0.5, runif(sum(ii),min = -1,max = 1)*0.5)


# neg val ?
pt1 <- pt0
m <- 11
for (j in seq(1,60,by = 3)) {
  pt1[c(j,j+m)%%64+1,] <- pt1[c(j+m,j)%%64+1,]
}
m <- 7
for (j in seq(1,60,by = 3)) {
  pt1[c(j,j+m)%%64+1,] <- pt1[c(j+m,j)%%64+1,]
}


pt1 <- 1/(sweep(pt0,2,c(1,2),'+'))

# RI
RI_2D_avg(pt0, pt1)

RI_2D_proj_avg(pt0, pt1)


# -- plot
getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
col_64Blue <- rev(rev(getPalette(40))[1:(Nr*Nc/2)])
# col_64Blue <- rev(rev(getPalette(40))[1:32])
getPalette <- colorRampPalette(brewer.pal(9, "Reds"))
col_64Blue <- c(col_64Blue, rev(getPalette(40))[1:(Nr*Nc/2)])


df0 <- pt0 %>% cbind(.,seq(1,Nr*Nc)) %>% data.frame()
colnames(df0) <- c('x','y', 'no')

dev.new()
ggplot() +
  geom_point(data = df0, aes(x,y, colour = factor(no)), size =10,  shape = 16) + 
  scale_color_manual(values = col_64Blue) +
  # labs(title = mat_names[j]) +
  theme_void() +
  coord_fixed(ratio = 1) +
  theme(legend.position = "none")



df1 <- pt1 %>% cbind(.,seq(1,Nr*Nc)) %>% data.frame()
colnames(df1) <- c('x','y', 'no')

dev.new()
ggplot() +
  geom_point(data = df1, aes(x,y, colour = factor(no)), size =10,  shape = 16) + 
  scale_color_manual(values = col_64Blue) +
  # labs(title = mat_names[j]) +
  theme_bw() +
  coord_fixed(ratio = 1) +
  theme(legend.position = "none")



# -- rot sweep 

# -- RI for all rot ang
RI_2D_rot <- matrix(ncol = 18, nrow = Nr*Nc)
RI_1D_x_rot <- matrix(ncol = 18, nrow = Nr*Nc)
RI_1D_y_rot <- matrix(ncol = 18, nrow = Nr*Nc)
RI_12_rot <- matrix(ncol = 18, nrow = Nr*Nc)
nr <- nrow(pt0)
swapMean <- (nr-1)*(nr-2)/4
for (k in 1:18) {
  ang <- k*10/180*pi
  M_rot <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol = 2, byrow = T) # cw
  pt1 <- pt0 %*% M_rot
  
  # - 2D
  N_swap <- c()
  for (j in 1:nr) {
    o_pt0 <- sweep(pt0, 2, pt0[j,],'-')^2 %>% rowSums() %>% order()
    N <- sweep(pt1, 2, pt1[j,],'-')^2 %>% rowSums() %>% .[o_pt0] %>% order() %>% bubble_sort_swapCount()
    N_swap <- c(N_swap, N)
  }
  
  #  - 1D
  x0 <- pt0[,1]
  y0 <- pt0[,2]
  
  x1 <- pt1[,1]
  y1 <- pt1[,2]
  
  N_swap_x <- c()
  N_swap_y <- c()
  N_swap_12 <- c()
  for (j in 1:nr) {
    o_x0 <- (x0 - x0[j])^2 %>% order()
    N <- (x1 - x1[j])^2 %>% .[o_x0] %>% order() %>% bubble_sort_swapCount()
    N_swap_x <- c(N_swap_x, N)
    
    o_y0 <- (y0 - y0[j])^2 %>% order()
    N <- (y1 - y1[j])^2 %>% .[o_y0] %>% order() %>% bubble_sort_swapCount()
    N_swap_y <- c(N_swap_y, N)
    
    o_pt0 <- sweep(pt0, 2, pt0[j,],'-')^2 %>% rowSums() %>% order()
    N <- (x1 - x1[j])^2 %>% .[o_pt0] %>% order() %>% bubble_sort_swapCount()
    N_swap_12 <- c(N_swap_12, N)
  }
  
  RI_2D_rot[,k] <- 1 - N_swap/swapMean
  RI_1D_x_rot[,k] <- 1 - N_swap_x/swapMean
  RI_1D_y_rot[,k] <- 1 - N_swap_y/swapMean
  RI_12_rot[,k] <- 1 - N_swap_12/swapMean
}


# -- 1D vs glo sim

ang <- 3*10/180*pi
M_rot_30 <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol = 2, byrow = T) # cw
pt2_sj <- pt0 + cbind(runif(Nr*Nc,min = -1,max = 1)*0, runif(Nr*Nc,min = -1,max = 1)*3)
pt2_sj <- pt2_sj %*% M_rot_30 

RI_12_rot_jitter <- matrix(ncol = 18, nrow = Nr*Nc)

nr <- nrow(pt0)
swapMean <- (nr-1)*(nr-2)/4
for (k in 1:18) {
  ang <- k*10/180*pi
  M_rot <- matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), ncol = 2, byrow = T) # cw
  pt1 <- pt0 %*% M_rot
  
  x1 <- pt1[,1]
  y1 <- pt1[,2]
  
  N_swap_12 <- c()
  for (j in 1:nr) {
    o_pt0 <- sweep(pt2_sj, 2, pt2_sj[j,],'-')^2 %>% rowSums() %>% order()
    N <- (x1 - x1[j])^2 %>% .[o_pt0] %>% order() %>% bubble_sort_swapCount()
    N_swap_12 <- c(N_swap_12, N)
  }
  RI_12_rot_jitter[,k] <- 1 - N_swap_12/swapMean
}


# PLOT
df <- data.frame(cbind(seq(1,Nr*Nc), RI_12_rot))
df <- data.frame(cbind(seq(1,Nr*Nc), RI_2D_rot))
df <- data.frame(cbind(seq(1,Nr*Nc), RI_1D_x_rot))
# df <- data.frame(cbind(seq(1,Nr*Nc), RI_1D_y_rot)) 
df <- data.frame(cbind(seq(1,Nr*Nc), RI_12_rot_jitter))

colnames(df) <- c('neu', seq(0,170,by = 10))
dfm <- melt(df, id = 'neu')

dev.new()
# pdf(file = "index rotation.pdf", width = 10, height = 4,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  # geom_boxplot() +
  # geom_point() +
  # geom_jitter(colour = 'black', size = 4, height = 0, width = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', size = 5,position = position_nudge(x = 0, y = 0) ) +
  coord_cartesian(ylim = c(-.5, 1)) +
  theme_minimal_grid() +
  labs(title = "grid rotation 12 jitter_x wide", x='Angle', y='RI') 

# dev.off()



# -- stats
# as.matrix(rperm(3,6))

# - 1d single point switch
Np <- 64




swap_stats <- matrix(ncol = 2, nrow = 65)
# - variance in num of swaps
for (j in 1:8) { #exact
  pp <- gtools::permutations(j,j)
  nn <- apply(pp, 1, function(x) bubble_sort_swapCount(x))
  swap_stats[j,] <- c(mean(nn), sd(nn))
}

for (j in 9:65) { #100k sample
  pp <- as.matrix(permutations::rperm(1e5, j))
  nn <- apply(pp, 1, function(x) bubble_sort_swapCount(x))
  swap_stats[j,] <- c(mean(nn), sd(nn))
}
swap_stats[1,] <- c(0,0)

# plot mean and sd swaps
df <- as.data.frame(cbind(seq(1,65),swap_stats))
colnames(df) <- c('N', 'mean','sd')

dev.new()
ggplot(df) + 
  geom_ribbon(aes(x = N, ymin = mean-sd, ymax = mean+sd), fill = "grey") +
  geom_line(aes(x=N, y=mean), colour = 'black', lwd = 2) +
  xlab('num of elements') +
  ylab('num of swaps') +
  labs(title = paste("swaps number and uncertainty", sep = ',')) 

# plot, sd vs index
df <- data.frame(N = seq(1,65), sd_ind = 2 * swap_stats[,2] / (65*64/2))
dev.new()
ggplot(df)+
  geom_line(aes(x=N, y=sd_ind), colour = 'black', lwd = 2) +
  xlab('num of mis-ordered elements') +
  ylab('index') +
  labs(title = paste("effect on index due to order ambiguity", sep = ',')) 

# -- retinotopy index vs num swaps, 1 - 2*(swaps / max swaps)
smax <- 65*(65-1)/2
x <- seq(1, smax)
y <- 1 - 2*x/smax

dev.new()
plot(x,y)

# T-bar target dist vs num synp -----------------------------------------------------------------------------------

neu_target_rs <- resample(neu_target, stepsize = 400)


n_d_ls <- list() # number of syn and avg dist to all LC6 pre-syn
p_r <- matrix(ncol = 4, nrow = length(neu_target_rs))
for (t in 1:length(neu_target_rs)) {
  n_d <- matrix(ncol = 3, nrow = 65)
  txyz <- xyzmatrix(neu_target_rs[[t]]$d)
  for (n in 1:length(LC6_pre_glo)) {
    lxyz <- as.matrix(LC6_pre_glo[[n]])
    # lxyz <- as.matrix(LC6_post_glo[[n]])
    dd_c <- 0
    for (m in 1:dim(lxyz)[1]) {
      dd <- min(sqrt(rowSums(sweep(txyz, 2, lxyz[m,], '-')^2)))
      dd_c <- c(dd_c, dd)
    }
    qq <- quantile(dd_c, probs = c(0.25, 0.5, 1))
    # dd_q <- dd_c[dd_c < qq[3]]
    # dd_avg <- median(dd_q)
    n_d[n,] <- c(conn_target[[t]]$tofrom_glu[n],
                 median(dd_c[dd_c < qq[1]]),
                 median(dd_c[dd_c < qq[3]])  )
  }
  n_d <- as.data.frame(n_d)
  colnames(n_d) <- c('num', 'med25', 'med100')
  linReg <- lm(med25 ~ num, data = n_d)
  n_d$fitted25 <- linReg$fitted.values
  pval <- summary(linReg)$coefficients['num', 'Pr(>|t|)'] #p-value
  rsq <- summary(linReg)$r.squared # R^2
  p_r[t,1:2] <- c(pval, rsq)
  
  linReg <- lm(med100 ~ num, data = n_d)
  n_d$fitted100 <- linReg$fitted.values
  n_d_ls[[t]] <- n_d
  pval <- summary(linReg)$coefficients['num', 'Pr(>|t|)'] #p-value
  rsq <- summary(linReg)$r.squared # R^2
  p_r[t,3:4] <- c(pval, rsq)
  
  # pval <- summary(linReg)$coefficients['num', 'Pr(>|t|)'] #p-value
  # rsq <- summary(linReg)$r.squared # R^2
  # p_r[t,] <- c(pval, rsq)
}


# plot
# mycols <- adjustcolor(palette(), alpha.f = 0.3)
# opal <- palette(mycols)

dev.new()
par(mfrow=c(2,5))
for (j in 1:9) {
  # plot(n_d_ls[[j]][,1:2], pch=16, cex=2, col=4)
  plot(n_d_ls[[j]][,1:2], pch=16, cex=2, col=4,
       xlim = c(0,25), ylim = c(0,4000), main = j,
       cex.axis=2, cex.lab=2, cex.main=2)
  lines(n_d_ls[[j]][,c(1,4)], lwd=2, col='black')
  text(15, 2000, labels = bquote(paste("p-value = ", .(signif(p_r[j,1],2)))),
       adj = 0, srt = 0, cex = 2)
  text(15, 2500, labels = bquote(paste(R^2, " = ", .(round(p_r[j,2],2)))),
       adj = 0, srt = 0, cex = 2)
  
  points(n_d_ls[[j]][,c(1,3)], pch=15, cex=2, col=4,
         xlim = c(0,25), ylim = c(0,4000), main = j,
         cex.axis=2, cex.lab=2, cex.main=2)
  lines(n_d_ls[[j]][,c(1,5)], lwd=2, col='black')
  text(15, 1000, labels = bquote(paste("p-value = ", .(signif(p_r[j,3],2)))),
       adj = 0, srt = 0, cex = 2)
  text(15, 500, labels = bquote(paste(R^2, " = ", .(round(p_r[j,4],2)))),
       adj = 0, srt = 0, cex = 2)
}
plot(c(0,0))
title("median vs median of 25 percentile")


dev.new()
par(mfrow=c(2,5))
for (j in 1:9) {
  # plot(n_d_ls[[j]][,1:2], pch=16, cex=2, col=4)
  plot(n_d_ls[[j]][,1:2], pch=16, cex=2, col=4,
       xlim = c(0,25), ylim = c(0,2000), main = j,
       cex.axis=2, cex.lab=2, cex.main=2)
  lines(n_d_ls[[j]][,c(1,3)], lwd=2, col='black')
  text(15, 1500, labels = bquote(paste("p-value = ", .(signif(p_r[j,1],2)))),
       adj = 0, srt = 0, cex = 2)
  text(15, 1400, labels = bquote(paste(R^2, " = ", .(round(p_r[j,2],2)))),
       adj = 0, srt = 0, cex = 2)
}

# plot(n_d_ls[[1]][,1:2], pch=16, cex=2, col=4)
# lines(n_d_ls[[1]][,c(1,3)], lwd=2, col='black')
# 
# plot(n_d_ls[[2]])
# plot(n_d_ls[[3]])
# plot(n_d_ls[[4]])
# plot(n_d_ls[[5]])
# plot(n_d_ls[[6]])
# plot(n_d_ls[[7]])
# plot(n_d_ls[[8]])
# plot(n_d_ls[[9]])


# map delta RF to glo ---------------------------------------------------------------------------------------------

ind_bi_ipsi <- c(48,1,4,39,44,14,19,8)
LC6_bi_ipsi <- LC6_glo[ind_bi_ipsi]

syn <- conn_LC6_tar[[3]][conn_LC6_tar[[3]]$ID %in% ind_bi_ipsi, ]

nopen3d()
points3d(xyzmatrix(syn), size = 10, col = 'cyan')
shade3d(glo.msh, col='grey',alpha = 0.2)
plot3d(LC6_bi_ipsi, col ='blue',WithNodes = F)


conn_btw <- catmaid_get_connectors_between(pre_skids = neu[[1]]$skid, post_skids = neu_target[[3]]$skid)
conn_btw[, c('connector_x','connector_y','connector_z')]

conn <- connectors(LC6[[1]])
conn_pre <- conn[conn$prepost==0, c('x','y','z')]
conn_pre_glo <- conn_pre[pointsinside(conn_pre, surf=glo.msh, rval='distance') > - 1000, ] 


nopen3d()
points3d(conn_pre_glo, size =20, col = 'green', alpha = 0.5)
# points3d(conn_pre, size =13, col = 'red', alpha = 0.4)
points3d()


# --distance betw pre-synap connectors, is clustering ?
# LC6_pre_glo
dd <- c()
for (j in ind_bi_ipsi) {
  ddj <- c()
  for (k in 1:nrow(LC6_pre_glo[[j]])) {
    ddk <- c()
    for (m in ind_bi_ipsi) {
      if (m != j) {
        ddm <- sweep(as.matrix(LC6_pre_glo[[m]]), 2, as.matrix(LC6_pre_glo[[j]][k,]))^2 %>% rowSums() %>% sqrt() %>% min()
        ddk <- c(ddk, ddm)
      }
    }
    ddj <- c(ddj, mean(ddk))
  }
  dd <- c(dd, mean(ddj))
}

# all neurons
dd65 <- c()
for (j in 1:65) {
  ddj <- c()
  for (k in 1:nrow(LC6_pre_glo[[j]])) {
    ddk <- c()
    for (m in 1:65) {
      if (m != j) {
        ddm <- sweep(as.matrix(LC6_pre_glo[[m]]), 2, as.matrix(LC6_pre_glo[[j]][k,]))^2 %>% rowSums() %>% sqrt() %>% min()
        ddk <- c(ddk, ddm)
      }
    }
    ddj <- c(ddj, mean(ddk))
  }
  dd65 <- c(dd65, mean(ddj))
}

ks.test(dd,dd65)

# random 8
ddran_sample <- c()
for (r in 1:100) {
  ind_ran <- sample(x = 65, size = 8)
  ddran <- c()
  for (j in ind_ran) {
    ddj <- c()
    for (k in 1:nrow(LC6_pre_glo[[j]])) {
      ddk <- c()
      for (m in ind_ran) {
        if (m != j) {
          ddm <- sweep(as.matrix(LC6_pre_glo[[m]]), 2, as.matrix(LC6_pre_glo[[j]][k,]))^2 %>% rowSums() %>% sqrt() %>% min()
          ddk <- c(ddk, ddm)
        }
      }
      ddj <- c(ddj, mean(ddk))
    }
    ddran <- c(ddran, mean(ddj))
  }
  ddran_sample <- c(ddran_sample, mean(ddran))
}


# 8 nb
ddran_nb <- c()
for (r in 1:100) {
  ind_samp <- sample(x = 65, size = 1)
  ind_nb <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[ind_samp,])), decreasing = F)[1:8]
  ddran <- c()
  for (j in ind_nb) {
    ddj <- c()
    for (k in 1:nrow(LC6_pre_glo[[j]])) {
      ddk <- c()
      for (m in ind_nb) {
        if (m != j) {
          ddm <- sweep(as.matrix(LC6_pre_glo[[m]]), 2, as.matrix(LC6_pre_glo[[j]][k,]))^2 %>% rowSums() %>% sqrt() %>% min()
          ddk <- c(ddk, ddm)
        }
      }
      ddj <- c(ddj, mean(ddk))
    }
    ddran <- c(ddran, mean(ddj))
  }
  ddran_nb <- c(ddran_nb, mean(ddran))
}



# nb
ind_nb <- c(30,36,53,32,21,5,35,60)
ind_nb <- c(30,36,53,21,61,58,33,7)
dd <- c()
for (j in ind_nb) {
  ddj <- c()
  for (k in 1:nrow(LC6_pre_glo[[j]])) {
    ddk <- c()
    for (m in ind_nb) {
      if (m != j) {
        ddm <- sweep(as.matrix(LC6_pre_glo[[m]]), 2, as.matrix(LC6_pre_glo[[j]][k,]))^2 %>% rowSums() %>% sqrt() %>% min()
        ddk <- c(ddk, ddm)
      }
    }
    ddj <- c(ddj, mean(ddk))
  }
  dd <- c(dd, mean(ddj))
}

mean(dd)


# PLOT
range(c(ddran_sample, ddran_nb))
x <- seq(2000,6500, length.out = 15)
y_nb <- hist(ddran_nb, breaks = x, plot = F)$counts
y_sample <- hist(ddran_sample, breaks = x, plot = F)$counts

dev.new()
plot(x[-1], y_nb, type = 'l', col='red')
points(x[-1], y_sample, type = 'l', col='black')

# pairwise distance pre-syn node, vs arcLength --------------------------------------------------------------------

dd_pair <- matrix(ncol = 3, nrow = 65*64/2) 
n <- 1
for (j in 1:64) {
  for (k in (j+1):65) {
    dd <- c()
    for (m in 1:nrow(LC6_pre_glo[[j]])) {
      ddm <- sweep(as.matrix(LC6_pre_glo[[k]]), 2, as.matrix(LC6_pre_glo[[j]][m,]))^2 %>%
        rowSums() %>%
        sqrt() %>%
        min()
      dd <- c(dd, ddm)
    }
    for (m in 1:nrow(LC6_pre_glo[[k]])) {
      ddm <- sweep(as.matrix(LC6_pre_glo[[j]]), 2, as.matrix(LC6_pre_glo[[k]][m,]))^2 %>%
        rowSums() %>%
        sqrt() %>%
        min()
      dd <- c(dd, ddm)
    }
    dd_pair[n,] <- c(j,k, mean(dd))
    n <- n + 1
  }
}
    
df <- as.data.frame(dd_pair)
colnames(df) <- c('j','k', 'dd')
mat_names <- seq(1,65)

dev.new()
ggplot(df, aes(x = j, y = k)) + 
  geom_raster(aes(fill=dd)) + 
  # geom_point(aes(size=value)) +
  scale_fill_gradient(low="grey90", high="red", trans = "log", breaks = c(1,10,100,1000,10000), labels = c(1,10,100,1000,10000)) +
  guides(fill = guide_colourbar(title = "dist")) +
  # guides(fill = guide_colourbar(label = c("1","10","100","1000"))) +
  # geom_text(aes(label = dd)) +
  labs(x="j", y="k", title="dist_glo Matrix") +
  scale_x_continuous(breaks = seq(1,65),labels = mat_names, position = "top", expand = c(0,0)) +
  scale_y_reverse(breaks = seq(1,65), labels = mat_names, expand = c(0,0)) +
  # theme_bw() +
  theme(axis.text.x = element_text(face="bold", color="#993333", size=10, angle=90, hjust = 0, vjust = 0.5),
        axis.text.y = element_text(face="bold", color="#993333", size=10, angle=0),
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "black"))


# -- graph
library(igraph)
df <- as.data.frame(dd_pair)
colnames(df) <- c('from','to', 'weight')
g <- graph.data.frame(df, directed = F)
adjM <- igraph::get.adjacency(g, sparse = F, attr = 'weight')

lc <- cluster_louvain(g)



# - clustering in LO
windows(record = F, width = 10.5, height = 10.5)

bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
bd_theta <- seq(1, 180, by = 1)
xy_bd <- matrix(ncol = 2)
bd_grid <- expand.grid(bd_phi, bd_theta)
# range(bd_grid$Var1)
plot(bd_grid, xlim = c(-20,215), ylim = c(190,-50), type = "n", axes = FALSE, ann = F)
for (j in 1:length(xy_poly)) {
  xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
  if (lc$membership[j] == 1) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[1], cex = 3, pch = 16)
  }
  else if (lc$membership[j] == 2) {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[2], cex = 3, pch = 16)
  }
  else {
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="grey", cex = 2, pch = 16)
    # points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="grey", type = '.') 
  }
  text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = round(cable_glo_min1_x[j],0), pos = 1, offset = 0.3)
}
xy_bd <- xy_bd[-1,]
hpts_2 <- chull(xy_bd)
hpts_2 <- c(hpts_2, hpts_2[1])
xy_bd_chull <- xy_bd[hpts_2,] # hull edge points
polygon(xy_bd_chull)
# lines(rbind(c(-11,180), c(-2,180)), lwd = 1)
# text(-5, 180, labels = expression(paste("9",degree)), pos = 1, offset = 0.3)
# lines(rbind(c(-11,180), c(-11,171)), lwd = 1)
# text(-11, 175, labels = expression(paste("9",degree)), pos = 2, offset = 0.2)
# lines(rbind(c(0,0), c(0,180)), lwd = 1) 
# text(0, -5, labels = "front", pos = 1, offset = 0)
# lines(rbind(c(90,0), c(90,180)), lwd = 1)
# text(90, -5, labels = expression(paste("side 90",degree)), pos = 1, offset = 0)
# lines(rbind(c(-12,90), c(162,90)), lwd = 1)
# text(-17, 90, labels = "equator", pos = 1, offset = 0, srt = 90)



modularity(gc)

# -- hierarchical clustering
# hclust()

gc <- cluster_fast_greedy(g)
dev.new()
plot(gc, g)
plot_dendrogram(gc)

# dist in LO
d_LO <- dist(sph2cartZ(cbind(1,xy_com_m_tp)))
hc <- hclust(d_LO, method = "ward.D2")
hc <- hclust(d_LO)
dev.new()
plot(hc)


# make dist from avg dist
mat <- matrix(ncol = 65, nrow = 65)
ii <- 1
for (j in 1:64) {
  for (k in (j+1):65) {
    mat[k,j] <- dd_pair[ii,3]
    ii <- ii + 1
  }
}
mat <- as.dist(mat)
hc2 <- hclust(mat, method = "ward.D2")
dev.new()
plot(hc2)


# - LC6 connections
conn <- conn_LC6LC6[, c(1,2)]
conn_g <- graph.data.frame(conn, directed = T)
adjM <- get.adjacency(conn_g, sparse = F)
mat3 <- as.dist(adjM)
attributes(mat3)$Labels <- NULL
hc3 <- hclust(mat3, method = "complete")
hc4 <- hclust(mat3, method = "ward.D2")
hc4 <- hclust(mat3)

# - cp 
library(dendextend)

dend1 <- as.dendrogram(hc)
dend2 <- as.dendrogram(hc2)
dend3 <- as.dendrogram(hc3)
dend4 <- as.dendrogram(hc4)

dev.new()
dendlist(dend1, dend1) %>% 
  untangle(method = "step1side") %>% 
  tanglegram()

entanglement(dendlist(dend1, dend4))

# target neurite distr --------------------------------------------------------------------------------------------

# -- syn spread and median along glo
tar_syn_range <- list()
tar_syn_med <- list()
for (j in 1:9) {
  conn <- conn_LC6_tar[[j]]
  conn_range <- matrix(NA, ncol = 2, nrow = 65)
  conn_med <- matrix(NA, ncol = 1, nrow = 65)
  for (k in 1:65) {
    conn_xyz <- xyzmatrix(conn[conn$ID == k, ])
    if (nrow(conn_xyz) > 0) {
      conn_range[k,] <- range(conn_xyz[,'X'])
      conn_med[k] <- median(conn_xyz[,'X'])
    }
  }
  tar_syn_range[[j]] <- conn_range
  tar_syn_med[[j]] <- conn_med
}

# PLOT
df <- data.frame(neu = rep(1:9, times = 1, each = 65), med = matrix(unlist(tar_syn_med), ncol = 1, byrow = T))

dev.new()
ggplot(df, aes(factor(neu), med) ) + 
  # geom_boxplot() +
  # geom_point() +
  geom_jitter(colour = 'black', size = 4, height = 0, width = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', size = 5,position = position_nudge(x = 0.2, y = 0) ) +
  # coord_cartesian(ylim = c(-.5, 1)) +
  theme_minimal_grid() +
  labs(title = "target syn med", x='neu', y='x') 

# PLOT
df <- matrix(ncol = 2, nrow = 0)
for (j in 1:9) {
  df <- rbind(df, tar_syn_range[[j]])
}
df <- data.frame(cbind(rep(1:9, times = 1, each = 65), df))
colnames(df) <- c('neu', 'min','max')

dev.new()
ggplot(df, aes(factor(neu), min) ) + 
  # geom_boxplot() +
  # geom_point() +
  geom_jitter(colour = 'black', size = 4, height = 0, width = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = median, geom = "point",  colour = 'black', size = 5,position = position_nudge(x = 0.2, y = 0) ) +
  # coord_cartesian(ylim = c(-.5, 1)) +
  theme_minimal_grid() +
  labs(title = "target syn range min", x='neu', y='x') 



# -- neurite in cross sections

# DEBUG
nopen3d()
# plot3d(nrefrs_sp, WithNodes = F)
points3d(xyzmatrix(nrefrs_sp$d[ref_ind_sp,]), col = 'red',size = 10)
points3d(o_xyz_m_rs_sp, col = 'grey',size = 5)
points3d(o_xyz_m_rs_sp[ref_ind_sp,], col = 'red',size = 10)
arrow3d(o_xyz_m_rs_sp[ref_ind_sp[j],], o_xyz_m_rs_sp[ref_ind_sp[j],] + xp*5e3, theta = pi / 9, n = 8, col = "red", type = "rotation")
points3d((z_xyz_m_sp[z_ind,]), size =5, col='blue')
arrow3d(o_xyz, o_xyz + zp, theta = pi / 9, n = 8, col = "red", type = "rotation")
planes3d(xp[1],xp[2],xp[3],-d1, alpha=0.3)
planes3d(xp[1],xp[2],xp[3],-d2, alpha=0.3)
points3d(xyz_m[xyz_m_ind,], size =10, col = 'cyan')


yz_ls_target <- list()
ref_d <- c()
dL <- 1000 # half thickness of cross section
for (j in 1:N_ind_sp) {
  # use pc1 to def x-axis
  xp_1 <- o_xyz_m_rs_sp[ref_ind_sp[j],] - o_xyz_m_rs_sp[ref_ind_sp[j]-1,] # +x axis along cable
  xp_pc <- o_xyz_m_rs_sp[(ref_ind_sp[j]-3):(ref_ind_sp[j]+3),] %>%
    prcomp() %>% 
    .[["rotation"]] %>% 
    .[,"PC1"]
  if (xp_1 %*% xp_pc > 0) {
    xp <- xp_pc
  } else {
    xp <- - xp_pc
  }
  
  o_xyz <- o_xyz_m_rs_sp[ref_ind_sp[j],] #origin
  d1 <- (o_xyz + xp*dL) %*% xp #distance to plane
  d2 <- (o_xyz - xp*dL) %*% xp #distance to plane
  
  # use 2 planes
  z_ind <- apply(z_xyz_m_sp, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
  if (sum(z_ind) == 1 ) {
    z_xyz <- z_xyz_m_sp[z_ind,] 
  } else {
    z_xyz <- colMeans(z_xyz_m_sp[z_ind,]) 
  }
  zp <- z_xyz - o_xyz
  zp <- zp - c(zp %*% xp)*xp #+z
  
  ref_d <- c(ref_d, sqrt(sum(zp^2))) #unit length
  
  
  pt_yz <- matrix(NA, ncol = 2, nrow = length(neu_target_rs))
  for (k in 1:nrow(pt_yz)) {
    xyz_m <- xyzmatrix(neu_target_rs[[k]]$d)
    # xyz_m <- xyz_m[xyz_m[,1] > glu_div[j] & xyz_m[,1] < glu_div[j+1],]
    xyz_m_ind <- apply(xyz_m, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0) &
      pointsinside(xyz_m, surf=glo.msh, rval='distance') > -100
                         
    if (sum(xyz_m_ind) > 0) {
      if (sum(xyz_m_ind) == 1) { #if single pt
        xyz <- xyz_m[xyz_m_ind,] 
      } else {
        xyz <- colMeans(xyz_m[xyz_m_ind,]) 
      }
      # dd <- sqrt(rowSums(sweep(xyz_m, MARGIN = 2, STATS = o_xyz, '-')^2))
      # ind <- which.min(dd)
      # xyz <- xyz_m[ind,]
      v <- xyz - o_xyz
      v <- v - c(v %*% xp)*xp 
      ang <- acos(v %*% zp / sqrt(sum(v^2)) / sqrt(sum(zp^2)))
      ang <- ang * sign(v %*% zp)
      mag <- sqrt(sum(v^2)) / sqrt(sum(zp^2)) #normalize to +z
      pt_yz[k,] <- c(cos(ang)*mag, sin(ang)*mag)
    }
  }
  # pt_yz[i_AP[2],] <- c(0,1)
  # pt_yz[i_AP[1],] <- c(0,0)
  yz_ls_target[[j]] <- pt_yz
}

# make df for plotting
df <- matrix(ncol = 3, nrow = 0)
for (j in 1:length(neu_target) ) {
  df <- rbind(df, as.numeric(c(yz_eye[j,], j)) )
  for (k in 1:N_ind_sp) {
    df <- rbind(df, c(yz_ls[[k]][j,], j) )
  }
}
df <- as.data.frame(df)
colnames(df) <- c('y','z','n')
df$n <- factor(df$n)

df_cable_tar <- df

# all 65
m <- 3 #5, 6
df_plt <- df_cable_tar[((N_ind_sp+1)*(seq(1,9)-1)+m),]
df_plt <- cbind(df_plt, factor(seq(1,9)) )
colnames(df_plt)[4] <- 'colcol'
gpl <- ggplot(data = df_plt, aes(x=y, y=z)) +
  geom_point(aes(colour = colcol), size = 15, shape = 20 ) +
  scale_colour_manual(values = pal_tar,
                      breaks = seq(1,9),
                      labels = seq(1,9) ) +
  guides(colour = guide_legend("groups") )+
  geom_text(hjust=0, vjust=0, label=seq(1,9) ) +
  # xlim(xrange) +
  # ylim(yrange) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  # theme_void() +
  # geom_segment(aes(x = 0, y = -2.5, xend = 1/ref_d[m]*1e3, yend = -2.5), size=2, lineend = "round") +
  # annotate("text", x = 0.05, y = -2.8, label = '1 um') +
  # annotate("text", x = -.5, y = -2.8, label = paste('retino-index = ', mean(rt_ind[,m]))) +
  # labs(title = paste('retino-index = ', round(mean(rt_ind[,m]),digits = 2)))
  labs(title = m-1)

windows(record = F, width = 8, height = 8)
# pdf(file = paste("65 - ", m, ".pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl

dev.off()

# target with 4 group
m <- 5 #5, 6

df_plt <- rbind(df_cable_sp[((N_ind_sp+1)*(c(i1,i2,i3,i4)-1)+m),],
                df_cable_tar[((N_ind_sp+1)*(seq(1,9)-1)+m),])
df_plt <- cbind(df_plt, factor(c(rep(1,length(i1)),
                                 rep(2,length(i2)),
                                 rep(3,length(i3)),
                                 rep(4,length(i4)),
                                 seq(101,109))) )
colnames(df_plt)[4] <- 'colcol'
gpl <- ggplot(data = df_plt, aes(x=y, y=z)) +
  geom_point(aes(colour = colcol), size = 15, shape = 20 ) +
  scale_colour_manual(values = c(col4, pal_tar[1:9]),
                      breaks = c(seq(1,4),seq(101,109)),
                      labels = c(seq(1,4),seq(101,109)) ) +
  guides(colour = guide_legend("groups") )+
  # geom_text(hjust=0, vjust=0, label=seq(1,9) ) +
  # xlim(xrange) +
  # ylim(yrange) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  # theme_void() +
  # geom_segment(aes(x = 0, y = -2.5, xend = 1/ref_d[m]*1e3, yend = -2.5), size=2, lineend = "round") +
  # annotate("text", x = 0.05, y = -2.8, label = '1 um') +
  # annotate("text", x = -.5, y = -2.8, label = paste('retino-index = ', mean(rt_ind[,m]))) +
  # labs(title = paste('retino-index = ', round(mean(rt_ind[,m]),digits = 2)))
  labs(title = m-1)

windows(record = F, width = 8, height = 8)
# pdf(file = paste("65 - ", m, ".pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl

dev.off()


# cross section occoupancy  ---------------------------------------------------------------------------------------

# -- x min 1 in LO
ind_bi_ipsi <- c(48,1,4,39,44,14,19,8) #front

ind_bi_ipsi <- c(24,25,56,55,57,9,26,2) #lateral

ind_bi_ipsi <- sample(65, 8) # random


# windows(record = F, width = 10.5, height = 10.5)
# 
# bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
# bd_theta <- seq(1, 180, by = 1)
# xy_bd <- matrix(ncol = 2)
# bd_grid <- expand.grid(bd_phi, bd_theta)
# # range(bd_grid$Var1)
# plot(bd_grid, xlim = c(-20,215), ylim = c(190,-50), type = "n", axes = FALSE, ann = F)
# for (j in 1:length(xy_poly)) {
#   xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
#   if (j %in% i1) {
#     points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[1], cex = 3, pch = 16)
#   }
#   else if (j %in% i2) {
#     points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[2], cex = 3, pch = 16)
#   }
#   else if (j %in% i3) {
#     points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[3], cex = 3, pch = 16)
#   }
#   else if (j %in% i4) {
#     points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=col4_a1[4], cex = 3, pch = 16)
#   }
#   else {
#     points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="grey", cex = 2, pch = 16)
#     # points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="grey", type = '.') 
#   }
#   text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = round(cable_glo_min1_x[j],0), pos = 1, offset = 0.3)
# }
# xy_bd <- xy_bd[-1,]
# hpts_2 <- chull(xy_bd)
# hpts_2 <- c(hpts_2, hpts_2[1])
# xy_bd_chull <- xy_bd[hpts_2,] # hull edge points
# polygon(xy_bd_chull)
# # lines(rbind(c(-11,180), c(-2,180)), lwd = 1)
# # text(-5, 180, labels = expression(paste("9",degree)), pos = 1, offset = 0.3)
# # lines(rbind(c(-11,180), c(-11,171)), lwd = 1)
# # text(-11, 175, labels = expression(paste("9",degree)), pos = 2, offset = 0.2)
# # lines(rbind(c(0,0), c(0,180)), lwd = 1) 
# # text(0, -5, labels = "front", pos = 1, offset = 0)
# # lines(rbind(c(90,0), c(90,180)), lwd = 1)
# # text(90, -5, labels = expression(paste("side 90",degree)), pos = 1, offset = 0)
# # lines(rbind(c(-12,90), c(162,90)), lwd = 1)
# # text(-17, 90, labels = "equator", pos = 1, offset = 0, srt = 90)




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
range_x <- c(356000, 440000)
N_cs <- diff(range_x) / 2 / dL

xp <- c(-1, 0, 0)

yz_ls <- list()
yz_ls_tar <- list()
dd_ls_tar <- list()
ref_d <- c()

for (j in 1:N_cs) {
  
  ii <- o_xyz_m[,1] < (range_x[2] - (j-1)*2*dL) & o_xyz_m[,1] > (range_x[2] - (j)*2*dL)
  o_xyz <- colMeans(o_xyz_m[ii,])
  d1 <- (o_xyz + xp*dL) %*% xp #distance to plane
  d2 <- (o_xyz - xp*dL) %*% xp #distance to plane
  
  # use 2 planes
  z_ind_1 <- apply(z_xyz_m, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
  # y_ind_2 <- sqrt(rowSums(sweep(y_xyz_m_sp, MARGIN = 2, STATS = o_xyz, '-')^2)) < 20000
  z_ind <- z_ind_1
  if (sum(z_ind) == 1 ) {
    z_xyz <- z_xyz_m[z_ind,] 
  } else {
    z_xyz <- colMeans(z_xyz_m[z_ind,]) 
  }
  zp <- z_xyz - o_xyz
  zp <- zp - c(zp %*% xp)*xp # +z
  
  ref_d <- c(ref_d, sqrt(sum(zp^2))) #unit length
  
  pt_yz <- matrix(NA, ncol = 2, nrow = length(LC6))
  for (k in 1:nrow(pt_yz)) {
    xyz_m <- xyzmatrix(LC6_sp_glo_rs[[k]]$d)
    # xyz_m <- xyz_m[xyz_m[,1] > glu_div[j] & xyz_m[,1] < glu_div[j+1],]
    xyz_m_ind_1 <- apply(xyz_m, 1, function(x) (x %*% xp - d1)*(x %*% xp - d2) < 0 )
    # xyz_m_ind_2 <- sqrt(rowSums(sweep(xyz_m, MARGIN = 2, STATS = o_xyz, '-')^2)) < 10000
    xyz_m_ind <- xyz_m_ind_1 #& xyz_m_ind_2 # within 2 planes and 20 um within ref pt
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



# - distance
dd_tar <- matrix(unlist(dd_ls_tar), ncol = 9, byrow = T)
dd_tar <- dd_tar[seq(nrow(dd_tar),1), ]

df <- data.frame(dd_tar[1:80, c(3,4,7)])
df <- cbind(seq(1,80), df)
colnames(df) <- c('pos', 3,4,7)

df <- data.frame(dd_tar[1:80, 3:9])
df <- cbind(seq(1,80), df)
colnames(df) <- c('pos', 3,4,5,6,7,8,9)

dfm <- melt(df, id = 'pos')

tick_x <- rep("", 81)
tick_x[seq(0,80,by = 20)+1] <- seq(0,80,by = 20)


dev.new()
# pdf(file = paste("avg dist to lat LC.pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(x=pos, y=value) ) + 
  geom_point( ) +
  geom_line(aes(colour = factor(variable), group = variable), lwd = 3) +
  # coord_cartesian(ylim = c(-.5, 1)) +
  # theme_minimal_grid() +
  # theme_bw() +
  # scale_colour_manual(values = pnw_palette("Bay",9)[c(3,4,7)] ) +
  scale_colour_manual(values = pnw_palette("Bay",9)[3:9] ) +
  scale_x_continuous(breaks = seq(0,80, by=1), labels = tick_x, expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,8000, by=2000), labels = seq(0,8, by=2)) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color='black'),
        # panel.grid.minor.x = element_line(colour = 'black', size = 0.2),
        axis.ticks.length.x = unit(-.25, "cm"),
        axis.text.x = element_text(margin = margin(t = .3, unit = "cm"))) +
  labs(title = "avg dist to all LC", x='glomerulus axis [um]', y='avg dist [um]') 

dev.off()


# PLOT cartoon
nopen3d()
par3d('windowRect' = c(100,100,1700,1700))
shade3d(glo.msh, col='grey',alpha = 0.5, lit=F)
rgl.viewpoint(userMatrix = rotationMatrix(0/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0))
par3d("FOV" = 0)
segments3d(rbind(c(range_x[1],210000,150000), c(range_x[1]+80*2*dL,210000,150000)), lwd = 5, col = 'black')
segments3d(rbind(c(range_x[1]+80*2*dL,210000,150000),
                 c(range_x[1]+80*2*dL,210000-cos(30/180*pi)*10000,150000-sin(30/180*pi)*10000)), lwd = 5)
segments3d(rbind(c(range_x[2]-20*2*dL,230000-cos(30/180*pi)*12000,130000-sin(30/180*pi)*12000),
                 c(range_x[2]-20*2*dL,230000-cos(30/180*pi)*32000,130000-sin(30/180*pi)*32000)),
           lwd = 5, col = 'black')
# axes3d(c('x','y','z')); title3d('','','x','y','z')

# rgl.snapshot(filename = "80 sectors.png",fmt = "png")


# PLOT cable
m <- 30

df_plt <- rbind(yz_ls[[m]][c(ind_sep,ind_bi_ipsi), ], yz_ls_tar[[m]][c(3,4,7),]) %>% 
  as.data.frame()
df_plt$colcol <- factor(c(rep(1,2),
                          rep(2,length(ind_bi_ipsi)),
                          rep(3,2),
                          rep(4,1)))
colnames(df_plt) <- c('y','z','colcol')
gpl <- ggplot() +
  geom_point(data = df_plt, aes(x=y, y=z, colour = colcol), size = 9, shape = 20 ) +
  scale_colour_manual(values = col4,
                      breaks = c("1", "2", "3","4"),
                      labels = c("ref", "front", "bi", "ipsi") ) +
  guides(colour = guide_legend("groups") )+
  # scale_y_reverse() +
  # xlim(xrange) +
  # ylim(yrange) +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  # theme_void() +
  # geom_segment(aes(x = 0, y = -0.5, xend = 1/ref_d[m]*1e3, yend = -0.5), size=2, lineend = "round") +
  # annotate("text", x = 0.05, y = -0.8, label = '1 um') +
  labs(title = m)

windows(record = F, width = 8, height = 8)
# pdf(file = paste(m, ".pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl

dev.off()



# -- PLOt cross section
j <- 20

ii <- o_xyz_m[,1] < (range_x[2] - (j-1)*2*dL) & o_xyz_m[,1] > (range_x[2] - (j)*2*dL)
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

k <- 7
ii <- pointsinside(xyzmatrix(neu_target_rs[[k]]$d), surf=glo.msh, rval='logical') 
xyz_m <- xyzmatrix(neu_target_rs[[k]]$d)[ii,]
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
df <- cbind(df, c(rep(1,65), rep(2,nrow(pt_yz_k))))
df <- as.data.frame(df)
colnames(df) <- c('x','y','type')

gpl <- ggplot(df) +
  geom_point(aes(x=x, y=y, colour = factor(type), size = type), shape = 16, alpha=0.5 ) +
  scale_colour_manual(values = c('black', 'red'),
                      breaks = c("1", "2"),
                      labels = c("LC6", "LC6G2") ) +
    scale_size(range = c(4,6)) +
  coord_fixed(ratio = 1) +
  theme_void() +
  geom_segment(aes(x = 0, y = -0.5, xend = 1/ref_d[m]*1e3, yend = -0.5), size=2, lineend = "round") +
  annotate("text", x = 0.05, y = -0.8, label = '1 um') +
  labs(title = "cross section occupancy")

windows(record = F, width = 8, height = 8)
# pdf(file = paste("cross section occ.pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl

dev.off()
  
# target neuite vs syn --------------------------------------------------------------------------------------------

break_x <- seq(340000, 450000, length.out = 12)

break_x <- glu_div
break_x <- glo_div_pre

tar_syn_xbin <- matrix(ncol = 9, nrow = length(break_x)-1)
for (j in 1:9) {
  tar_syn_xbin[,j] <- hist(conn_LC6_tar[[j]]$x, breaks = break_x, plot = F)$counts
}
tar_syn_xbin <- cbind(hist(conn_LC6_tar[[1]]$x, breaks = break_x, plot = F)$mids, tar_syn_xbin)

df <- data.frame(tar_syn_xbin)
colnames(df) <- c('mid', seq(1,9))
dfm <- melt(df, id = 'mid')

dev.new()

ggplot(dfm, aes(x=mid, y=value) ) + 
  geom_point(size = 5 ) +
  geom_line(aes(colour = factor(variable), group = variable), lwd = 3) +
  # coord_cartesian(ylim = c(-.5, 1)) +
  theme_minimal_grid() +
  # scale_colour_manual(values = brewer.pal(10, 'RdYlBu')[-5] ) +
  scale_colour_manual(values = pnw_palette("Bay",9) ) +
  labs(title = "target syn histo", x='x', y='count') 


# node hist
neu_target_rs <- resample(neu_target, stepsize = 400)

tar_node_xbin <- matrix(ncol = 9, nrow = length(break_x)-1)
for (j in 1:9) {
  ii <- pointsinside(xyzmatrix(neu_target_rs[[j]]$d), surf=glo.msh, rval='distance') > -500
  xx <- xyzmatrix(neu_target_rs[[j]]$d)[ii,'X']
  tar_node_xbin[,j] <- hist(xx, breaks = break_x, plot = F)$counts
}
tar_node_xbin <- cbind(hist(xx, breaks = break_x, plot = F)$mids, tar_node_xbin)

df <- data.frame(tar_node_xbin)
colnames(df) <- c('mid', seq(1,9))
dfm <- melt(df, id = 'mid')

dev.new()

ggplot(dfm, aes(x=mid, y=value) ) + 
  geom_point( ) +
  geom_line(aes(colour = factor(variable), group = variable), lwd = 3) +
  # coord_cartesian(ylim = c(-.5, 1)) +
  theme_minimal_grid() +
  scale_colour_manual(values = brewer.pal(10, 'RdYlBu')[-5] ) +
  labs(title = "target node histo", x='x', y='count') 


# ratio
df <- data.frame(tar_syn_xbin / tar_node_xbin *2.5) #per 1um
colnames(df) <- c('mid', seq(1,9))
df$mid <- hist(xx, breaks = break_x, plot = F)$mids
dfm <- melt(df, id = 'mid')

dev.new()

# pdf(file = paste("target syn density.pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(x=mid, y=value) ) + 
  geom_point(size = 5 ) +
  geom_line(aes(colour = factor(variable), group = variable), lwd = 3) +
  # coord_cartesian(ylim = c(-.5, 1)) +
  # theme_minimal_grid() +
  scale_colour_manual(values = pnw_palette("Bay",9) ) +
  scale_x_continuous(breaks = seq(range_x[1],range_x[1]+80*2*dL, by=2*dL), limits = c(340000,450000),labels = tick_x, expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0,8000, by=2000), labels = seq(0,8, by=2)) +
  theme(panel.background = element_blank(), 
        axis.line = element_line(color='black'),
        # panel.grid.minor.x = element_line(colour = 'black', size = 0.2),
        axis.ticks.length.x = unit(-.25, "cm"),
        axis.text.x = element_text(margin = margin(t = .3, unit = "cm"))) +
  labs(title = "target syn density", x='glomerulus axis[um]', y='syn count per um') 

dev.off()

