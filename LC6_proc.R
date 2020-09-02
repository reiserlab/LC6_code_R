# analyze LC6 morphology and connectivity, esp retinotopy

# make a lobula layer using LC6 dendrites -----------------------------------------------------------------

# define a plane separating lobula portion of the LC6
a = -0.6; b = 1; c = 1.3; d = -230000

# - get all end points and branch points
neu <- LC6
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

# - fit curved surface
polyfitorder <- 2 # 2nd order surface
gridMargin <- 20000 #add margin on the edge
xyz_node <- data.matrix(xyz_node)
X <- xyz_node[, 1];  Y <- xyz_node[, 2]; Z <- xyz_node[, 3]
fitlm <- lm(Z ~ poly(X, Y, degree = polyfitorder, raw = TRUE)) #linear model fit

# make an underlying grid to interpret the fit as a point set
dx2 <- 500 # grid spacing
dy2 <- 500
xx <- seq(range(X)[1] - gridMargin, range(X)[2] + gridMargin, by = dx2)
yy <- seq(range(Y)[1] - gridMargin, range(Y)[2] + gridMargin, by = dy2)
xygrid <- expand.grid(xx, yy)
xygrid <- setNames(data.frame(xygrid), c('X', 'Y'));
valfit <- predict(fitlm, xygrid) #generate values from the fit
xyz_lm <- cbind(xygrid, valfit)
dist_min <- apply(xyz_lm, 1, function(pt) {min(rowSums(sweep(xyz_node, 2, pt) ^ 2))}) 

ii <- dist_min > 25e6 # for visulization, to make 'm_surf'
valfit_sup <- valfit
valfit_sup[ii] <- NA
m_surf <-  matrix(valfit_sup, nrow = length(xx), ncol = length(yy)) # surface

ii <- dist_min > 12e6 # for later analysis
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

# Figure 5B, LC6 and TM5 with lobula layer and glomerulus
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


# 2d projections -------------------------------------------------------------------
# LC6 dendrites onto the layer, both contour and center-of-mass

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


# - use TM5 to determine the center and meridian  

# too rough
# TM5_lo <- nlapply(TM5, subset, function(x) pointsinside(x, msh))
# TM5_u <- colMeans(xyzmatrix(TM5_lo[[1]]$d)) #upper pt
# TM5_c <- colMeans(xyzmatrix(TM5_lo[[2]]$d)) #upper pt

# tar <- TM5[[1]]
# ng <- as.ngraph(tar)
# distal_points <- igraph::graph.dfs(ng, root=1237, unreachable=FALSE, neimode='out')$order
# proximal_points <- setdiff(igraph::V(ng), distal_points)
# # points3d(xyzmatrix(tar$d[proximal_points,]), col = 'red', size = 5)
# TM5_u <- colMeans(xyzmatrix(tar$d[proximal_points,])) #upper pt
# 
# tar <- TM5[[2]]
# ng <- as.ngraph(tar)
# distal_points <- igraph::graph.dfs(ng, root=30, unreachable=FALSE, neimode='out')$order
# proximal_points <- setdiff(igraph::V(ng), distal_points)
# TM5_c <- colMeans(xyzmatrix(tar$d[proximal_points,])) # central pt

# for old neuron data
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


# - map equator and coord system 
ymax <- max(xy_edge_grid[,2])
ymin <- min(xy_edge_grid[,2])
center_new <- xyz_layer[grid_c,1:2]
x_med_new <- as.numeric(center_new[1])
y_eq_new <- as.numeric(center_new[2])
angR_new <- acos((x_med_new - xyz_layer[grid_u,1])/sqrt(sum((xyz_layer[grid_c,1:2] - xyz_layer[grid_u,1:2])^2)))

# 2d projection with dendrites
ang_2 <- pi/2 - angR_new
rot_2 <- matrix(c(cos(ang_2), sin(ang_2), -sin(ang_2), cos(ang_2)), ncol = 2)
xy_layer_coarse_rot <- sweep(xy_layer_coarse, MARGIN = 2, STATS = c(x_med_new, y_eq_new))
xy_layer_coarse_rot <- sweep(t(rot_2 %*% t(xy_layer_coarse_rot)), 2, c(x_med_new, y_eq_new), '+')

# boundary
xy_edge_grid_rot <- sweep(xy_edge_grid, MARGIN = 2, STATS = c(x_med_new, y_eq_new))
xy_edge_grid_rot <- sweep(t(rot_2 %*% t(xy_edge_grid_rot)), 2, c(x_med_new, y_eq_new), '+')

# Figure 5C, 2d proj of lo and LC6, TM5
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
twig <- neu_lo[[10]]
pp <- sweep(twig$d[,c("X","Y")], MARGIN = 2, STATS = c(x_med_new, y_eq_new))
twig$d[,c("X","Y")] <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new), '+')
plot(twig,  col= "#d7191c", lwd = 2, soma=T, WithNodes = F, add = T)
twig <- neu_lo[[51]]
pp <- sweep(twig$d[,c("X","Y")], MARGIN = 2, STATS = c(x_med_new, y_eq_new))
twig$d[,c("X","Y")] <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new), '+')
plot(twig,  col= "#2c7bb6", lwd = 2, soma=T, WithNodes = F, add = T)
# lines(rbind(xyz_layer[grid_c,1:2], xyz_layer[grid_c,1:2]+(xyz_layer[grid_u,1:2]-xyz_layer[grid_c,1:2])*2), lwd = 3, col = 'cyan')
pp <- sweep(xyz_layer[grid_u,1:2], MARGIN = 2, STATS = c(x_med_new, y_eq_new))
pp <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new), '+')
lines(rbind(c(x_med_new,y_eq_new), c(x_med_new,y_eq_new)+(pp-c(x_med_new,y_eq_new))*2), lwd = 3, col = 'cyan')
points(matrix(c(x_med_new,y_eq_new), ncol =2), pch = 18, col = 'brown', cex = 2)
points(c(x_med_new,y_eq_new)+(pp-c(x_med_new,y_eq_new)), pch = 18, col = 'gold4', cex = 2)
lines(c(300000,310000), c(300000, 300000), col = 'black', lwd = 3)
text(x = 305000, 295000, labels = "10 µm")
# polygon(xy_edge_grid_rot)



# transform to eye coordiante -------------------------------------------------------------------------------------------
# wrap the lobula layer onto a hemisphere, use line-polygon boundary to calculate radial distance

# get edge points for each LC6 projection
xy_poly <- list()
for (j in 1:length(ind_pj)) {
  ls_poly <- mkpoly(xy_pj[[j]]$ashape$edges)
  xy_poly[[j]] <- ls_poly[[1]][,3:4]
}

grid_bdpt <- xy_edge_grid
grid_bdpt <- rbind(grid_bdpt, grid_bdpt[1,])

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


# Figure 5D, LC6 in eye coord
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
    polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "cyan", border = 'black', density = 20, angle = j*2, lwd = 2)
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


# Figure S6A, LC6 overlap
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


# some prep -------------------------------------------------------------------------------------------------------

xy_com_m <- matrix(unlist(xy_com), ncol = 4, byrow = T)[,c(4,3)] #[phi_deg, theta_deg]
xy_com_m_tp <- xy_com_m[,c(2,1)] / 180 * pi # [theta_rad, phi_rad]


#  approx retinotopy along axon  -------------------------------------------------------------------------------

i_AP <- c(4, 10) #long but not i_excl, 4 ante, 10 post, 4 -> 10 as y-axis

# - first get the backbone, and resample

# nopen3d()
# plot3d(LC6[[i_AP[1]]], WithNodes = F)
# points3d(xyzmatrix(LC6[[i_AP[1]]]$d), size = 10)
# plot3d(LC6)
# identify3d(xyzmatrix(LC6[[i_AP[1]]]$d))

# start_xyz <- xyzmatrix(LC6[[i_AP[1]]]$d[4727,]) #starting node for all neu
start_xyz <- c(349390, 182064, 164040) #reloading changes node index

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

# choose ref pt, mid point of glo_div
ref_ind_sp <- c() 
for (j in 1:10) {
  mid_x <- (glo_div[j] + glo_div[j+1])/2
  ii <- which.min((o_xyz_m_rs_sp[175:323,1] - mid_x)^2) + 174 # 175:323 use pts within glo
  ref_ind_sp <- c(ref_ind_sp, ii)
}
# ref_ind_sp <- c(7, 99, ref_ind_sp[c(4,6)]) #all 65
ref_ind_sp <- c(7, 175, ref_ind_sp[c(9, 5)]) #better positions

N_ind_sp <- length(ref_ind_sp)

# cross section projection
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
    # xyz_m <- xyz_m[xyz_m[,1] > glo_div[j] & xyz_m[,1] < glo_div[j+1],]
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
i1 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[30,], '-')^2), decreasing = F)[1:5] #dorsal
i2 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[31,], '-')^2), decreasing = F)[1:(5+0)] #back
i3 <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[12,], '-')^2), decreasing = F)[1:5] #ventral
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

# rgl.snapshot(filename = "corssSection.png",fmt = "png")


# cable plot axon -----------------------------------------------------------------------------------------------

# - Figure 5E, chosen neurons i1, i2 etc

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

# -- 20 chosen
m <- 2 #5, 6
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

windows(record = F, width = 8, height = 8)
# pdf(file = paste(m, ".pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl

dev.off()


# Figure 5F, all 65
m <- 2
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

windows(record = F, width = 8, height = 8)
# pdf(file = paste("65 - ", m, ".pdf",sep = ''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
gpl

dev.off()

# Figure 5H, ordering 65 axon ------------------------------------------------------------------------------------

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


# LC6 glo median, etc --------------------------------------------------------------------------------------------------

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


# # - PLOT, order by median x in LO, 65 color
# getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
# synpc_order_pal <-  getPalette(65)
# cable_glo_x_o <- order(cable_glo_med_x)
# # cable_glo_x_o <- order(cable_glo_mean_x)
# # cable_glo_x_o <- order(cable_glo_mid_x)
# # cable_glo_x_o <- order(cable_glo_min10_x)
# 
# windows(record = F, width = 8, height = 8)
# bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
# bd_theta <- seq(1, 180, by = 1)
# xy_bd <- matrix(ncol = 2)
# bd_grid <- expand.grid(bd_phi, bd_theta)
# plot(bd_grid, ylim = rev(range(bd_grid$Var2)), type = "n", axes = FALSE, ann = F)
# for (j in 1:length(xy_poly)) {
#   xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
#   # points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=synpc_order_pal[cable_glo_x_o[j]], cex = 1.5, pch = 20)
#   points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col=synpc_order_pal[match(j,cable_glo_x_o)], cex = 1.5, pch = 20)
#   # text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = cable_glo_x_o[j], pos = 1, offset = 0.3)
# }
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
# title('x median')


# - Figure 6A, order by median x in LO, 3 colors
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

# - Figure 6B

# 3 color glo, median
nopen3d()
par3d('windowRect' = c(100,100,3100,1100))
rgl.viewpoint(userMatrix = rotationMatrix(30/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0)
              , zoom = 0.3)
# points3d(cable_glo_avg, col = synpc_order_pal_3[cable_glo_x_o[,'gp']], size = 30)
points3d(cable_glo_med, col = synpc_order_pal_3[cable_glo_x_o[,'gp']], size = 30)
# points3d(cable_glo_min1, col = synpc_order_pal_3[cable_glo_x_o[,'gp']], size = 10)
# shade3d(glo.msh, alpha = 0.1)
shade3d(glo.msh, col='cyan',alpha = 0.05)

rgl.pop(type = "light")
rgl.light(theta = 30, phi = 30)
title3d('x median')

# 3 color glo, tip
nopen3d()
par3d('windowRect' = c(100,100,3100,1100))
rgl.viewpoint(userMatrix = rotationMatrix(30/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0)
              , zoom = 0.3)
points3d(cable_glo_min10, col = synpc_order_pal_3[cable_glo_x_o[,'gp']], size = 30)
shade3d(glo.msh, col='cyan',alpha = 0.05)

title3d('cable min 10')


# Figure 6D, index for LC6 median ----------------------------------------------------------------------------------------------------

# -- mean and x-median syn position
# syn_xyz_avg <- matrix(ncol = 3, nrow = length(LC6))
syn_xyz_med_x <- matrix(ncol = 1, nrow = length(LC6))
for (j in 1:length(LC6)) {
  conn <- connectors(LC6[[j]])
  conn_pre <- conn[conn$prepost==0, c('x','y','z')]
  inglo <- pointsinside(conn_pre, glo.msh)
  # inglo <- rowSums(sweep(conn_pre,2,c(a,b,c), '*')) + d < 0
  # syn_xyz_avg[j,] <- colMeans(conn_pre[inglo,])
  syn_xyz_med_x[j,] <- median(conn_pre[inglo,1])
}

# -- 90 deg rot then proj on equator
vr <- c(cos(sum(range(xy_com_m[,1]))/2/180*pi), sin(sum(range(xy_com_m[,1]))/2/180*pi),0)
R_mat <- quaternion3D(vr, ang = 90)

xy_com_m_rot_90 <- sph2cartZ(cbind(1, xy_com_m_tp)) %*% t(R_mat)
xy_com_m_rot_90 <- cart2sphZ(xy_com_m_rot_90)[, c(3,2)]
xy_com_m_rot_90[,1] <- if_else(xy_com_m_rot_90[,1] > pi, xy_com_m_rot_90[,1] - 2*pi, xy_com_m_rot_90[,1])

# -- make df of indices
df <- cbind(seq(1,65),
            RI(sph2cartZ(cbind(1,xy_com_m_tp)), cable_glo_med_x),
            RI(sph2cartZ(cbind(1,xy_com_m_tp)), syn_xyz_med_x),
            RI(sph2cartZ(cbind(1,xy_com_m_tp)), xy_com_m[,1]),
            RI(sph2cartZ(cbind(1,xy_com_m_tp)), xy_com_m_rot_90[,1])
            )

df <- data.frame(df) 
colnames(df) <- c('neu', 'skeleton median', 'syn median', 'AP', 'DV')
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
               geom = "linerange", position = position_nudge(x = 0.2, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun.min = median, geom = "errorbar",  colour = 'black', size = 1, width=.2,
               position = position_nudge(x = 0.2, y = 0) ) +
  coord_cartesian(ylim = c(-.5, 1)) +
  theme_minimal_grid() +
  labs(title = "4 mapping", x='quantity', y='RI') 

dev.off()


ks.test(df[,3], df[,2])

# Figure 6E, 2d LO, with projection onto 1D -----------------------------------------------------------------------

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

# rotation proj ------------------------------------------------------------------------------------------

# -- lin reg for rotations
k_fit <- c()
pval <- c()
rsq <- c()
xy_com_m_tp <- xy_com_m[,c(2,1)] / 180 * pi
cable_glo_med_x_norm <- (cable_glo_med_x - glo_div[1])/ diff(range(glo_div)) 
# rot axis
vr <- c(cos(sum(range(xy_com_m[,1]))/2/180*pi), sin(sum(range(xy_com_m[,1]))/2/180*pi),0)
for (j in 1:36) {
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
# cable_glo_med_x_norm <- (cable_glo_med_x - glo_div[1])/ diff(range(glo_div)) 

# -- Figure 6F, RI for all rot ang
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


df <- data.frame(cbind(seq(1,65), N_swap_rot_index, N_swap_rot_index[,1])) 
colnames(df) <- c('neu', seq(0,180,by = 10))
dfm <- melt(df, id = 'neu')

dev.new()
# pdf(file = "index rotation.pdf", width = 10, height = 4,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  # geom_boxplot() +
  # geom_point() +
  # geom_jitter(colour = 'black', size = 4, height = 0, width = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "linerange", position = position_nudge(x = 0, y = 0),
               colour = 'black', size = 1) +
  stat_summary(fun.min = median, geom = "errorbar",  colour = 'black', size = 1, width=.2,
               position = position_nudge(x = 0.2, y = 0) ) +
  coord_cartesian(ylim = c(-.5, 1)) +
  theme_minimal_grid() +
  labs(title = "index rotation LO", x='Angle', y='RI') 

dev.off()

# -- Figure 6F, 1D vs glo median
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


df <- data.frame(cbind(seq(1,65), N_swap_1D_rot_index, N_swap_1D_rot_index[,1])) 
colnames(df) <- c('neu', seq(0,180,by = 10))
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


# plot 2d
# dev.new()
# plot(xy_com_m_rot)
# points(xy_com_m_rot[c(1,2),], cex = 2, col ='red', pch = 16)
# points(xy_com_m_rot[c(3,4),], cex = 2, col ='blue', pch = 16)

# # PLOT 3d eg. neuron
# which.min(cable_glo_min10_x_norm)

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


# # -- plot regression
# dev.new()
# # par(mfrow=c(4,2))
# 
# # pdf(file = "median reg 2.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
# 
# df <- cbind.data.frame(com_rot_norm, cable_glo_med_x_norm)
# colnames(df) <- c('x','y')
# linReg <- lm(y ~ x, data = df)
# df$fitted <- linReg$fitted.values
# pval <- summary(linReg)$coefficients['x', 'Pr(>|t|)'] #p-value
# rsq <- summary(linReg)$r.squared # R^2
# plot(df[,1:2], pch=16, cex=2, col='black', xaxt='none', yaxt='n',
#      # xlim = c(0,.6), ylim = c(0.3,0.8), 
#      xlim = c(-0.5,0.5), ylim = c(0.3,0.8), 
#      # xlim = c(-1,1), ylim = c(0,1), 
#      cex.axis=2, cex.lab=2, cex.main=2, 
#      main = "lin reg", xlab = "position in LO", ylab = "position in glo")
# # axis(1, at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5), las=1) 
# # axis(2, at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5), las=2) 
# # plot(df[,1:2], pch=16, cex=2, col=4,
# #      xlim = c(20,250), ylim = c(320000,430000), main = j,
# #      cex.axis=2, cex.lab=2, cex.main=2, xaxt='n',yaxt='n', ann=F)
# lines(df[,c(1,3)], lwd=2, col='black')
# text(0, 0.8, labels = bquote(paste("slope = ", .(signif(linReg$coefficients[2],2)))),
#      adj = 0, srt = 0, cex = 1)
# # text(100, 420000, labels = bquote(paste("p-value = ", .(signif(pval,2)))),
# #      adj = 0, srt = 0, cex = 1)
# # text(100, 420000, labels = bquote(paste("p-value = ", .(signif(pval,2)))),
# #      adj = 0, srt = 0, cex = 1)
# 
# dev.off()
# 
# df <- cbind.data.frame(xy_com_m_rot[,1], cable_glo_min10_x)
# colnames(df) <- c('x','y')
# linReg <- lm(y ~ x, data = df)
# df$fitted <- linReg$fitted.values
# pval <- summary(linReg)$coefficients['x', 'Pr(>|t|)'] #p-value
# rsq <- summary(linReg)$r.squared # R^2
# points(df[,1:2], pch=16, cex=2, col=5,
#        xlim = c(20,250), ylim = c(320000,430000), main = j,
#        cex.axis=2, cex.lab=2, cex.main=2,xaxt='n',yaxt='n', ann=F)
# lines(df[,c(1,3)], lwd=2, col='black')
# text(100, 350000, labels = bquote(paste("p-value = ", .(signif(pval,2)))),
#      adj = 0, srt = 0, cex = 1)
# title('median and min 10 reg')


# index of 3 glo sectors--------------------------------------------------------------

# -- Figure 6G
# Function for desaturating colors by specified proportion
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
ii <- LC6_pre_glo_all[,1] > glo_div[4] & LC6_pre_glo_all[,1] < glo_div[5]
points3d(LC6_pre_glo_all[ii,], col = dolphin_col_456[1], size = 10)
ii <- LC6_pre_glo_all[,1] > glo_div[5] & LC6_pre_glo_all[,1] < glo_div[6]
points3d(LC6_pre_glo_all[ii,], col = dolphin_col_456[2], size = 10)
ii <- LC6_pre_glo_all[,1] > glo_div[6] & LC6_pre_glo_all[,1] < glo_div[7]
points3d(LC6_pre_glo_all[ii,], col = dolphin_col_456[3], size = 10)
# ii <- LC6_pre_glo_all[,1] > glo_div[7] & LC6_pre_glo_all[,1] < glo_div[8]
# points3d(LC6_pre_glo_all[ii,], col ='brown')

shade3d(glo.msh, alpha = 0.1,lit = F)
title3d('sectors in synapo glo')
rgl.snapshot(filename = "sector.png",fmt = "png")


# Figure 6G, eg neuron
# sec 4 mesh
ii <- conn_LC6_gloHull[,1] > glo_div[4] & conn_LC6_gloHull[,1] < glo_div[5]
glo4.as <- ashape3d(as.matrix(conn_LC6_gloHull[ii,]), alpha = 10000)
glo4.msh <- as.mesh3d(glo4.as)


nopen3d()
par3d('windowRect' = c(100,100,2100,1100))
rgl.viewpoint(userMatrix = rotationMatrix(30/180*pi,0,1,0) %*% rotationMatrix(150/180*pi,1,0,0)
              , zoom = 0.3)

# shade3d(glo4.msh, alpha = 0.1, col = dolphin_col_456[1], lit = F)
# planes3d(1,0,0, -glo_div[4], col = dolphin_col_456[1], lit=F)
# planes3d(1,0,0, -glo_div[5], col = dolphin_col_456[1], lit=F)

j <- 5
conn <- connectors(LC6[[j]])
conn_pre <- conn[conn$prepost==0, c('x','y','z')]
inglo <- rowSums(sweep(conn_pre,2,c(a,b,c), '*')) + d < 0
conn_pre <- conn_pre[inglo,]
insec <- conn_pre[,1] > glo_div[4] & conn_pre[,1] < glo_div[5]
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
insec <- conn_pre[,1] > glo_div[4] & conn_pre[,1] < glo_div[5]
conn_pre <- conn_pre[insec,]
points3d(conn_pre, size = 40, col = 'red', alpha =0.8)
plot3d(LC6_glo[[j]], col = 'red', lwd = 10, WithNodes = F, alpha =0.5)
ii <- which.min(rowSums(sweep(as.matrix(conn_pre),2,bdot,'-')^2))
segments3d(rbind(bdot, conn_pre[ii,]), lwd = 5, col = 'purple')

# -- Figure 6H

# -- pre-synapse in glo
LC6_pre_glo <- list() #pre-synapse xyz 
LC6_post_glo <- list() #post-synapse xyz 
for (j in 1:length(LC6)) {
  conn <- connectors(LC6[[j]])
  conn_pre <- conn[conn$prepost== 0, c('x','y','z')]
  # inglo <- rowSums(sweep(conn_pre,2,c(a,b,c), '*')) + d < 0
  inglo <- pointsinside(conn_pre, glo.msh)
  LC6_pre_glo[[j]] <- conn_pre[inglo,]
  
  # ii_pre <-  tar$d$PointNo %in% conn[conn$prepost==0,]$treenode_id
  # conn_pre <- xyzmatrix(tar$d[ii_pre,])
  
  conn_post <- conn[conn$prepost== 1, c('x','y','z')]
  # inglo <- rowSums(sweep(conn_post,2,c(a,b,c), '*')) + d < 0
  inglo <- pointsinside(conn_pre, glo.msh)
  LC6_post_glo[[j]] <- conn_post[inglo,]
}

LC6_pre_glo_all <- do.call(rbind.data.frame, LC6_pre_glo)

# check numbers within each compartment
N_pre_comp <- matrix(ncol = 10, nrow = 65)
for (j in 1:65) {
  xx <- LC6_pre_glo[[j]][,1]
  for (k in 1:10) {
    N_pre_comp[j,k] <- sum(xx > glo_div[k] & xx < glo_div[k+1])
  }
}

colSums(N_pre_comp == 0)

# -- order by dist betw pre-syn
# # 4, 5, 6, 8 and combined
# N_swap_pre_ls <- list()
# N_swap_pre_3_ls <- list()
# 
# for (s in c(4,5,6,8)) {
#   LC6_pre_glo_comp <- list() # selected compartment
#   jxyz_3 <- list()
#   for (j in 1:65) {
#     xx <- LC6_pre_glo[[j]][,1]
#     ii <- xx > glo_div[s] & xx < glo_div[s+1] # USE old
#     LC6_pre_glo_comp[[j]] <- LC6_pre_glo[[j]][ii,]
#     # i3 <- xx > glo_div[s-1] & xx < glo_div[s+1+1]
#     # jxyz_3[[j]] <- LC6_pre_glo[[j]][ii,] #check nb compartment
#   }
#   
#   N_swap_pre <- matrix(ncol = 2, nrow = 65)
#   for (j in 1:65) {
#     # oj <- order(rowSums(sweep(xy_com_m, 2, xy_com_m[j,], '-')^2), decreasing = F) 
#     oj <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[j,])), decreasing = F)
#     jxyz <- LC6_pre_glo_comp[[j]] 
#     jswap <- c()
#     for (k in 1:nrow(jxyz)) {
#       dd <- c() #distance
#       for (m in 1:65) {
#         # mdd <- min(rowSums(sweep(jxyz_3[[m]], 2, as.numeric(jxyz[k,]), '-')^2)) #shortes to a neuron
#         mdd <- min(rowSums(sweep(LC6_pre_glo_comp[[m]], 2, as.numeric(jxyz[k,]), '-')^2)) #shortes to a neuron
#         dd <- c(dd, mdd)
#       }
#       
#       order_k <- order(dd[oj])  # re-order
#       # order_k <- order(dd[sample(65)])  # random
#       
#       jswap <- c(jswap, bubble_sort_swapCount(order_k)) #num of swaps
#     }
#     N_swap_pre[j,] <- c(mean(jswap), sd(jswap))
#   }
#   
#   
#   N_swap_pre_ls[[s]] <- N_swap_pre  
#   # N_swap_pre_3_ls[[s]] <- N_swap_pre
# }
# 
# # combine 4,5,6
# LC6_pre_glo_comp <- list() # selected compartment
# for (j in 1:65) {
#   xx <- LC6_pre_glo[[j]][,1]
#   ii <- xx > glo_div[4] & xx < glo_div[7]
#   LC6_pre_glo_comp[[j]] <- LC6_pre_glo[[j]][ii,]
# }
# N_swap_pre_456 <- matrix(ncol = 2, nrow = 65)
# N_swap_pre_ran <- matrix(ncol = 2, nrow = 65)
# for (j in 1:65) {
#   oj <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[j,])), decreasing = F)
#   jxyz <- LC6_pre_glo_comp[[j]] 
#   jswap <- c()
#   jswap_ran <- c()
#   for (k in 1:nrow(jxyz)) {
#     dd <- c() #distance
#     for (m in 1:65) {
#       mdd <- min(rowSums(sweep(LC6_pre_glo_comp[[m]], 2, as.numeric(jxyz[k,]), '-')^2)) #shortes to a neuron
#       dd <- c(dd, mdd)
#     }
#     
#     order_k <- order(dd[oj])  # re-order
#     jswap <- c(jswap, bubble_sort_swapCount(order_k)) #num of swaps
#     
#     order_k_ran <- order(dd[sample(65)])  # random
#     jswap_ran <- c(jswap_ran, bubble_sort_swapCount(order_k_ran)) #num of swaps
#     
#   }
#   N_swap_pre_456[j,] <- c(mean(jswap), sd(jswap))
#   N_swap_pre_ran[j,] <- c(mean(jswap_ran), sd(jswap_ran))
# }
# 
# 
# # - sectors
# swapMax <- 64*(64-1)/2
# df <- data.frame(cbind(seq(1,65), 1 - 2*N_swap_pre_ls[[4]][,1]/swapMax,
#                        1 - 2*N_swap_pre_ls[[5]][,1]/swapMax,
#                        1 - 2*N_swap_pre_ls[[6]][,1]/swapMax,
#                        1 - 2*N_swap_pre_456[,1]/swapMax,
#                        1 - 2*N_swap_pre_ran[,1]/(65*(65-1)/2)) )
# colnames(df) <- c('neu', 'sec4', 'sec5','sec6', 'sec456', 'sec456 random')
# dfm <- melt(df, id = 'neu')
# 
# # Figure 6H
# dev.new()
# # pdf(file = "sector RI.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
# 
# ggplot(dfm, aes(factor(variable), value) ) + 
#   geom_jitter(aes(colour=variable),size = 4,height = 0, width = 0.1) +
#   scale_colour_manual(values = c(dolphin_col_456, 'black', 'black'), guide=FALSE) +
#   stat_summary(fun.min = function(z) { quantile(z,0.25) },
#                fun.max = function(z) { quantile(z,0.75) },
#                geom = "errorbar", position = position_nudge(x = 0.2, y = 0),
#                colour = 'black', size = 1, width = .2 ) +
#   stat_summary(fun = median, geom = "point",  colour = 'black', position = position_nudge(x = 0.2, y = 0),size = 5 ) +
#   coord_cartesian(ylim = c(-0.5, 1)) +
#   theme_minimal_grid() +
#   labs(title = "sector RI", x='quantity', y='RI') 
# 
# dev.off()



# -- dist betw LC6 based on pre-syn
# 4, 5, 6, 8 and combined
N_swap_pre_ls <- list()

for (s in c(4,5,6,8)) {
  LC6_pre_glo_comp <- list() # selected compartment
  for (j in 1:65) {
    xx <- LC6_pre_glo[[j]][,1]
    ii <- xx > glo_div[s] & xx < glo_div[s+1] # USE old
    LC6_pre_glo_comp[[j]] <- LC6_pre_glo[[j]][ii,]
  }
  
  dist_pre <- matrix(0, ncol = 65, nrow = 65)
  for (j in 1:65) {
    jxyz <- LC6_pre_glo_comp[[j]] 
    jdd <- matrix(ncol = 1, nrow = 65)
    for (k in 1:65) {
      if (k != j) {
        dd <- c()
        for (m in 1:nrow(jxyz)) {
          ddm <- sweep(as.matrix(LC6_pre_glo[[k]]), 2, as.matrix(jxyz[m,]),'-')^2 %>%
            rowSums() %>%
            sqrt() %>%
            min()
          dd <- c(dd, ddm)
        }
        dist_pre[k,j] <- mean(dd)
      }
    }
  }
  
  N_swap_pre <- matrix(ncol = 1, nrow = 65)
  for (j in 1:65) {
    oj <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[j,])), decreasing = F)
    order_pre <- order(dist_pre[,j][oj])  # re-order
    N_swap_pre[j] <- bubble_sort_swapCount(order_pre) #num of swaps
  }
  N_swap_pre_ls[[s]] <- N_swap_pre
}


# combine 4,5,6
LC6_pre_glo_comp <- list() # selected compartment
for (j in 1:65) {
  xx <- LC6_pre_glo[[j]][,1]
  ii <- xx > glo_div[4] & xx < glo_div[7]
  LC6_pre_glo_comp[[j]] <- LC6_pre_glo[[j]][ii,]
}

for (j in 1:65) {
  dist_pre <- matrix(0, ncol = 65, nrow = 65)
  for (j in 1:65) {
    jxyz <- LC6_pre_glo_comp[[j]] 
    jdd <- matrix(ncol = 1, nrow = 65)
    for (k in 1:65) {
      if (k != j) {
        dd <- c()
        for (m in 1:nrow(jxyz)) {
          ddm <- sweep(as.matrix(LC6_pre_glo[[k]]), 2, as.matrix(jxyz[m,]),'-')^2 %>%
            rowSums() %>%
            sqrt() %>%
            min()
          dd <- c(dd, ddm)
        }
        dist_pre[k,j] <- mean(dd)
      }
    }
  }
  
  N_swap_pre <- matrix(ncol = 1, nrow = 65)
  N_swap_pre_ran <- matrix(ncol = 1, nrow = 65)
  for (j in 1:65) {
    oj <- order(apply(xy_com_m_tp, 1, function(x) arcLength(x, xy_com_m_tp[j,])), decreasing = F)
    order_pre <- order(dist_pre[,j][oj])  # re-order
    N_swap_pre[j] <- bubble_sort_swapCount(order_pre) #num of swaps
    
    order_pre_ran <- order(dist_pre[,j][sample(65)])  # random
    N_swap_pre_ran[j] <- bubble_sort_swapCount(order_pre_ran) #num of swaps
  }
}


# - sectors
swapMax <- 64*(64-1)/2
df <- data.frame(cbind(seq(1,65), 1 - 2*N_swap_pre_ls[[4]][,1]/swapMax,
                       1 - 2*N_swap_pre_ls[[5]][,1]/swapMax,
                       1 - 2*N_swap_pre_ls[[6]][,1]/swapMax,
                       1 - 2*N_swap_pre[,1]/swapMax,
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

# LC6-LC6 syn-weighted distance -----------------------------------------------------------------------------------

# - dist to a neuron avg over all partners weighted by num of conn
dist_Nconn %<>% as_tibble() %>%
  mutate(distxN = dist_com*Nconn_glo) %>%
  as.data.frame()
LC6_avgdist_wt <- matrix(ncol = 2, nrow = length(LC6))
LC6_avgdist <- c()
LC6_nnbdist <- c()
# N_nb <- c()
for (j in 1:length(LC6)) {
  tib_tmp <- dist_Nconn %>%
    as_tibble() %>%
    mutate(bool = from == j | to == j) %>%
    filter(bool == TRUE) %>%
    # transmute(Nconn_glo, distxN) %>%
    data.matrix()
  # N_nb <- c(N_nb, sum(tib_tmp[,"distxN"] != 0))
  LC6_avgdist[j] <- sum(tib_tmp[,"dist_com"])/dim(tib_tmp)[1]
  LC6_nnbdist[j] <- min(tib_tmp[,"dist_com"])
  LC6_avgdist_wt[j,] <- c(j, sum(tib_tmp[,"distxN"])/sum(tib_tmp[,"Nconn_glo"]))
}
# avg_nb <- mean(N_nb)


# dist_nonzero <- dist_Nconn %>%
#   filter(distxN != 0) %>%
#   transmute(Nconn_glo, distxN) %>%
#   data.matrix()

# - randomize connection partners
mat_avgdist_rand <- matrix(ncol = 2, nrow = length(LC6))
for (j in 1:length(LC6)) {
  #  maintain connection num for given LC6
  ii_rand <- sample.int(length(LC6)-1)
  tib_tmp <- dist_Nconn %>%
    as_tibble() %>%
    mutate(bool = from == j | to == j) %>%
    filter(bool == TRUE) %>%
    select(dist_com, Nconn_glo) %>%
    mutate(dist_com = dist_com[ii_rand]) %>%
    mutate(distxN = dist_com*Nconn_glo) %>%
    as.data.frame()
  mat_avgdist_rand[j,] <- c(j, sum(tib_tmp[,"distxN"])/sum(tib_tmp[,"Nconn_glo"]))
}


# Figure 5SB, syn-wt LC6 distance
dat_ggplot <- data.frame(rbind(LC6_avgdist_wt, mat_avgdist_rand))
gp  <- factor(c(rep(1,length(LC6)),rep(2,length(LC6))), labels = c("EM data","Randomized"))
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



# glomerulus dolphin compartment ----------------------------------------------------------------------------------

N_gp <- 11
getPalette <- colorRampPalette(brewer.pal(9, "RdYlBu"))
dolphin_col <- getPalette(N_gp - 1)


# 
# # for making ahull
# conn_LC6 <- LC6_pre # only pre-synapses to divide the dolphin
# glo_div <- quantile(conn_LC6$x, probs = seq(0,1,length.out = N_gp)) #separate into divisions of equal pts based on x
# glo_div[1] <- glo_div[1]-1

conn_LC6 <- LC6_pre

conn_LC6 <- conn_LC6 %>%
  as_tibble() %>%
  mutate(gp_x = 0)
for (j in 1:length(glo_div)) {
  conn_LC6 %<>% as_tibble() %>%
    mutate(gp_x = gp_x + (x>glo_div[j]))
}

# assign group
conn_LC6LC6 <- LC6LC6
conn_LC6LC6 <- conn_LC6LC6 %>%
  as_tibble() %>%
  mutate(gp_x = 0)
for (j in 1:length(glo_div)) {
  conn_LC6LC6 %<>% as_tibble() %>%
    mutate(gp_x = gp_x + (x>glo_div[j]))
}

ID_wt <- list() #neuron skid with synapse weight in each division
for (j in 1:length(glo_div)) {
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
for (j in 1:(length(glo_div)-1)) {
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$Z = 0
  colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(LC6)) {
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


# # plot, all panels, 
# dev.new()
# plot_grid(plvl[[1]],
#           plvl[[2]],
#           plvl[[3]],
#           plvl[[4]],
#           plvl[[5]],
#           plvl[[6]],
#           plvl[[7]],
#           plvl[[8]],
#           plvl[[9]],
#           plvl[[10]])


# - Figure S6C, LC6 syn in glo 10 compartments
# color contours
xy_bd_chull_df <- as.data.frame(xy_bd_chull)
gg_cont <- ggplot(grid_Gaussian_sum, aes(X, Y, z = Z)) +
  coord_fixed(ratio = 1) +
  scale_y_reverse()
for (j in 1:(length(glo_div)-1)) {
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
# conn_LC6_hull <- as.data.frame(conn_LC6)[-c(5646,5679,5640,5631,5654,5652,1642, 1641),]
conn_LC6_hull <- as.data.frame(conn_LC6)[-c(1642, 1641),]
gp_x <- factor(conn_LC6_hull[,"gp_x"], labels = seq(1,N_gp-1, by = 1))
conn_LC6_hull[,"gp_x"] <- gp_x
ii_ahull <- ashape(unique(conn_LC6_hull[, c("x","y")]), alpha = 7000)$edges[, 3:6]
ii_ahull <- as.data.frame(ii_ahull)
# make separate alpha hulls
ii_ahull_ls <- list()
for (j in 1:(length(glo_div)-1)) {
  ii_ahull_ls[[j]] <- ashape(unique(conn_LC6_hull[conn_LC6_hull$gp_x == j, c("x","y")]), alpha = 7000)$edges[, 3:6]
  ii_ahull_ls[[j]] <- as.data.frame(ii_ahull_ls[[j]])
}

# com of each division
range_syn <- matrix(nrow = length(glo_div)-1, ncol = 4)
for (j in 1:(length(glo_div)-1)) {
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
for (j in 1:(length(glo_div)-1)) {
  N_syn_LC6[j] <- sum(conn_LC6LC6$gp_x == j)
}

# dolphin in colors
conn_LC6LC6 <- as.data.frame(conn_LC6LC6)
gp_x <- factor(conn_LC6LC6[,"gp_x"], labels = seq(1,N_gp-1, by = 1))
conn_LC6LC6[,"gp_x"] <- gp_x

pglo <- ggplot(conn_LC6LC6) + 
  geom_point(aes(x = x, y = y, colour = gp_x), shape = 16)
for (j in 1:(length(glo_div)-1)) { #add ahulls
  pglo <- pglo + 
    geom_segment(data = ii_ahull_ls[[j]], aes(x = x1, y = y1, xend = x2, yend = y2)) +
    annotate("text", x = mean(range_syn[j,c(1,2)]), y = mean(range_syn[j,c(3,4)]), label = paste(sprintf("%.1f %%", 100*N_syn_LC6[j]/sum(N_syn_LC6))))
}
pglo <- pglo +
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
  draw_plot(pglo + theme(legend.justification = "bottom"), 0, 0.51, 1, 0.5) +
  draw_plot(gg_cont + theme(axis.title = element_text()), 0, 0, 0.5, 0.5)



# Figure S6D, 10 individual compartments
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
  for (k in 1:length(LC6)) {
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



# stats test ------------------------------------------------------------------------------------------------------

# KS test
ks.test(LC6_avgdist_wt, mat_avgdist_rand)

# Mann-Whitney
wilcox.test(LC6_avgdist_wt, mat_avgdist_rand, alternative = 'less')
