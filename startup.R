# codes used for analyze LC6 and downstream neurons
# !!! open with encoding UTF-8


#  load libraries and startup -----------------------------------------------------------------------------------------

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
library(Hmisc)
library(colorspace)   ## hsv colorspace manipulations


setwd("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/LC6_code_R")
# setwd("C:/Users/zhaoa/Dropbox (HHMI)/LC6 downstream paper new/Code")
# setwd("./") # set to this code's directory
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



# functions  ---------------------------------------------------------------------

# -- generate polygon from set of points
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

# - cross product
cross3D <- function(a, b){
  prod <- matrix(c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1]), ncol = 1)
  return(prod)
}

# - arclength
arcLength <- function(p1, p2) {
  p1_xyz <- c(sin(p1[1])*cos(p1[2]), sin(p1[1])*sin(p1[2]), cos(p1[1]))
  p2_xyz <- c(sin(p2[1])*cos(p2[2]), sin(p2[1])*sin(p2[2]), cos(p2[1]))
  c2 <- sum((p1_xyz - p2_xyz)^2)
  arc_ang <- acos((2-c2)/2)
  return(arc_ang*1)
}

# - bubble sort function, x is member of Sn, symmetric group of n letters from an ordered set
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

# - alt, find the nth letter and insert at nth position, works for set(seq(1,n))
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

# - position after rotation:  pt1 <- pt0 %*% t(R_mat)
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

# - theta is from z-axis
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

# - RI, use Euclidean dist
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

# # load neurons ----------------------------------------------------------------------------------------------------
# 
# # -- LC6 neurons
# anno_LC6 <- catmaid_query_by_annotation("^LC6 neuron$")
# anno_RHS <- catmaid_query_by_annotation("^LO-R")
# anno_rightLC6 <- anno_LC6[anno_LC6[,1] %in% anno_RHS[,1],]
# LC6_skid <- anno_rightLC6[,"skid"]
# neu <-  read.neurons.catmaid(LC6_skid, .progress='text')
# LC6 <- neu
# 
# altTract <- c(8,14,19,39)

# -- target neurons
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

# neu_target <- c(neu_biL, neu_biR, neu_ipsi)
# 
# # -- TM5
# neu_JSON <- fromJSON(file = "data/Tm5_LC6 mapping.json")
# neu_skid <- c()
# for (j in 1:length(neu_JSON)) {
#   neu_skid[j] <- neu_JSON[[j]]$skeleton_id
# }
# TM5 = read.neurons.catmaid(neu_skid, .progress='text')
# 
# # -- load glomerulus volumne mesh 
# glo_vol <- catmaid_get_volume("v14.LC6_glomerulus_preSyn_R")
# 
# # -- whole brain mesh
# v14 <- catmaid_get_volume(439, rval = 'mesh3d')  


# load old neuron data --------------------------------------------------------------------------------------------

load("neurons_20200828.RData")

# -- load glomerulus volumne mesh 
glo_vol <- catmaid_get_volume("v14.LC6_glomerulus_preSyn_R")

# -- whole brain mesh
v14 <- catmaid_get_volume(439, rval = 'mesh3d')  


# connections -----------------------------------------------------------------------------------------------------

# - all the connectors on LC6
neu <- LC6
LC6_pre <- data.frame()
LC6_post <- data.frame()
conn_glo <- list()
conn_pre <- list()
conn_post <- list()
for (j in 1:length(neu)) {
  tar <- neu[[j]]
  conn_glo[[j]] <- tar$connectors %>% 
    as_tibble() %>%
    mutate(glo = a*x + b*y + c*z + d) %>%
    filter(glo < 0)
  conn_pre[[j]] <- conn_glo[[j]] %>% 
    filter(prepost == 0) %>%
    select(treenode_id, connector_id, x, y, z) %>%
    as.data.frame()
  conn_post[[j]] <- conn_glo[[j]] %>% 
    filter(prepost == 1) %>%
    select(treenode_id, connector_id, x, y, z) %>%
    as.data.frame()
  LC6_pre <- rbind(LC6_pre, cbind(rep(j, dim(conn_pre[[j]])[1]),conn_pre[[j]]))
  LC6_post <- rbind(LC6_post, cbind(rep(j, dim(conn_post[[j]])[1]),conn_post[[j]]))
}
colnames(LC6_pre)[1] <- "ID"
colnames(LC6_post)[1] <- "ID"


# - betw LC6
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
LC6LC6 <- as_tibble(LC6LC6) %>%
  mutate(glo = a*x + b*y + c*z + d) %>%
  filter(glo < 0) %>%
  select(-glo) %>%
  as.data.frame()
colnames(LC6LC6) <- c('pre_ind','post_ind','pre_skid','post_skid','conn_id','pre_node_id','post_node_id','x','y','z')
rownames(LC6LC6) <- seq(1,dim(LC6LC6)[1])

# - look at pairwise distance vs pairwise connection
a1 = -0.6; b1 = 1; c1 = 1.3; d1 = -200000 #separate LO and glomerulus
dist_Nconn <- matrix(ncol = 7, nrow = 0) #pairwise dist and connections
for (j in 1:(length(conn_pre)-1)) {
  for (k in (j+1):length(conn_pre)) {
    dist_d <- dist(rbind(xyz_com[[j]], xyz_com[[k]]))
    conn_fromto <- catmaid_get_connectors_between(pre_skids = neu[[j]]$skid, post_skids = neu[[k]]$skid)
    if (is.null(conn_fromto)) {
      fromto_LO <- 0
      fromto_glo <- 0
    } else {
      mat_conn <- data.matrix(conn_fromto[,c("connector_x", "connector_y", "connector_z")])
      fromto_LO <- sum(mat_conn %*% c(a1,b1,c1) + d1 > 0)
      fromto_glo <- dim(conn_fromto)[1] - fromto_LO
    }
    conn_tofrom <- catmaid_get_connectors_between(pre_skids = neu[[k]]$skid, post_skids = neu[[j]]$skid)
    if (is.null(conn_tofrom)) {
      tofrom_LO <- 0
      tofrom_glo <- 0
    } else {
      mat_conn <- data.matrix(conn_tofrom[,c("connector_x", "connector_y", "connector_z")])
      tofrom_LO <- sum(mat_conn %*% c(a1,b1,c1) + d1 > 0)
      tofrom_glo <- dim(conn_tofrom)[1] - tofrom_LO
    }
    dist_Nconn_tmp <- matrix(c(j,k,dist_d, fromto_LO, fromto_glo, tofrom_LO, tofrom_glo), nrow = 1)
    dist_Nconn <- rbind(dist_Nconn, dist_Nconn_tmp)
  }
}
colnames(dist_Nconn) <- c("from","to","dist_com","fromto_LO","fromto_glo","tofrom_LO","tofrom_glo")
dist_Nconn %<>% 
  as_tibble() %>%
  mutate(Nconn_glo = tofrom_glo + fromto_glo) %>%
  as.data.frame()
colSums(dist_Nconn)
dist_com_mean <- mean(dist_Nconn$dist_com)
Nconn_glo_tot <- sum(dist_Nconn$Nconn_glo)


# glo sectors and mesh,  glo LC6 --------------------------------------------------------------------------------------------

N_gp <- 11
glo_div <- quantile(LC6_pre$x, probs = seq(0,1,length.out = N_gp)) #separate into divisions of equal pts based on x
glo_div[1] <- glo_div[1]-1

# make a glo mesh
# conn_LC6_gloHull <- LC6_pre[-c(1642, 1641, 6649, 6643, 6647, 6686, 6679, 6681, 6685,
#                                 2942, 2943, 6650, 6651, 6652, 6655, 6646, 4482, 4483,
#                                 4485, 4512, 4480, 4479),
#                              c('x','y','z')]
conn_LC6_gloHull <- LC6_pre[-c(1644, 1643, 6649, 6658, 6655, 6653, 6654, 6652, 6650,
                               6646, 6689, 2946, 2947, 6682, 6684, 6688, 4485, 4486,
                               4489, 4488, 4491, 4405),
                            c('x','y','z')]
glo.as <- ashape3d(as.matrix(conn_LC6_gloHull), alpha = 8000)
glo.msh <- as.mesh3d(glo.as)

# filter for neurites within the mesh
LC6_glo <- nlapply(LC6, subset, function(x) pointsinside(x, surf=glo.msh, rval='distance') > -5000)


nopen3d()
plot3d(LC6_glo)
shade3d(glo.msh, alpha = 1)
# points3d(conn_LC6_gloHull)

# resample
LC6_glo_rs <- resample(LC6_glo, stepsize = 400)


#  target syn -----------------------------------------------------------------------------------------------------

a1 = -0.6; b1 = 1; c1 = 1.3; d1 = -200000 #separate LO and glomerulus
conn_target <- list()
for (j in 1:length(neu_target)) {
  tb_conn <- matrix(ncol = 7)
  for (k in 1:length(neu)) {
    conn_fromto <- catmaid_get_connectors_between(pre_skids = neu_target[[j]]$skid, post_skids = neu[[k]]$skid)
    if (is.null(conn_fromto)) {
      fromto_LO <- 0
      fromto_glo <- 0
    } else {
      mat_conn <- data.matrix(conn_fromto[, c("connector_x", "connector_y", "connector_z")])
      fromto_LO <- sum(mat_conn %*% c(a1, b1, c1) + d1 > 0)
      fromto_glo <- dim(conn_fromto)[1] - fromto_LO
    }
    conn_tofrom <- catmaid_get_connectors_between(pre_skids = neu[[k]]$skid, post_skids = neu_target[[j]]$skid)
    if (is.null(conn_tofrom)) {
      tofrom_LO <- 0
      tofrom_glo <- 0
    } else {
      mat_conn <- data.matrix(conn_tofrom[, c("connector_x", "connector_y", "connector_z")])
      tofrom_LO <- sum(mat_conn %*% c(a1, b1, c1) + d1 > 0)
      tofrom_glo <- dim(conn_tofrom)[1] - tofrom_LO
    }
    Nconn_glo <- fromto_glo + tofrom_glo
    tb_conn_tmp <- matrix(c(j, k, fromto_LO, fromto_glo, tofrom_LO, tofrom_glo, Nconn_glo), nrow = 1)
    tb_conn <- rbind(tb_conn, tb_conn_tmp)
  }
  conn_target[[j]] <- as.data.frame(tb_conn[-1,])
  colnames(conn_target[[j]]) <- c("target","LC6","fromto_LO","fromto_glo","tofrom_LO","tofrom_glo", "Nconn_glo")
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


# - LC6 and target
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
  
  for (k in 1:(length(glo_div)-1)) {
    tmp %<>% as_tibble() %>%
      mutate(gp_x = gp_x + (x>glo_div[k]))
  }
  
  ID_wt <- list() #neuron skid with synapse weight in each division
  for (k in 1:length(glo_div)) {
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



# ephy data -------------------------------------------------------------------------------------------------------

library(cowplot)
library(R.matlab)
expBi3 <- readMat("data/ss825_bi_RF_indiv_tseries.mat") #bi
expBi3 <- expBi3[[1]]
expBi <- readMat("data/ss825_mean_RFmap_n=4.mat") #bi
# expBi <- readMat("C:/Users/zhaoa/Dropbox (HHMI)/sync_userA/Documents/ReiserGroup/p_LC6/dataFromMai/ss825_mean_AUC_RFmap_n=4.mat") #bi
expBi <- as.matrix(expBi[[1]])


expIpsi3 <- readMat("data/ss2036_ipsi_RF_indiv_tseries.mat") #ipsi
expIpsi3 <- expIpsi3[[1]]
expIpsi <- readMat("data/ss2036_mean_RFmap_n=7.mat") #ipsi
expIpsi <- as.matrix(expIpsi[[1]])

# indexing
ind_mai <- c(t(matrix(seq(1:98), byrow = F, ncol = 14)))
tar_pal <- brewer.pal(4,"RdYlBu")[c(1,3,4,2)]


# x2 <- seq(-18 - 2.25, 99 - 2.25, by = 9) # -18 to 99, with sim data
# y2 <- seq(54, 108,by = 9)
# xygrid2 <- expand.grid(x2, y2) # looming center position 

# loom position from Matt
loom_theta_mat <- read.csv('data/loom_center_theta.txt', sep = ',', header = F) / pi * 180
loom_theta_mat <- 180 - loom_theta_mat
loom_theta_mat <- loom_theta_mat[seq(7,1),]
loom_phi_mat <- read.csv('data/loom_center_phi.txt', sep = ',', header = F) / pi * 180
loom_phi_mat <- loom_phi_mat[seq(7,1),]

loom_theta <- melt(t(as.matrix(loom_theta_mat)))$value 
loom_phi <- melt(t(as.matrix(loom_phi_mat)))$value 

xygrid2 <- cbind(loom_phi, loom_theta) 

# shim
shim_xy <- read.csv('data/shim.csv', sep = ',', header = F)
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



# user stats ------------------------------------------------------------------------------------------------

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

