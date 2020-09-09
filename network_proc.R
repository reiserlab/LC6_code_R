# pairwise pre-syn distance --------------------------------------------------------------------

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

dd_pair <- as.data.frame(dd_pair)
colnames(dd_pair) <- c('j','k', 'dd')

mat_names <- seq(1,65)
df <- dd_pair

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



#  graph ----------------------------------------------------------------------------------------------------------

library(igraph)
library(dendextend)


df <- dd_pair
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





#  hierarchical clustering ----------------------------------------------------------------------------------------


# - pre-syn dist
df <- dd_pair
colnames(df) <- c('from','to', 'weight')
g <- graph.data.frame(df, directed = F)
gc <- cluster_fast_greedy(g)
modularity(gc)

dev.new()
plot(gc, g)
plot_dendrogram(gc)


# - user hclust()

# -- dist in LO
d_LO <- dist(sph2cartZ(cbind(1,xy_com_m_tp)))
hc_LO <- hclust(d_LO, method = "ward.D2")

dev.new()
plot(hc_LO)


# -- pre-syn dist
# # make dist from pre-syn dist
# mat <- matrix(ncol = 65, nrow = 65)
# ii <- 1
# for (j in 1:64) {
#   for (k in (j+1):65) {
#     mat[k,j] <- dd_pair[ii,3]
#     ii <- ii + 1
#   }
# }
# mat <- as.dist(mat)
df <- dd_pair
colnames(df) <- c('from','to', 'weight')
g <- graph.data.frame(df, directed = F)
adjM <- igraph::get.adjacency(g, sparse = F, attr = 'weight')
mat <- as.dist(adjM)
attributes(mat)$Labels <- NULL
hc_glo <- hclust(mat, method = "ward.D2")

dev.new()
plot(hc_glo)


# - LC6 connections
LC6LC6 <- matrix(ncol = 10, nrow = 0)
for (j in 1:length(neu)) {
  # tar_pre <- neu[[j]]
  for (k in 1:length(neu)) {
    # tar_post <- neu[[k]]
    cft <- catmaid_get_connectors_between(pre_skids = neu[[j]]$skid, post_skids = neu[[k]]$skid)
    if (!is.null(cft)) {
      cft_comb <- cbind(j,k, cft[,1:8])
      LC6LC6 <- rbind(LC6LC6,cft_comb)
    }
  }
}
LC6LC6 <- as.data.frame(LC6LC6)
colnames(LC6LC6) <- c('pre_ind','post_ind','pre_skid','post_skid','conn_id','pre_node_id','post_node_id','x','y','z')
rownames(LC6LC6) <- seq(1,dim(LC6LC6)[1])

conn <- LC6LC6[, c(1,2)]
conn_g <- graph.data.frame(conn, directed = F)
# conn_g2 <- graph_from_edgelist(as.matrix(conn), directed = T)
# adjM <- get.adjacency(conn_g, sparse = F)
adjM <- as_adjacency_matrix(conn_g, sparse = F)

# conn_g3 <- graph_from_adjacency_matrix(adjM, weighted = T, mode = "plus")
# adjM3 <- as_adjacency_matrix(conn_g3, sparse = F)

# dev.new()
# plot(conn_g3, edge.label=E(conn_g3)$weight, vertex.size =4)
     # layout = layout_in_circle)

mat3 <- as.dist(adjM)
attributes(mat3)$Labels <- NULL

# mat3 <- 1/(mat3+1) # 1/(x+1)

mat3 <- max(mat3) - mat3

  # hc_conn <- hclust(mat3, method = "complete")
hc_conn <- hclust(mat3, method = "ward.D2")

dev.new()
plot(hc_conn)

# - cp 
dend1 <- as.dendrogram(hc_LO)
dend2 <- as.dendrogram(hc_glo)
dend3 <- as.dendrogram(hc_conn)


dev.new()
dendlist(dend2, dend3) %>% 
  untangle(method = "step1side") %>% 
  tanglegram()

entanglement(dendlist(dend1, dend4))
