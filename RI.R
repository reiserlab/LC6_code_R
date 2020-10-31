# Kendall tau distance, https://en.wikipedia.org/wiki/Kendall_tau_distance
# No. swaps = No. inversion

# Kendall rank corr coeff, https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient
# ordinary generating func, https://math.stackexchange.com/questions/245018/bubble-sorting-question

# on a grid ----------------------------------------------------------------------------------------------------

# Figure 5H and 5S1A
Nr <- 11
Nc <- 6
y <- seq(-Nr/2+0.5, Nr/2-0.5, length.out = Nr)
x <- seq(-Nc/2+0.5, Nc/2-0.5, length.out = Nc)
pt0 <- expand.grid(x, y) + cbind(runif(Nr*Nc,min = -1,max = 1)*0.1, runif(Nr*Nc,min = -1,max = 1)*0.1)
pt0 <- as.matrix(pt0)

# -- pal
getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
col_66 <- rev(rev(getPalette(38))[1:(Nr*Nc/2)])
getPalette <- colorRampPalette(brewer.pal(9, "Reds"))
col_66 <- c(col_66, rev(getPalette(38))[1:(Nr*Nc/2)])

#  mapping --------------------------------------------------------------------------------------------------------

# - swap points
pt1 <- pt0
pt1[c(63,39),] <- pt1[c(39,63),]
pt1[c(63,12),] <- pt1[c(12,63),]

# - random
pt1 <- pt0[sample(66),]


# - swap rows
pt1 <- pt0
pt1[(1:8)+2*8,] <- pt0[(1:8)+8*7,]
pt1[(1:8)+7*8,] <- pt0[(1:8)+2*8,]


# - 1D proj
pt1 <- pt0
pt1[,1] <- pt1[,1] / 2
pt1 <- pt1 + cbind(runif(Nr*Nc,min = -1,max = 1)*0.1, runif(Nr*Nc,min = -1,max = 1)*0.1)


#  - cos
pt1 <- pt0
pt1[,2] <- cos(pt1[,2]/2)


# - cylinder
pt1 <- cbind(cos(pt0[,1]/8*2*pi), sin(pt0[,1]/8*2*pi), pt0[,2]/7*2*pi)

nopen3d()
points3d(pt1, col = col_66, size = 10)
axes3d(c('x','y','z')); title3d('','','x','y','z')
# rgl.bbox(color="gray90", alpha = 1, 
#          # emission="#333377", specular="#333377", shininess=5,
#          xlen = 0, ylen = 0, zlen = 0)


# - Gaussian
z <- 10 * exp( - rowSums(pt0^2) / 2 / 4^2)
pt1 <- cbind(pt0,z)
colnames(pt1) <- c('x','y','z')

nopen3d()
points3d(pt1, col = col_66, size = 10)
rgl.bbox(color="gray90", alpha = 1, 
         # emission="#333377", specular="#333377", shininess=5,
         xlen = 0, ylen = 0, zlen = 0)


# - jittering
pt1 <- expand.grid(x, y) + cbind(runif(Nr*Nc,min = -1,max = 1)*1, runif(Nr*Nc,min = -1,max = 1)*1)
pt1 <- as.matrix(pt1)


# - triangle
# pt1 <- expand.grid(x, y) + cbind(runif(64,min = -1,max = 1)*0.5, runif(64,min = -1,max = 1)*0.5)
pt1 <- pt0
for (j in 1:Nr) {
  pt1[seq(1,Nc)+(j-1)*Nc,1] <- pt1[seq(1,Nc)+(j-1)*Nc,1]*0.05*j
}


# - local jitter
pt1 <- pt0
ii <- pt0[,1] > -1.5 & pt0[,1] < 1.5 & pt0[,2] > -2.5 & pt0[,2] < 2.5
pt1[ii,] <- pt0[ii,][sample(sum(ii)),] + cbind(runif(sum(ii),min = -1,max = 1)*1, runif(sum(ii),min = -1,max = 1)*1)


# - block jitter
pt1 <- pt0
for (j in seq(-2.5, 6.5, by = 3)) {
  ii <- pt0[,2] > (j-3) & pt0[,2] < j & pt0[,1] < 0
  pt1[ii,] <- pt0[ii,] + cbind(runif(sum(ii),min = -1,max = 1)*1, runif(sum(ii),min = -1,max = 1)*1)
  ii <- pt0[,2] > (j-3) & pt0[,2] < j & pt0[,1] > 0
  pt1[ii,] <- pt0[ii,] + cbind(runif(sum(ii),min = -1,max = 1)*1, runif(sum(ii),min = -1,max = 1)*1)
}


# - block random
pt1 <- pt0
for (j in seq(-2.5, 6.5, by = 3)) {
  ii <- pt0[,2] > (j-3) & pt0[,2] < j & pt0[,1] < 0
  pt1[ii,] <- pt0[ii,][sample(sum(ii)),]
  ii <- pt0[,2] > (j-3) & pt0[,2] < j & pt0[,1] > 0
  pt1[ii,] <- pt0[ii,][sample(sum(ii)),]
}



# plot ------------------------------------------------------------------------------------------------------------

df0 <- pt0 %>% cbind(.,seq(1,Nr*Nc)) %>% data.frame()
colnames(df0) <- c('x','y', 'no')

dev.new()
# pdf(file = "orig.pdf", width = 8, height = 8,pointsize=12,family="Courier", useDingbats = F)
ggplot() +
  geom_point(data = df0, aes(x,y, colour = factor(no)), size =10,  shape = 16) + 
  scale_color_manual(values = col_66) +
  # labs(title = mat_names[j]) +
  theme_void() +
  coord_fixed(ratio = 1) +
  theme(legend.position = "none")

dev.off()


df1 <- pt1 %>% cbind(.,seq(1,Nr*Nc)) %>% data.frame()
colnames(df1) <- c('x','y', 'no')

dev.new()
# pdf(file = paste("block jitter, ", round(mean(RI(pt0, pt1)),2), " .pdf"), width = 8, height = 8,pointsize=12,family="Courier", useDingbats = F)
ggplot() +
  geom_point(data = df1, aes(x,y, colour = factor(no)), size =10,  shape = 16) + 
  scale_color_manual(values = col_66) +
  # labs(title = mat_names[j]) +
  theme_void() +
  coord_fixed(ratio = 1) +
  theme(legend.position = "none")

dev.off()

mean(RI(pt0, pt1))

# rot sweep  ------------------------------------------------------------------------------------------------------
# Figure 5S1B

# -- RI for all rot ang
RI_2D_rot <- matrix(ncol = 19, nrow = Nr*Nc)
RI_1D_x_rot <- matrix(ncol = 19, nrow = Nr*Nc)
RI_1D_y_rot <- matrix(ncol = 19, nrow = Nr*Nc)
RI_12_rot <- matrix(ncol = 19, nrow = Nr*Nc)
nr <- nrow(pt0)
swapMean <- (nr-1)*(nr-2)/4
for (k in 1:19) {
  ang <- (k-1)*10/180*pi
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
pt2_sj <- pt0 + cbind(runif(Nr*Nc,min = -1,max = 1)*0, runif(Nr*Nc,min = -1,max = 1)*3) #jitter
pt2_sj <- pt2_sj %*% M_rot_30 

RI_12_rot_jitter <- matrix(ncol = 19, nrow = Nr*Nc)

nr <- nrow(pt0)
swapMean <- (nr-1)*(nr-2)/4
for (k in 1:19) {
  ang <- (k-1)*10/180*pi
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


# -- PLOT
df <- data.frame(cbind(seq(1,Nr*Nc), RI_12_rot))
df <- data.frame(cbind(seq(1,Nr*Nc), RI_2D_rot))
df <- data.frame(cbind(seq(1,Nr*Nc), RI_1D_x_rot))
# df <- data.frame(cbind(seq(1,Nr*Nc), RI_1D_y_rot)) 
df <- data.frame(cbind(seq(1,Nr*Nc), RI_12_rot_jitter))

colnames(df) <- c('neu', seq(0,180,by = 10))
dfm <- melt(df, id = 'neu')

dev.new()
# pdf(file = "synth index rotation.pdf", width = 10, height = 4,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  stat_summary(fun.data = mean_se,
               geom = "errorbar", position = position_nudge(x = 0, y = 0),
               colour = 'black', size = 1, width = .2 ) +
  stat_summary(fun = mean, geom = "point",  colour = 'black', size = 1, position = position_nudge(x = 0, y = 0)) +
  coord_cartesian(ylim = c(-.5, 1)) +
  theme_minimal_grid() +
  labs(title = "grid rotation 12 ", x='Angle', y='RI') 

dev.off()


