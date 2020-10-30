# Kendall tau distance, https://en.wikipedia.org/wiki/Kendall_tau_distance
# No. swaps = No. inversion

# Kendall rank corr coeff, https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient
# ordinary generating func, https://math.stackexchange.com/questions/245018/bubble-sorting-question

# on a grid ----------------------------------------------------------------------------------------------------

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

# getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
# col_66 <- (rev(getPalette(25))[1:(Nr*Nc/3)])
# getPalette <- colorRampPalette(brewer.pal(9, "Reds"))
# col_66 <- c(col_66, rev(getPalette(25))[1:(Nr*Nc/3)])
# getPalette <- colorRampPalette(brewer.pal(9, "Purples"))
# col_66 <- c(col_66, rev(getPalette(25))[1:(Nr*Nc/3)])

# col_66 <- c('#deebf7', '#d8e5f3', '#d2deee', '#ccd8ea', '#c5d2e5', '#bfcce1', '#b9c5dc', '#b3bfd8', '#adb9d4', '#a7b3cf', '#a1adcb', '#9ba7c7', '#95a1c2', '#8f9bbe', '#8995ba', '#8390b6', '#7d8ab1', '#7784ad', '#727ea9', '#6c79a5', '#6673a0', '#606e9c', '#5a6898', '#546394', '#4e5d90', '#48588c', '#415387', '#3b4e83', '#34497f', '#2d437b', '#263e77', '#1e3a73', '#15356f', '#6c0b13', '#711518', '#761d1d', '#7b2422', '#802b27', '#85322c', '#8a3831', '#8f3e37', '#93453c', '#984b42', '#9d5147', '#a2574d', '#a65e53', '#ab6459', '#af6a5f', '#b47165', '#b9776b', '#bd7d71', '#c28477', '#c68a7d', '#cb9083', '#cf978a', '#d39d90', '#d8a496', '#dcaa9d', '#e0b1a3', '#e5b8aa', '#e9beb0', '#edc5b7', '#f1ccbe', '#f6d2c4', '#fad9cb', '#fee0d2')
# 
# col_66 <- lacroix_palette("Pamplemousse", n = 66, type = "continuous")


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
# pdf(file = "index rotation.pdf", width = 10, height = 4,pointsize=12,family="Helvetica", useDingbats = F)

ggplot(dfm, aes(factor(variable), value) ) + 
  # geom_boxplot() +
  # geom_point() +
  # geom_jitter(colour = 'black', size = 4, height = 0, width = 0.1) +
  stat_summary(fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               geom = "linerange", position = position_nudge(x = 0, y = 0),
               colour = 'black', size = 1 ) +
  stat_summary(fun.min = median, geom = "errorbar",  colour = 'black', size = 1, width = .2, position = position_nudge(x = 0, y = 0) ) +
  coord_cartesian(ylim = c(-.5, 1)) +
  theme_minimal_grid() +
  labs(title = "grid rotation 12 ", x='Angle', y='RI') 

dev.off()


# stats\ ----------------------------------------------------------------------------------------------------------

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



# 1D and 2D RI range test ---------------------------------------------------------------------------------------------------------

# original
N <- 4
pt <- cumsum(seq(0, by =0.003, length.out =N)) + seq(0,by=1, length.out=N) 
# pt <- rev(pt)

sum(duplicated(dist(pt)))
N*(N-1)/2

# all
RI_pt <- c()
RI_pt_avg <- c()
newOrder <- gtools::permutations(N,N) #N <=8
for (j in 1:nrow(newOrder)) {
  ptn <- pt[newOrder[j,]]
  RI_pt <- c(RI_pt, RI(pt, ptn))
  RI_pt_avg <- c(RI_pt_avg, mean(RI(pt, ptn)) )
}

dev.new()
hist(RI_pt,50)
dev.new()
hist(RI_pt_avg,50)

# random sample
RI_pt <- c()
RI_pt_avg <- c()
for (j in 1:500) {
  ptn <- pt[sample(N)]
  RI_pt <- c(RI_pt, RI(pt, ptn))
  RI_pt_avg <- c(RI_pt_avg, mean(RI(pt, ptn)) )
}


# PLOT scatter
df <- as.data.frame(RI_pt)
df <- as.data.frame(RI_pt_avg)


# - 2D

# original
Nr <- 4
Nc <- 1
N <- Nc*Nr
y <- seq(-Nr/2+0.5, Nr/2-0.5, length.out = Nr)
x <- seq(-Nc/2+0.5, Nc/2-0.5, length.out = Nc)
pt0 <- expand.grid(x, y) + cbind(runif(Nr*Nc,min = -1,max = 1)*0.1, runif(Nr*Nc,min = -1,max = 1)*0.1)
pt0 <- as.matrix(pt0)

sum(duplicated(dist(pt0)))
N*(N-1)/2

pt1 <- pt0
# for (j in 1:Nr) {
#   pt1[seq(1,Nc)+(j-1)*Nc,1] <- pt1[seq(1,Nc)+(j-1)*Nc,1]*0.05*j
# }

# -- all
RI_pt <- c()
RI_pt_avg <- c()
newOrder <- gtools::permutations(N,N) #N <=8
for (j in 1:nrow(newOrder)) {
  pt2 <- pt1[newOrder[j,],]
  RI_pt <- c(RI_pt, RI(pt0, pt2))
  RI_pt_avg <- c(RI_pt_avg, mean(RI(pt0, pt2)) )
}

# # -- random sample
# RI_pt <- c()
# RI_pt_avg <- c()
# for (j in 1:1000) {
#   pt2 <- pt1[sample(N),]
#   RI_pt <- c(RI_pt, RI(pt0, pt2))
#   RI_pt_avg <- c(RI_pt_avg, mean(RI(pt0, pt2)) )
# }

dev.new()
hist(RI_pt,breaks = seq(-1,1,length.out = 51), main = paste(Nr,'x',Nc,sep=''))
dev.new()
hist(RI_pt_avg,breaks = seq(-1,1,length.out = 51),main = paste(Nr,'x',Nc,sep=''))


# - 2D random stats

# original
Nr <- 6
Nc <- 11
N <- Nc*Nr

pt0 <- cbind(runif(N,min = -1,max = 1), runif(N,min = -1,max = 1)) %>% as.matrix()
sum(duplicated(dist(pt0)))
# N*(N-1)/2
# factorial(N-1)

# -- random sample
RI_pt <- c()
RI_pt_avg <- c()
for (j in 1:10000) {
  pt1 <- cbind(runif(N,min = -1,max = 1), runif(N,min = -1,max = 1)) %>% as.matrix()
  RI_pt <- c(RI_pt, RI(pt0, pt1))
  RI_pt_avg <- c(RI_pt_avg, mean(RI(pt0, pt1)) )
}

dev.new()
hist(RI_pt,breaks = seq(-.5,.5,length.out = 51), main = paste(Nr,'x',Nc,sep=''))
dev.new()
hist(RI_pt_avg,breaks = seq(-.1,.1,length.out = 51),main = paste(Nr,'x',Nc,sep=''))

# distr
sigma2 <- 2*(2*N+5)/(9*N*(N-1))
x <- seq(-1,1,length.out = 200)
y <- 1/sqrt(2*pi*sigma2)*exp(-x^2/2/sigma2)


cc <- hist(RI_pt,breaks = seq(-1,1,length.out = 51), plot = F)
dev.new()
plot(x,y/max(y),type='l')
points(cc$mids, cc$counts/max(cc$counts), cex = 1, pch = 16)


# - check var ratio
varm <- matrix(ncol = 4, nrow = 0)
for (N in seq(10,100,by = 10)) {
  pt0 <- cbind(runif(N,min = -1,max = 1)) %>% as.matrix()
  RI_pt <- c()
  RI_pt_avg <- c()
  for (j in 1:500) {
    pt1 <- cbind(runif(N,min = -1,max = 1)) %>% as.matrix()
    RI_pt <- c(RI_pt, RI(pt0, pt1))
    RI_pt_avg <- c(RI_pt_avg, mean(RI(pt0, pt1)) )
  }
  varm <- rbind(varm, c(N,var(RI_pt), var(RI_pt_avg), var(RI_pt)/var(RI_pt_avg)))
  
}



