if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MassSpecWavelet")

pck <-
  c("polynom",
    "orthopolynom",
    "pracma",
    "MassSpecWavelet",
    "qboxplot",
    "latex2exp",
    "RcppRoll",
    "ggplot2",
    "dplyr",
    "ggpubr"
  )

lapply(pck, library, character.only = TRUE)

set.seed(754321)

M <- 10
k <- 100
fN <- 7

#Timing code
ptm <- proc.time()

delta1_mat1 <- matrix(ncol = fN, nrow = k*M)
delta2_mat1 <- matrix(ncol = fN - 1, nrow = k*M)



c2 <- 6.283185
c3 <- 31.41593

# Calculating constants c2 and c3
# library(wmtsa)
# w <- seq(-100, 100, length = 10000)
# c2 <-  Re(sum(diff(w) * (
#     wavCWTFilters(wavelet = "gaussian2", frequency = w)[-1]
#   ) ^ 2))
# c3 <- 31.41593
#   Re(2 * sum(diff(w) * ((
#     w * wavCWTFilters(wavelet = "gaussian2", frequency = w)
#   )[-1]) ^ 2))

for (m in 1:M){
  print(m)
  #Values used for simulations
  Ts <- 10000000
  kt <- 100
  k <- 100
  n0 <- 100
  
  
  #Parameters of Gegenbauer time series
  u1 <- 0.3
  d1 <- 0.1
  alpha <- d1
  s0 <- acos(u1)
  
  #Multiplier to normalize h(0)=1
  A <- sqrt((2 * pi * 4 ^ (2 * alpha) * (sin(s0 / 2)) ^ (4 * alpha)) / (s0^(4 * alpha)))
  
  #Gegenbauer polynomial
  gcoef <-
    function(i) {
      polynomial.values(gegenbauer.polynomials(i - 1, d1), u1)[[i]]
    }
  gcoef <- Vectorize(gcoef)
  coef <- seq(n0 + 1)
  coef <- gcoef(coef)
  coef1 <- rev(coef)
  n <- Ts / kt + 200 + n0
  
  #Creating epsilon:
  epsil0 <- rnorm(n, mean = 0, sd = 1)
  
  #Simulating Gengenbauer random process
  ts1 <- seq(Ts / kt + 200)

  ts1 <- A * roll_sum(epsil0, weights = coef1, normalize = FALSE,align = "left")
  
  #Range of scales and the number of scales
  ns <- 11
  nsc <- 20
  
  #Calculating the Continuous wavelet transform (CWT) of Gengenbauer random process
  ts1.cwt <- cwt(ts1, scales = seq(1,11,1), wavelet = "mexh")
  a1 <- as.matrix(ts1.cwt[101:(Ts / kt + 100)])
  #a1-a2-to see the difference in cwt of two packages
  ################
  # fN <- size(ts1.cwt)[2]
  # fN
  #attr(ts2.cwt, "scale")
  ################

  #aj <- attr(ts1.cwt, "scale")
  aj <- seq(1,11,1)
  squ_aj <- (aj) ^ (-2)
  
  #Calculating aj^(-2)- aj+1^(-2)
  denominator <- squ_aj[1:(fN - 1)] - squ_aj[2:(fN)]
  
  ###################
  delta1_mat <- matrix(ncol = fN, nrow = k)
  
  delta2j <- rep(0, fN)
  
  delta2_mat <- matrix(ncol = fN - 1, nrow = k)
  ####################
  
  ts1cwtmatr <- matrix(0, Ts - 100, 11)
  #Creating matrix with wavelet coefficients
  ts10 <- seq(Ts / kt + 100)
  for (N in 1:k) {
    epsil <- rnorm(Ts / kt + 100 + n0, mean = 0, sd = 1)
    ts10 <- A * roll_sum(epsil, weights = coef1, normalize = FALSE)
    ts10.cwt <-
      as.matrix(cwt(ts10, scales = seq(1,11,1), wavelet = "mexh"))
    ts1cwtmatr[1:(Ts / kt - 100), ] <- ts10.cwt[101:(Ts/kt), ]
    ts1 <- ts10
    
    for (j in 1:(kt - 1)) {
      epsil <-
        c(epsil[(length(epsil) - n0 + 1):length(epsil)], rnorm(Ts /kt, mean =0, sd = 1))
      ts10 <- A * roll_sum(epsil, weights = coef1, normalize = FALSE)
      ts1 <- c(ts1[(length(ts1) - 200 + 1):length(ts1)], ts10[1:(Ts / kt)])
      ts10.cwt <-
        as.matrix(cwt(ts1, scales = seq(1,11,1), wavelet = "mexh"))
      ts1cwtmatr[(Ts / kt - 100 + 1 + (j - 1) * Ts / kt):(Ts / kt - 100 +
                                                            j * Ts / kt), ] <- ts10.cwt[101:(Ts / kt + 100), ]
    }
    
    ts1cwtmatr <- (ts1cwtmatr * ts1cwtmatr)
    
    #Calculating deltaj^(2)- deltaj+1^(2)
    for (j in 1:(fN)) {
      delta2j[j] <-
        round(mean(ts1cwtmatr[1:min(Ts - 100, 200 * (j < 3) + (aj[j]) ^ 9), j]), 5)
    }
    delta1_mat[N, ] <- delta2j
    numerator <- delta2j[1:(fN - 1)] - delta2j[2:fN]
    
    #Calculating Delta(deltaj.^(2))
    delat_result <- numerator / denominator
    delta2_mat[N, ] <- delat_result
    
  delta1_mat1[(1+k*(m-1)):(m*k),] <- delta1_mat
  delta2_mat1[(1+k*(m-1)):(m*k),]<-delta2_mat
  }
  rm(list=setdiff(ls(), c("ptm","m","M", "k", "fN", "c2","c3", "delta1_mat1", "delta2_mat1", "s_mat1", "alpha_mat1")))
 
}

timetaken1 <- (proc.time() - ptm) / 60

save.image(file = "Longmemory_estimates1.RData")

# load("Longmemory_estimates1.RData")
######################################################################################################################
u1 <- 0.3
d1 <- 0.1
alpha <- d1
s0 <- acos(u1)
delta1_mat <- delta1_mat1
delta2_mat<-delta2_mat1

#Box Plots of bar(delta)j^(2),   bar(Delta(deltaj.^(2))),   s^hat_0j   and   alpha^hat_j
q1 <- 1

#Box Plot of bar(delta)j^(2) for j=1,2,3,4,5,6,7
qboxplot(
  data.frame(delta1_mat[, q1:(fN)]),
  probs = c(0.25, 0.5, 0.75),
  range = 0.0,
  col = "bisque",
  xaxt = "n",
  medcol = "bisque",
  main = TeX("Boxplot of  $\\bar{\\delta}_{j\\cdot}^{(2)}$")
)
m1 <- apply(delta1_mat[, q1:(fN)], 2, mean)
points((q1:(fN)), m1, col = "red", pch = 18)
axis(1, at = (q1:(fN)), labels = (q1:(fN)))
abline(h = c2 * s0 ^ (-4 * alpha),
       col = "red",
       lty = 3)

#####################################################################################################################
#Box Plot of bar(Delta(deltaj.^(2))) for j=1,2,3,4,5,6
qboxplot(
  data.frame(delta2_mat[,q1:(fN - 1)]),
  probs = c(0.25, 0.5, 0.75),
  range = 0.0,
  col = "bisque",
  xaxt = "n",
  medcol = "bisque",
  main = TeX("Boxplot of  $\\Delta \\bar{\\delta}_{j\\cdot}^{(2)}$")
)
m1 <- apply(delta2_mat[, q1:(fN - 1)], 2, mean)
points((q1:(fN - 1)), m1, col = "red", pch = 18)
axis(1, at = (q1:(fN - 1)), labels = (q1:(fN - 1)))
abline(h = c3 * alpha * s0 ^ (-4 * alpha - 2),
       col = "red",
       lty = 3)

####################################################################################################################

Ts <- 10000000

shapiro.test(sqrt(Ts - 100)*((delta1_mat1[, (fN)])/c2-mean((delta1_mat1[, (fN)])/c2)))
shapiro.test(sqrt(Ts - 100)*((delta2_mat1[, (fN-1)])/c3-mean((delta2_mat1[, (fN-1)])/c3)))

####### Figure 4(a) - This figure gives the Q-Q plot of S1
ggqqplot(sqrt(Ts - 100)*((delta1_mat1[, (fN)])/c2-mean((delta1_mat1[, (fN)])/c2)))

####### Figure 4(b) - This figure gives the Q-Q plot of S2
ggqqplot(sqrt(Ts - 100)*((delta2_mat1[, (fN-1)])/c3-mean((delta2_mat1[, (fN-1)])/c3)))

p <- cbind((delta1_mat1[, (fN-3)])/c2,(delta2_mat1[,(fN-1)])/c3)

sum(p[,2]>0)
sum(p[,2]<0.5*p[,1]^2)
sum(p[,1]>0)
sum(p[,1]<1)

p <- p[p[,2]>0 & p[,2]<1 & p[,2]<0.5*p[,1]^2,]


m10 <- mean(p[,1])
m20 <- mean(p[,2])

pnorm <- cbind( sqrt(Ts - 100)*(p[,1]-m10),
                sqrt(Ts - 100)*(p[,2]-m20))

sigma <- cov(pnorm)

#Computing the correlation matrix
cor(pnorm)

sigma.inv <- solve(sigma, matrix(c(1,0,0,1),2,2))


ellipse <- function(s,t) {u<-c(s,t)-c(0, 0); u %*% sigma.inv %*% u / 2}


######## Figure 4(c) - This figure gives the density ellipsoid of (S1,S2).
plot(pnorm, pch=20, xlab=TeX('$S_1$'),ylab=TeX('$S_2$'))
points(0,0, col = "red", lty = 2, pch = 18, cex = 3)

n <- 200
x <- seq(min(pnorm[,1]), max(pnorm[,1]),length.out=n)
y <-  seq(min(pnorm[,2]), max(pnorm[,2]),length.out=n)
z <- mapply(ellipse, as.vector(rep(x,n)), as.vector(outer(rep(0,n), y, `+`)))
contour(x,y,matrix(z,n,n), levels=(0:10), col = terrain.colors(11), add=TRUE)


q <- p[,2] / p[,1]
#To find s^hat0j and alpha^hatj; equation (5.5)
hat_s0j <-
  exp(0.5 * lambertWp(
    log(1 / p[,1]) / (2 * q)))
str(hat_s0j)

hat_alphaj <-
  q * (exp(lambertWp(log(1 / p[,1]) / (2 * q))))
str(hat_alphaj)


m30 <- mean(hat_s0j)
m40 <- mean(hat_alphaj)

pnorm1 <- cbind(sqrt(Ts - 100)*(hat_s0j-m30),
                sqrt(Ts - 100)*(hat_alphaj-m40))

sigma1 <- cov(pnorm1)
sigma.inv1 <- solve(sigma1, matrix(c(1,0,0,1),2,2))


ellipse <- function(s,t) {u<-c(s,t)-c(0, 0); u %*% sigma.inv1 %*% u / 2}


######## Figure 4(d) - This figure gives the density ellipsoid of sqrt{m_j}((hat{s_0,alpha}_j)_(s_0,alpha))
plot(pnorm1, pch=20, xlab=TeX('$\\sqrt{m_j}(\\hat{s}_{0,j}-s_0)$'),ylab=TeX('$\\sqrt{m_j}(\\hat{\\alpha}_j-\\alpha)$'),mgp=c(2.6,1,0))
points(0,0, col = "red", lty = 2, pch = 18, cex = 3)

n <- 200
x <- seq(min(pnorm1[,1]), max(pnorm1[,1]),length.out=n)
y <-  seq(min(pnorm1[,2]), max(pnorm1[,2]),length.out=n)
z <- mapply(ellipse, as.vector(rep(x,n)), as.vector(outer(rep(0,n), y, `+`)))
contour(x,y,matrix(z,n,n), levels=(0:10), col = terrain.colors(11), add=TRUE)



