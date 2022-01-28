library(rgl)
library(latex2exp)

#Let s be the seasonality parameter.
s <- seq(1, 2, length.out = 100)

#Let a be the long-memory parameter.
a <- seq(0, 1/2, length.out = 100)

#This function computes the correlation of the components of the asymptotic vector.
zcoor <-function(s,a){
  (((((1-(4*a*log(s)))*a*(4*a+2)*(s^(-1)))/(4*(pi^2))) 
      - ((8*a*s^3*log(s))/(((8*pi*(pi^2-2))/9)^2)))/(sqrt(((((1-(4*a*log(s)))^2)/(4*pi^2)) + ((8*s^4*(log(s))^2)/(((8*pi*(pi^2-2))/9)^2)))*(((a^2*((4*a+2)^2)*(s^(-2)))/(4*pi^2))+ ((8*a^2*s^2)/(((8*pi*(pi^2-2))/9)^2))))))  
}


rho <- outer(s,a, zcoor)

#This jet.colors function interpolates "red" and "yellow" colors and creates a new color palette.
jet.colors <- colorRampPalette( c("red", "yellow") ) 
pal <- jet.colors(100)
col.ind <- cut(rho,100) 


um <- matrix(c(-0.5988361,  0.8008499, 0.005933772 ,    0,
               -0.1430069, -0.1142176, 0.983109236,    0,
               0.7880003,  0.5878722, 0.182924747,    0,
               0.0000000,  0.0000000,  0.00000000,    1), 4,4,byrow =TRUE)
view3d(userMatrix = um)

######## Figure 3(b) - This figure gives the plot of the asymptotic correlation of hat{s_{0}} and hat{\alpha} based on the Meyer father wavelet.
persp3d(s,a,rho,col=pal[col.ind],xlab=" ",box = FALSE, ylab= " ",zlab=" ")

#These commands add axes labels of {s_{0}}, {\alpha} and {\rho} to their corresponding axes.
mtext3d(TeX("s_0"), "x-+", line = 3)
mtext3d(TeX('$\\alpha$'), "y+-", line = 3)
mtext3d(TeX('$\\rho$'), "z+-", line = 3)

#This gives Figure 3(b) in the PNG format.
rgl.snapshot("Fig3b.png")




