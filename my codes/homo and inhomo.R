#########################################
##### Example of homogeneous smoothness and inhomogeneous smoothness
#########################################

n <- 100
# generate the homogeneous smooth data
x.h <- 1:n/n
u.h <- gen.splines.signals(n, 3, n.knots = 1, noise = 0.5)
u.h <- as.vector(u.h$means)


#generate the "hill" data
x = 1:n/n
u = smoothwiggly(n,x)

set.seed(0)
y = u + rnorm(n,sd=0.005)
miny = min(y); maxy = max(y)
y = (y-miny)/(maxy-miny)*10
u = (u-miny)/(maxy-miny)*10


pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/homo.pdf", w = 10, h = 6)
par(mfrow = c(1, 2), mar = c(2,2,2,2) + 1)
plot(x.h, u.h, type = "l", main = "Homogeneous Smoothness", col = "blue", lwd = 1.5)
plot(x, u, type = "l", main = "Inhomogeneous Smoothness", col = "blue", lwd = 1.5)
dev.off()

