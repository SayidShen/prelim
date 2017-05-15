####################################################
### Reproduce the Figure 4 in TF
####################################################

source("~/Documents/UW/course work/17Spring/STAT 572/my codes/prelim_extra_functions.R")

n <- 128

#generate the "hill" data
x = 1:n/n
u = smoothwiggly(n,x)

set.seed(1)
y = u + rnorm(n,sd=0.005)
miny = min(y); maxy = max(y)
y = (y-miny)/(maxy-miny)*10
u = (u-miny)/(maxy-miny)*10
#plot(cubic.signal$x, cubic.signal$means, type = "l")

#fit cubic trend filtering
k <- 3; index <- 208;
D.cubic <- difference.matrix(n, k + 1)
tf.fit <- my.trendfilter(y, D.cubic)
tf <- tf.fit$primal.fit[, index]

tf.dof <- tf.fit$df[index]
#fit smoothing spline
ss.fit1 <- smooth.spline(x, y, df = tf.dof)$y
ss.fit2 <- smooth.spline(x, y, df = tf.dof + 10)$y

#plot figure 4
pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/ssvslars.pdf", w = 10, h = 10)
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
#plot the true function
plot(x, u, type = "l", lwd = 1.5, main = "True Function")
points(x, y, col = "grey60", cex = 0.8)
#plot the fit of Trend Filtering
#plot(x, tf, type = "l", col = "red", lwd = 1.5, main = paste("Trend Filtering, df=", round(tf.dof), sep=""))
plot(x, tf, type = "l", col = "red", lwd = 1.5, main = "Locally Adaptive Splines")
points(x, y, col = "grey60", cex = 0.8)
#plot the fit of Smoothing Spline with same dof
plot(x, ss.fit1, type = "l", col = "blue", lwd = 1.5, main = paste("Smoothing Splines, df=", round(tf.dof), sep=""))
points(x, y, col = "grey60", cex = 0.8)
#Plot the fit of smoothing spline with larger dof
plot(x, ss.fit2, type = "l", col = "blue", lwd = 1.5, main = paste("Smoothing Splines, df=", round(tf.dof + 10), sep=""))
points(x, y, col = "grey60", cex = 0.8)

par(mfrow = c(1,1))
dev.off()

