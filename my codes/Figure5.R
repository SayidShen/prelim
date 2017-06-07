####################################################
### Reproduce the Figure 5 in TF
####################################################

source("~/Documents/UW/course work/17Spring/STAT 572/my codes/prelim_extra_functions.R")

n <- 1000

#generate the simulation data for the Doppler function
set.seed(0)
doppler.signal <- generate.doppler(n, noise = 0.1)
plot(doppler.signal$x, doppler.signal$means, type = "l")
points(doppler.signal$x, doppler.signal$y, cex=0.5, col = "grey60")
#fit of trend filtering
k <- 3; index <- 2599
D.cubic <- difference.matrix(n, k + 1)
tf.fit <- my.trendfilter(doppler.signal$y, D.cubic, verbose = F, maxsteps = 8000)
tf <- tf.fit$primal.fit[, index]
tf.df <- tf.fit$df[index] 
#fit smoothing spline
ss.fit1 <- smooth.spline(seq(n)/n, doppler.signal$y, df = 50)$y
ss.fit2 <- smooth.spline(seq(n)/n, doppler.signal$y, df = 90)$y




#plot figure 5
pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/Figure5.pdf", w = 10, h = 10)
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
#plot the true function
plot(doppler.signal$x, doppler.signal$y, col = "grey60", cex = 0.8, main = "True function" )
lines(doppler.signal$x, doppler.signal$means, lwd = 1.2)

#plot the fit of trend filtering
plot(doppler.signal$x, doppler.signal$y, col = "grey60", cex = 0.8, main = paste("Trend Filtering, df=", tf.df, sep = ""))
lines(doppler.signal$x, tf, lwd = 1.2, col = "red")

#plot the fit of smoothing splines with the same dof
plot(doppler.signal$x, doppler.signal$y, col = "grey60", cex = 0.8, main = paste("Smoothing Splines, df=", 50, sep = ""))
lines(doppler.signal$x, ss.fit1, col = "blue", lwd = 1.2)

#plot the fit of smoothing splines with different dof
plot(doppler.signal$x, doppler.signal$y, col = "grey60", cex = 0.8, main = paste("Smoothing Splines, df=", 90, sep = ""))
lines(doppler.signal$x, ss.fit2, col = "blue", lwd = 1.2)

par(mfrow = c(1,1))
dev.off()