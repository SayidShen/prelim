####################################################
### Reproduce the Figure 2 in TF
####################################################

source("~/Documents/UW/course work/17Spring/STAT 572/my codes/prelim_extra_functions.R")
#source("~/Documents/UW/course work/17Spring/STAT 572/my codes/lasso_solvers.R")

n <- 100; 
# Reproduce the left panel of Figure2
#generate piecewise linear signal

set.seed(5)
linear.signal <- gen.signal(seq(n)/n, M = 1, knots = c(0.2, 0.5, 0.8), noise = 0.1)
#fit linear trend filtering
D.linear <- difference.matrix(n, 2)
linear.tf <- my.trendfilter(linear.signal$y, D.linear)
linear.tf <- linear.tf$primal.fit[, 22]
linear.tf.pen <- D.linear %*% linear.tf
linear.tf.knots <- linear.signal$x[which(abs(linear.tf.pen) > 1e-8) + 1]
#fit linear regression spline
linear.spline <- splines(linear.signal$x, linear.signal$y, M = 1, knots = linear.tf.knots)

# Reproduce the right panel of Figure2
#generate piecewise cubic signal

set.seed(5)

cubic.signal <- gen.signal(seq(n)/n, M = 3, knots = c(0.2, 0.5), noise = 0.5)
D.cubic <- difference.matrix(n, 4)
cubic.tf <- my.trendfilter(cubic.signal$y, D.cubic)
cubic.tf <- cubic.tf$primal.fit[, 30]
cubic.tf.pen <- D.cubic %*% cubic.tf
cubic.tf.knots <- cubic.signal$x[which(abs(cubic.tf.pen) > 1e-8) + 1]
#fit cubic regression spline
cubic.spline <- splines(cubic.signal$x, cubic.signal$y, M = 3, knots = cubic.tf.knots)



pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/Figure2.pdf", w = 10, h = 4)
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
#plot the left panel
plot(linear.signal$x, linear.signal$y, col = "grey30", cex = 0.8, xlab = "x", ylab = "y", main = "Piecewise Linear Signal")
lines(linear.signal$x, linear.signal$means, col = "black")
abline(v = linear.tf.knots, lty = 2)
lines(linear.signal$x, linear.tf, col = "red")
#lines(linear.signal$x, linear.tf2, col = "green")
lines(linear.signal$x, linear.spline, col = "blue", lty = 2)
legend("topleft", c("Trend Filtering", "Regression Splines"), lty = c(1, 2), col = c("red", "blue"), cex = 0.5)
#plot the right panel
plot(cubic.signal$x, cubic.signal$y, col = "grey30", cex = 0.8, xlab = "x", ylab = "y", main = "Piecewise Cubic Signal")
lines(cubic.signal$x, cubic.signal$means, col = "black")
abline(v = cubic.tf.knots, lty = 2)
#lines(cubic.signal$x, ryan.tf$beta[, 80], col = "red")
lines(cubic.signal$x, cubic.tf, col = "red")
lines(cubic.signal$x, cubic.spline, col = "blue", lty = 2)
#legend("bottomleft", c("Trend Filtering", "Regression Splines"), lty = c(1, 2), col = c("red", "blue"), cex = 0.8)

par(mfrow = c(1,1))
dev.off()
