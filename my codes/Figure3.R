####################################################
### Reproduce the Figure 3 in TF
####################################################

source("~/Documents/UW/course work/17Spring/STAT 572/my codes/prelim_extra_functions.R")

n <- 100

#generate piecewise cubic signal
set.seed(5)
#cubic.signal <- gen.splines.signals(n, M = 4, n.knots = 5, knots = c(0.1, 0.3, 0.4, 0.7, 0.9), noise = 0.02)
cubic.signal <- gen.signal(seq(n)/n, M = 3, knots = c(0.2, 0.5), noise = 0.5)
#plot(cubic.signal$x, cubic.signal$y)
#fit cubic trend filtering
D.cubic <- difference.matrix(n, 4)
cubic.tf <- my.trendfilter(cubic.signal$y, D.cubic)
cubic.tf <- cubic.tf$primal.fit[, 30]
cubic.tf.pen <- D.cubic %*% cubic.tf
cubic.tf.knots <- cubic.signal$x[which(abs(cubic.tf.pen) > 1e-8) + 3]

pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/Figure3.pdf", w = 10, h = 10)
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
#plot the estimate
plot(cubic.signal$x, cubic.tf, col = "red", cex = 1.4, main = "Estimate")
abline(v = cubic.tf.knots, lty = 2)
#plot the 1st derivative
der1.mat <- difference.matrix(n, 1)
der1 <- der1.mat %*% cubic.tf
plot(cubic.signal$x[-1], der1, col = "red", cex = 1.4, main = "1st Derivative")
abline(v = cubic.tf.knots, lty = 2)
#plot the 2nd derivative
der2.mat <- difference.matrix(n, 2)
der2 <- der2.mat %*% cubic.tf
plot(cubic.signal$x[-c(1,2)], der2, col = "red", cex = 1.4, main = "2nd Derivative")
abline(v = cubic.tf.knots, lty = 2)
#plot the 3rd derivative
der3.mat <- difference.matrix(n, 3)
der3 <- der3.mat %*% cubic.tf
plot(cubic.signal$x[-c(1,2,3)], der3, col = "red", cex = 1.4, main = "3rd Derivative")
abline(v = cubic.tf.knots, lty = 2)

par(mfrow = c(1,1))
dev.off()
