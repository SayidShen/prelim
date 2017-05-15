####################################################
### Reproduce the Figure 1 in TF
####################################################

source("~/Documents/UW/course work/17Spring/STAT 572/my codes/prelim_extra_functions.R")


n <- 100

# Compute the piecewise-constant case in Figure 1 of TF paper
set.seed(10)
const.signal <- gen.splines.signals(n, 1, n.knots = 3, knots = c(0.2, 0.4, 0.7), noise = 0.05)
D.const <- difference.matrix(n, 1)
const.tf <- my.trendfilter(as.vector(const.signal$y), D.const)$primal.fit[, 10]
# Done computing the piecewise-constant case in Figure 1 of TF paper

# Compute the piecewise-linear case in Figure 1 of TF paper
set.seed(40)
linear.signal <- gen.signal(seq(n)/n, 1, knots = c(0.2, 0.6, 0.8), noise = 0.05)
D.linear <- difference.matrix(n, 2)
linear.tf <- my.trendfilter(linear.signal$y, D.linear)$primal.fit[, 60]
# Done computing the piecewise-constant case in Figure 1 of TF paper

# Compute the piecewise-quadratic case in Figure 1 of TF paper
set.seed(100)
quad.signal <- gen.signal(seq(n)/n, 2, knots = c(0.4, 0.8), noise = 0.05)
D.quad <- difference.matrix(n, 3)
quad.tf <- my.trendfilter(quad.signal$y, D.quad)$primal.fit[, 40]
# Done computing the piecewise-constant case in Figure 1 of TF paper

# Reproduce the Figure 1 of TF paper
pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/Figure1.pdf", w = 10, h = 4)
par(mfrow = c(1, 3), mar = c(2, 2, 2, 2))
#plot the piecewise-constant case
plot(const.signal$x, const.signal$y, col = "grey30", cex = 0.8, xlab = "x", ylab = "y")
#lines(const.signal$x, const.tf$fit[, 14], type = "l", lwd = 1.5, col = "red")
lines(const.signal$x, const.tf, lwd = 1.5, col = "red")
#plot the piecewise-linear case
plot(linear.signal$x, linear.signal$y, col = "grey30", cex = 0.8, xlab = "x", ylab = "y")
#lines(linear.signal$x, linear.tf$fit[, 45], type = "l", lwd = 1.5, col = "blue")
lines(linear.signal$x, linear.tf, lwd = 1.5, col = "blue")
#plot the piecewise-quadratic case
plot(quad.signal$x, quad.signal$y, col = "grey30", cex = 0.8, xlab = "x", ylab = "y")
#lines(quad.signal$x, quad.tf$fit[, 40], type = "l", lwd = 1.5, col = "green")
lines(quad.signal$x, quad.tf, lwd = 1.5, col = "green") 

par(mfrow = c(1,1))

dev.off()


