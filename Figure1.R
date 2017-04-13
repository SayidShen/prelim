####################################################
### Reproduce the Figure 1 in TF
####################################################

source("prelim_extra_functions.R")
library(genlasso)

n = 100

# Compute the piecewise-constant case in Figure 1 of TF paper
set.seed(10)
const.signal <- gen.splines.signals(n, 1, n.knots = 4, noise = 0.02)
const.pen.mat <- difference.matrix(n, 1)
const.tf <- genlasso(const.signal$y, diag(rep(1, n)), const.pen.mat)
# Done computing the piecewise-constant case in Figure 1 of TF paper

# Compute the piecewise-linear case in Figure 1 of TF paper
set.seed(40)
linear.signal <- gen.splines.signals(n, 2, n.knots = 3, noise = 0.01)
linear.pen.mat <- difference.matrix(n, 2)
linear.tf <- genlasso(linear.signal$y, diag(rep(1, n)), linear.pen.mat)
# Done computing the piecewise-constant case in Figure 1 of TF paper

# Compute the piecewise-quadratic case in Figure 1 of TF paper
set.seed(100)
quad.signal <- gen.splines.signals(n, 3, n.knots = 1, noise = 0.02, coefficients = c(0, 0, -3, 5))
quad.pen.mat <- difference.matrix(n, 3)
quad.tf <- genlasso(quad.signal$y, diag(rep(1, n)), quad.pen.mat)
# Done computing the piecewise-constant case in Figure 1 of TF paper

# Reproduce the Figure 1 of TF paper
pdf("../my presentation/1st presentation/figures/Figure1.pdf", w = 10, h = 4)
par(mfrow = c(1, 3), mar = c(2, 2, 2, 2))
#plot the piecewise-constant case
plot(const.signal$x, const.signal$y, col = "grey30", cex = 0.8, xlab = "x", ylab = "y")
lines(const.signal$x, const.tf$fit[, 14], type = "l", lwd = 1.5, col = "red")
#plot the piecewise-linear case
plot(linear.signal$x, linear.signal$y, col = "grey30", cex = 0.8, xlab = "x", ylab = "y")
lines(linear.signal$x, linear.tf$fit[, 45], type = "l", lwd = 1.5, col = "blue")
#plot the piecewise-quadratic case
plot(quad.signal$x, quad.signal$y, col = "grey30", cex = 0.8, xlab = "x", ylab = "y")
lines(quad.signal$x, quad.tf$fit[, 40], type = "l", lwd = 1.5, col = "green")

dev.off()


