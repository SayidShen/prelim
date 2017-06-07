#################################################
### Figures in the introductory presentation
#################################################

source("prelim_extra_functions.R")


# Generate a piecewise quadratic signal
n = 100
quad.signal <- gen.splines.signals(n, 3, n.knots = 2, knots = c(0.3, 0.7), noise = 0.05, coefficients = c(5, 1, -2, 5, -8))
knots <- quad.signal$knots
x <- quad.signal$x
# Done generating a piecewise quadratic signal

pdf("../my presentation/1st presentation/figures/splines.pdf", w = 10, h = 6)
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
# fit piecewise linear model
piecewise.linear.fit <- piecewise.polynomial(quad.signal, 1)
plot(x, quad.signal$means, type = "l", xlab = "x", ylab = "y", col = "blue", main = "Piecewise Linear Fit", lwd = 1.5)
points(x, quad.signal$y, col = "grey50", cex = 0.8)

piece.index <- piecewise.linear.fit$piece.index; fit <- piecewise.linear.fit$fit
for (j in 1 : length(knots))
  abline(v = knots[j], lty = 2)
for (j in 1 : (length(knots) + 1))
  lines(x[piece.index[[j]]], fit[[j]], col = "red", lwd = 1.5)
legend("bottomright", c("True Signal", "Piecewise Linear"), pch = 15, col = c("blue", "red"), cex = 0.78)

# fit linear splines
linear.splines <- splines(quad.signal, 1)
fit <- linear.splines$fit
plot(x, quad.signal$means, type = "l", xlab = "x", ylab = "y", col = "blue", main = "Linear Splines", lwd = 1.5)
points(x, quad.signal$y, col = "grey50", cex = 0.8)
lines(x, fit, col = "red", lwd = 1.5)
for (j in 1 : length(knots))
  abline(v = knots[j], lty = 2)
legend("bottomright", c("True Signal", "Linear Splines"), pch = 15, col = c("blue", "red"), cex = 0.8)

dev.off()

#fit smoothing splines
smooth.splines.fit <- smoothing.splines(quad.signal, df = 15)
fit <- smooth.splines.fit$fit

pdf("../my presentation/1st presentation/figures/smoothing_splines.pdf", w = 10, h = 6)
par(mfrow = c(1, 1), mar = c(2, 2, 2, 2))
plot(x, quad.signal$means, type = "l", xlab = "x", ylab = "y", col = "blue", main = "Cubic Smoothing Splines (df = 15)", lwd = 1.5)
points(x, quad.signal$y, col = "grey50", cex = 0.8)
lines(x, fit, col = "red", lwd = 1.5)
for (j in 1 : length(knots))
  abline(v = knots[j], lty = 2)
legend("bottomright", c("True Signal", "Cubic Smoothing Splines"), pch = 15, col = c("blue", "red"), cex = 0.8)
dev.off()




