####################################################
### Reproduce the Figure 7 in TF
####################################################

source("~/Documents/UW/course work/17Spring/STAT 572/my codes/prelim_extra_functions.R")
library(lars) # for locally adaptive regression splines

### Reproduce the hills data example in the left panel
n <- 128; k <- 3; tf.index <- 200; lars.index <- 183

#generate the "hill" data
x = 1:n/n
u = smoothwiggly(n,x)

set.seed(0)
y = u + rnorm(n,sd=0.005)
miny = min(y); maxy = max(y)
y = (y-miny)/(maxy-miny)*10
u = (u-miny)/(maxy-miny)*10

# fit trend filtering
tf.X <- generate.TF.H.design(n, k)
tf.sol <- non.stan.lasso(tf.X, y = y, k = k, tol = 1e-8)
tf.lambda <- tf.sol$lambda
tf.fit <- tf.sol$beta[, tf.index]

# fit locally adaptive regression splines
lars.X <- generate.LARS.design(n, k)
lars.sol <- non.stan.lasso(lars.X, y = y, k = k, tol = 1e-8)
lars.lambda <- lars.sol$lambda
lars.fit <- lars.sol$beta[, lars.index]

# plot the result
pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/Figure7.pdf", w = 10, h = 5)
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2) + 1)
#plot the true function
plot(x, y, type = "p", col = "grey60", cex = 0.8, main = "Hills example")
lines(x, tf.fit, col = "red")
lines(x, lars.fit, col = "blue", lty = 2)
legend("topleft", c("Trend Filtering", "Locally Adaptive Splines"), col = c("red", "blue"), lty = c(1, 2), cex = 0.7)
### End of reproducing the hills data example in the left panel


### Reproduce the Doppler data example in the right panel
n <- 1000; k <- 3; tf.index <- 3500; lars.index <- 3505;

#generate the simulation data for the Doppler function
set.seed(0)
doppler.signal <- generate.doppler(n, noise = 0.1)
x <- seq(n)/n
y <- doppler.signal$y
#plot(doppler.signal$x, doppler.signal$means, type = "l")
#points(doppler.signal$x, doppler.signal$y, cex=0.5, col = "grey60")

# fit trend filtering
tf.X <- generate.TF.H.design(n, k)
tf.sol <- non.stan.lasso(tf.X, y = y, k = k, tol = 1e-12)
tf.lambda <- tf.sol$lambda

D.cubic <- difference.matrix(1000, 3)
temp <- proc.time()
my.fit <- my.trendfilter(y, D = D.cubic, maxsteps = 5000)
proc.time() - temp

# fit locally adaptive regression splines
lars.X <- generate.LARS.design(n, k)
lars.sol <- non.stan.lasso(lars.X, y = y, k = k, tol = 1e-12)
lars.lambda <- lars.sol$lambda


tf.fit <- tf.sol$beta[, tf.index]
lars.fit <- lars.sol$beta[, lars.index]
# plot the result
#pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/Figure7.pdf", w = 10, h = 10)
#par(mfrow = c(1, 1), mar = c(2, 2, 2, 2) + 1)
#plot the true function
plot(x, y, type = "p", col = "grey60", cex = 0.8, main = "Doppler example")
lines(x, tf.fit, col = "red")
lines(x, lars.fit, col = "blue", lty = 2)
legend("topleft", c("Trend Filtering", "Locally Adaptive Splines"), col = c("red", "blue"), lty = c(1, 2), cex = 0.7)

dev.off()

### End of reproducing the Doppler data example in the right panel