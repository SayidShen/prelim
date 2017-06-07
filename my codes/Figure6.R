####################################################
### Reproduce the Figure 6 in TF
####################################################

source("~/Documents/UW/course work/17Spring/STAT 572/my codes/prelim_extra_functions.R")

### Reproduce the MSE result of the "hills" data

n <- 128; nsim <- 50; k <- 3; n.dfs <- 36;
df.seq <- seq(from = 5, to = 40)
hill.tf.MSE <<- hill.ss.MSE <- matrix(NA, nrow = n.dfs, ncol = nsim)
D.cubic <- difference.matrix(n, k + 1)

set.seed(1)
seeds <- sample(1:1000, nsim)

for (i in 1:nsim)
{
  set.seed(seeds[i])
  cat(i, fill = F)
  # generate the data
  x = 1:n/n
  u = smoothwiggly(n,x)
  y = u + rnorm(n,sd=0.004)
  miny = min(y); maxy = max(y)
  y = (y-miny)/(maxy-miny)*10
  u = (u-miny)/(maxy-miny)*10

  # fit trend filtering
  tf.fit <- my.trendfilter(y, D.cubic)
  #block <- floor(300 / n.lams)
  for (j in 1 : n.dfs)
  {
    index <- median(which(tf.fit$df == round(df.seq[j])))
    hill.tf.MSE[j, i] <- mse(tf.fit$primal.fit[, index], u)
    #hill.df[j, i] <- tf.fit$df[j * block]

    #fit smoothing splines
    ss.fit <- smooth.spline(x, y, df = df.seq[j])$y
    hill.ss.MSE[j, i] <- mse(ss.fit, u)
  }
}

# set.seed(1)
# x = 1:n/n
# u = smoothwiggly(n,x)
# y = u + rnorm(n,sd=0.005)
# miny = min(y); maxy = max(y)
# y = (y-miny)/(maxy-miny)*10
# u = (u-miny)/(maxy-miny)*10
# 
# tf.fit <- my.trendfilter(y, D.cubic)
# len <- length(tf.fit$lambda)
# tf.mse <- rep(NA, len)
# for (i in 1 : len)
#   tf.mse[i] <- mse(tf.fit$primal.fit[, i], u)
# cat("min mse:", min(tf.mse), "df:", tf.fit$df[which.min(tf.mse)])

# Compute the standard error 
hill.tf.MSE.mean <- apply(hill.tf.MSE, 1, mean)
hill.tf.MSE.sd <- apply(hill.tf.MSE, 1, sd)
hill.tf.MSE.u <- hill.tf.MSE.mean + hill.tf.MSE.sd
hill.tf.MSE.l <- hill.tf.MSE.mean - hill.tf.MSE.sd

hill.ss.MSE.mean <- apply(hill.ss.MSE, 1, mean)
hill.ss.MSE.sd <- apply(hill.ss.MSE, 1, sd)
hill.ss.MSE.u <- hill.ss.MSE.mean + hill.ss.MSE.sd
hill.ss.MSE.l <- hill.ss.MSE.mean - hill.ss.MSE.sd


# Plot the MSE
pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/Figure6.pdf", w = 10, h = 5)
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2) + 2)
#par(mar = c(2, 2, 2, 2) + 2)
plot(1, xlab = "Degrees of Freedom", ylab = "Square Error", main = "Hills Example", xlim = c(5, 40), ylim = c(0.01, 1.5), log = "y", type = "n")
lines(df.seq, hill.tf.MSE.mean, col = "red", type = "l")
points(df.seq, hill.tf.MSE.mean, col = "red", pch = 16, cex = 0.5)
lines(df.seq, hill.tf.MSE.l, col = "red", lty = 2)
lines(df.seq, hill.tf.MSE.u, col = "red", lty = 2)
lines(df.seq, hill.ss.MSE.mean, col = "blue")
points(df.seq, hill.ss.MSE.mean, col = "blue", pch = 16, cex = 0.5)
lines(df.seq, hill.ss.MSE.l, col = "blue", lty = 2)
lines(df.seq, hill.ss.MSE.u, col = "blue", lty = 2)
legend("topright", c("Trend Filtering", "Smoothing splines"), col = c("red", "blue"), lty = c(1, 1), cex = 0.9)
### End of reproducing the MSE result of the "hills" data

### Reproducing the MSE result of the Doppler example
n <- 1000; nsim <- 20; k <- 3; n.dfs <- 80
df.seq <- seq(from = 5, to = 100, length.out = n.dfs)

dop.tf.MSE <- dop.ss.MSE <- matrix(NA, nrow = n.dfs, ncol = nsim)
D.cubic <- difference.matrix(n, k + 1)
val.ind <- which(seq(n)/n >= 0.175)
set.seed(0)
seeds <- sample(1:1000, nsim)

for (i in 1:nsim)
{
  set.seed(seeds[i])
  cat(i, fill = F)
  # generate the data
  doppler.signal <- generate.doppler(n, noise = 0.1)


  # fit trend filtering
  tf.fit <- my.trendfilter(doppler.signal$y, D.cubic, maxsteps = 5000)

  for (j in 1 : n.dfs)
  {
    index <- median(which(abs(tf.fit$df - round(df.seq[j])) == min(abs(tf.fit$df - round(df.seq[j])))))
    dop.tf.MSE[j, i] <- mse(tf.fit$primal.fit[, index][val.ind], doppler.signal$means[val.ind])


    #fit smoothing splines
    ss.fit <- smooth.spline(doppler.signal$x, doppler.signal$y, df = df.seq[j])$y
    dop.ss.MSE[j, i] <- mse(ss.fit[val.ind], doppler.signal$means[val.ind])
  }
}

# Compute the standard error 
dop.tf.MSE.mean <- apply(dop.tf.MSE, 1, mean)
dop.tf.MSE.sd <- apply(dop.tf.MSE, 1, sd)
dop.tf.MSE.u <- dop.tf.MSE.mean + dop.tf.MSE.sd
dop.tf.MSE.l <- dop.tf.MSE.mean - dop.tf.MSE.sd

dop.ss.MSE.mean <- apply(dop.ss.MSE, 1, mean)
dop.ss.MSE.sd <- apply(dop.ss.MSE, 1, sd)
dop.ss.MSE.u <- dop.ss.MSE.mean + dop.ss.MSE.sd
dop.ss.MSE.l <- dop.ss.MSE.mean - dop.ss.MSE.sd


# Plot the MSE
#pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/Figure6.pdf", w = 10, h = 10)
#par(mfrow = c(1, 2), mar = c(2, 2, 2, 2) + 2)
#par(mar = c(2, 2, 2, 2) + 2)
plot(1, xlab = "Degrees of Freedom", ylab = "Square Error", main = "Doppler Example", xlim = c(0, 100), ylim = c(2e-4, 2e-1), log = "y", type = "n")
lines(df.seq, dop.tf.MSE.mean, col = "red", type = "l")
points(df.seq, dop.tf.MSE.mean, col = "red", pch = 16, cex = 0.5)
lines(df.seq, dop.tf.MSE.l, col = "red", lty = 2)
lines(df.seq, dop.tf.MSE.u, col = "red", lty = 2)
lines(df.seq, dop.ss.MSE.mean, col = "blue")
points(df.seq, dop.ss.MSE.mean, col = "blue", pch = 16, cex = 0.5)
lines(df.seq, dop.ss.MSE.l, col = "blue", lty = 2)
lines(df.seq, dop.ss.MSE.u, col = "blue", lty = 2)
legend("topright", c("Trend Filtering", "Smoothing splines"), col = c("red", "blue"), lty = c(1, 1), cex = 0.9)

dev.off()
### End of reproducing the MSE result of the Doppler example
