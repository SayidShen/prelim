####################################################
### Reproduce the Figure 8 in TF
####################################################

source("prelim_extra_functions.R")


n = 100

#generate the truncated power basis 
gen.power.basis <- function(n, k)
{
  N <- 300
  x <- seq(0, 1, length.out = N)
  basis <- matrix(NA, nrow = n, ncol = N)
  input <- seq(0, 1, length.out = n)
  if (k %% 2 == 0)
    knots <- input[(k/2+2) : (n-k/2)]
  else
    knots <- input[((k+1)/2 + 1) : (n - (k+1)/2)]
  for (j in 1 : (k + 1))
    basis[j, ] <- x^(j-1)
  for (j in 1 : (n - k - 1))
    basis[(j + k + 1), ] <- (x - knots[j])^k * as.numeric(x >= knots[j])
  return (basis)
}

#generate the falling factorial basis
gen.fall.basis <- function(n, k)
{
  N <- 300
  x <- seq(0, 1, length.out = N)
  basis <- matrix(NA, nrow = n, ncol = N)
  input <- seq(0, 1, length.out = n)
  for (j in 1 : (k+1))
    basis[j, ] <- x^(j-1)
  for (j in 1 : (n - k -1))
  {
    temp <- apply(as.matrix(x), 1, function(x) x - input[(j+1) : (j+k)])
    basis[(j + k + 1), ] <- apply(temp, 2, prod) * as.numeric(x >= input[j + k])
  }
  return (basis)
}

#plot the two basis
#plot the truncated power basis
n <- 22
k <- 3
x <- seq(0, 1, length.out = 300)

pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/Figure8.pdf", w = 10, h = 4)
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
a <- gen.power.basis(n, k)
plot(xlim = c(0, 1), ylim <- c(0, 0.6), xlab = "x", ylab = "y", type = "n", main = "Truncated Power Basis")
sapply((k + 2) : nrow(a), function(i) lines(x, a[i, ], col = "red", lty = 1))
#plot the falling factorial basis
a <- gen.fall.basis(n, k)
plot(xlim = c(0, 1), ylim <- c(0, 0.6), xlab = "x", ylab = "y", type = "n", main = "Falling Factorial Basis")
sapply((k + 2) : nrow(a), function(i) lines(x, a[i, ], col = "blue", lty = 1))
dev.off()

pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/Figure8b.pdf", w = 10, h = 4)
par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
a <- gen.power.basis(n, k)
plot(xlim = c(0.1, 0.3), ylim <- c(0, 1e-3), xlab = "x", ylab = "y", type = "n", main = "Zoomed in Truncated Power Basis")
sapply((k + 2) : nrow(a), function(i) lines(x, a[i, ], col = "red", lty = 1))
#plot the falling factorial basis
a <- gen.fall.basis(n, k)
plot(xlim = c(0.1, 0.3), ylim <- c(0, 1e-3), xlab = "x", ylab = "y", type = "n", main = "Zoomed in Falling Factorial Basis")
sapply((k + 2) : nrow(a), function(i) lines(x, a[i, ], col = "blue", lty = 1))
dev.off()


pdf("~/Documents/UW/course work/17Spring/STAT 572/my report/Figures/Figure8_beamer.pdf", w = 8, h = 8)
par(mfrow = c(2, 2), mar = c(2, 2, 2, 2))
a <- gen.power.basis(n, k)
plot(xlim = c(0, 1), ylim <- c(0, 0.6), xlab = "x", ylab = "y", type = "n", main = "Truncated Power Basis")
sapply((k + 2) : nrow(a), function(i) lines(x, a[i, ], col = "red", lty = 1))
#plot the falling factorial basis
a <- gen.fall.basis(n, k)
plot(xlim = c(0, 1), ylim <- c(0, 0.6), xlab = "x", ylab = "y", type = "n", main = "Falling Factorial Basis")
sapply((k + 2) : nrow(a), function(i) lines(x, a[i, ], col = "blue", lty = 1))

a <- gen.power.basis(n, k)
plot(xlim = c(0.1, 0.3), ylim <- c(0, 1e-3), xlab = "x", ylab = "y", type = "n", main = "Zoomed in Truncated Power Basis")
sapply((k + 2) : nrow(a), function(i) lines(x, a[i, ], col = "red", lty = 1))
#plot the falling factorial basis
a <- gen.fall.basis(n, k)
plot(xlim = c(0.1, 0.3), ylim <- c(0, 1e-3), xlab = "x", ylab = "y", type = "n", main = "Zoomed in Falling Factorial Basis")
sapply((k + 2) : nrow(a), function(i) lines(x, a[i, ], col = "blue", lty = 1))
dev.off()

