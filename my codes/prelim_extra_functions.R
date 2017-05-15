##################################################
### Functions used in the prelim
##################################################
source("~/Documents/UW/course work/17Spring/STAT 572/my codes/lasso_solvers.R")


### Generate k-th order difference matrix
difference.operator <- function(a)
{
  a.dim <- dim(a)
  for (j in a.dim[2] : 2)
    a[, j] <- a[, j] - a[, j - 1]
  return (a[-a.dim[1], ])
}

difference.matrix <- function(n, k)
{
  if (k >= n - 1)
    stop("Cannot take such high order of difference")
  a <- diag(rep(1, n))
  for (j in 1 : k) a <- difference.operator(a)
  return (a)
}
### End of generating k-th order difference matrix


### Truncated power basis for splines
#for single value x
trun.basis.splines.single <- function(M, knots, x, coefficients)
{
  n.knots <- length(knots)
  if (length(coefficients) != (M + n.knots))
    stop("Number of basis functions and coefficients don't match")
  mean <- 0
  for (i in 1 : M) mean = mean + coefficients[i] * x ^ (i - 1)
  for (i in 1 : n.knots) mean = mean + coefficients[i + M] * ifelse(x >= knots[i], (x - knots[i])^(M - 1), 0)
  return (mean)
}
#for a sequence of x
trun.basis.splines <- function(M, knots, x, coefficients)
{
  mean <- sapply(1 : length(x), function(i) trun.basis.splines.single(M, knots, x[i], coefficients))
  return (mean)
}

### Generate splines signals using truncated power basis
gen.splines.signals <- function(n, M, n.knots, knots = NULL, noise, coefficients = NULL)
{
  x <- seq(n) / n
  if (!is.null(knots))
  {
    if (n.knots != length(knots))
      stop("Inconsistency in number of knots")
  }
  if (is.null(knots)) knots <- runif(n.knots, min = 0, max = 1)
  if (is.null(coefficients)) coefficients <- rnorm(M + n.knots, sd = 1e9)
  means <- trun.basis.splines(M, knots, x, coefficients)
  means <- scale(means)
  y <- means+rnorm(n, mean = 0, sd = noise)
  return(list(x=x, y=y,means=means,knots=knots))
}
### End of generating splines signals using truncated power basis

### Generate signal using B-spline basis
gen.signal <- function(x, M, knots, noise)
{
  basis <- bs(x, degree = M, knots = knots, intercept = TRUE)
  coef <- rnorm(ncol(basis), sd = 500)
  means <- basis %*% coef
  means <- scale(means)
  y <- means + rnorm(length(x), mean = 0, sd = noise)
  return (list(x=x, y = y, means = means, knots = knots))
}
### End of generating signal using B-spline basis

### Fit piecewise polynomial models of order M
piecewise.polynomial <- function(signal, M)
{
  knots <- signal$knots
  x <- signal$x
  y <- signal$y
  n <- length(x)
  knots <- c(x[1], knots, x[n])
  n.knots <- length(knots)
  index <- fit <- vector("list", length(n.knots - 1))
  for (j in 1 : (n.knots - 1))
  {
    y.piece <- y[(knots[j] <= x) & (x <= knots[j+1])]
    x.piece <- x[(knots[j] <= x) & (x <= knots[j+1])]
    fit[[j]] <- lm(y.piece ~ poly(x.piece, M))$fitted.values
    index[[j]] <- seq(n)[(knots[j] <= x) & (x <= knots[j+1])]
  }
  signal$fit <- fit
  signal$piece.index <- index
  return (signal)
}
### End of fitting piecewise polynomial models of order M

### Fit splines of order M
splines <- function(x, y, M, knots)
{
  basis <- bs(x, degree = M, knots = knots, intercept = TRUE)
  est.coef <- solve(t(basis) %*% basis) %*% t(basis) %*% y
  fit <- basis %*% est.coef
  return (fit)
}
### End of fitting splines of order M


### My trend filtering solver
backsolveSparse <- function (QR, b) 
{
  R = qr.R(QR)
  x = solve(R, qr.qty(QR, b)[seq(1, nrow(R))])
  if (length(QR@q) == 0) 
    return(x)
  else return(x[order(QR@q + 1)])
}

my.trendfilter <- function(y, D, btol = 1e-7, maxsteps = 2000, verbose = FALSE)
{
  #initialization 
  m <- nrow(D)
  n <- ncol(D)
  
  # first variable to hit the Boundary
  buf <- min(maxsteps, 1000)
  dual.fit <- matrix(0, m, buf)
  lambda <- numeric(buf)
  if.hit <- logical(buf)
  df <- numeric(buf)
  
  D <- Matrix(D, sparse = T)
  x = qr(t(D))
  uhat <- backsolveSparse(x, y)
  dual.fit[, 1] <- uhat
  hit.lambda <- max(abs(uhat))
  hit.index <- which.max(abs(uhat))
  lambda[1] <- hit.lambda
  B <- hit.index
  I <- seq(1, m)[-hit.index]
  df[1] <- n - m
  sgn <- sign(uhat[hit.index])
  if.hit[1] <- T
  r <- 1
  k <- 1
  Ds = D[hit.index, ] * sgn
  D1 = D[-hit.index, , drop = FALSE]
  D2 = D[hit.index, , drop = FALSE]
  if (verbose) {
    cat(sprintf("\n%i. lambda=%.3f, adding coordinate %i, |B|=%i...", 
                k, hit.lambda, B[r], r))
  }
  k <- k + 1
  
  # Start down the path
  while (lambda[k-1] > 1e-5 && k <= maxsteps)
  {
    # out of space, enlarge the data structure
    if (k > length(lambda)) {
      buf <- length(lambda)
      lambda <- c(lambda, numeric(buf))
      if.hit <- c(if.hit, logical(buf))
      df <- c(df, numeric(buf))
      dual.fit <- cbind(dual.fit, matrix(0, m, buf))
    }
    x = qr(t(D1))
    # if the bounary set is already full, only need to check the leaving time
    if (r == m)
    {
      a = b = numeric(0)
      hit.lambda = 0
    }
    else {
      a <- backsolveSparse(x, y)
      b <- backsolveSparse(x, Ds)
      shits <- sign(a)
      hits <- a / (b + shits)
      hits[hits > lambda[k-1] + btol] <- 0
      hits[hits > lambda[k-1]] <- lambda[k-1]
      hit.index <- which.max(hits)
      hit.lambda <- max(hits)
      shit <- shits[hit.index]
    }
    
    # if the boundary is empty, only need to check the hitting time
    if (r == 0)
    {
      leave.lambda <- 0
    }
    else {
      cc = as.numeric(sgn * (D2 %*% (y - t(D1) %*% a)))
      dd = as.numeric(sgn * (D2 %*% (Ds - t(D1) %*% b)))
      leaves <- cc/dd
      leaves[cc >= 0] <- 0
      leaves[leaves > lambda[k-1] + btol] <- 0
      leaves[leaves > lambda[k-1]] <- lambda[k-1]
      leave.index <- which.max(leaves)
      leave.lambda <- leaves[leave.index]
    }
    
    # check hit/leave
    if (hit.lambda <= 0 && leave.lambda <= 0) 
      break
    # if hit
    if (hit.lambda > leave.lambda) {
      uhat = numeric(m)
      uhat[B] = hit.lambda * sgn
      uhat[I] = a - hit.lambda * b
      df[k] <- n - m + r
      lambda[k] <- hit.lambda
      r <- r + 1
      B <- c(B, I[hit.index])
      I <- I[-hit.index]
      if.hit[k] <- T
      Ds <- Ds + D1[hit.index, ] * shit
      sgn <- c(sgn, shit)
      D2 <- rbind(D2, D1[hit.index, ])
      D1 <- D1[-hit.index, , drop = FALSE]
      if (verbose) {
        cat(sprintf("\n%i. lambda=%.3f, adding coordinate %i, |B|=%i...", 
                    k, hit.lambda, B[r], r))
      }
    }
    else {
      uhat = numeric(m)
      uhat[B] = leave.lambda * sgn
      uhat[I] = a - leave.lambda * b
      df[k] <- n - m + r
      lambda[k] <- leave.lambda
      r <- r - 1
      I = c(I, B[leave.index])
      B <- B[-leave.index]
      if.hit[k] <- F
      Ds = Ds - D2[leave.index, ] * sgn[leave.index]
      sgn <- sgn[-leave.index]
      D1 <- rbind(D1, D2[leave.index, ])
      D2 <- D2[-leave.index, , drop = FALSE]
      if (verbose) {
        cat(sprintf("\n%i. lambda=%.3f, deleting coordinate %i, |B|=%i...", 
                    k, leave.lambda, I[m - r], r))
      }
    }
    dual.fit[, k] <- uhat
    k <- k + 1
  }
  
  # compute the primal estimate
  primal.fit <- y - t(D) %*% dual.fit[, 1:(k-1)]
  dual.fit <- as.matrix(dual.fit[, 1:(k-1)]); primal.fit <- as.matrix(primal.fit)
  lambda <- lambda[1:(k-1)]
  if.hit <- if.hit[1:(k-1)]
  df <- df[1:(k-1)]
  colnames(dual.fit) <- round(lambda, 4)
  colnames(primal.fit) <- round(lambda, 4)
  
  return(list(primal.fit = primal.fit, dual.fit = dual.fit, lambda = lambda, hit = if.hit, df = df, bls = y))
}
### End of trend filtering solver


### Test the path algorithm for trend filtering
set.seed(2)
n <- 1000; k <- 3
#cubic.signal <- gen.signal(seq(n)/n, M = 4, knots = c(0.1, 0.3, 0.4, 0.7, 0.9), noise = 0.02)
cubic.signal <- generate.hills(seq(n)/n, M = 3, knots = c(0.5, 0.69, 0.7, 0.9), noise = 0.5)
y <- cubic.signal$y; x <- diag(rep(1, n))
D <- difference.matrix(n, k + 1)

ptm <- proc.time()
my.path <- generalized.lasso.v3(y, D, verbose = F, maxsteps = 8000)
proc.time() - ptm
ptm <- proc.time()
ryan.path <- trendfilter(y, ord = 3, verbose = F, maxsteps = 8000)
proc.time() - ptm
ryan.lambda <- ryan.path$lambda
for (i in 1:length(ryan.lambda))
{
  kkt(x=x, y=y, beta=ryan.path$beta[, i], lambda=ryan.lambda[i])
}
### Finished testing the path algorithm


### Generate the design matrix for locally adaptive regression splines (with truncated power basis)
generate.LARS.design <- function(n, k)
{
  design <- matrix(NA, nrow = n, ncol = n)
  for (i in 1 : n)
    for (j in 1 : n)
    {
      if (k == 0)
      {
        if (i < j) design[i, j] <- 0
        if (i >= j) design[i, j] <- 1
      }
      else
        if (k %% 2 == 0)
        {
          if (j <= (k + 1)) design[i, j] <- (i/n)^(j-1)
          if (i <= (j - k/2) && j >= (k+2)) design[i, j] <- 0
          if (i > (j - k/2) && j >= (k+2)) design[i, j] <- ((i - j + k/2) / n)^k
        }
        else
        {
          if (j <= (k + 1)) design[i, j] <- (i/n)^(j-1)
          if (i <= j - (k+1)/2 && j >= (k+2)) design[i, j] <- 0
          if (i > j - (k+1)/2 && j >= (k+2)) design[i, j] <- (i - j + (k+1)/2)^k / n^k
        }
    }
  return (design)
}
### End of generating the design matrix for locally adaptive regression splines (with truncated power basis)


### Generate the design matrix for TF using formula (23) & (24)
cum.sum <- function(n, k)
{
  cs <- list()
  cs[[1]] <- rep(1, n)
  if (k == 0) return (cs)
  for (i in 1 : k)
    cs[[i + 1]] <- cumsum(cs[[i]])
  return (cs)
}

generate.TF.H.design <- function(n, k)
{
  cs <- cum.sum(n, k)
  design <- matrix(NA, nrow = n, ncol = n)
  for (i in 1 : n)
  {
    for (j in 1 : n)
    {
      if (j <= (k+1)) design[i, j] <- (i/n)^(j - 1)
      if (i <= (j - 1) && j >= (k+2)) design[i, j] <- 0
      if (i > (j - 1) && j >= (k+2)) design[i, j] <- cs[[k + 1]][i - j + 1] * factorial(k)/ n^k
    }
  }
  return (design)
}
### End of generating the design matrix for TF using formula (23) & (24)


### Generate the design matrix for TF using formula of M in the supplement
M.k <- function(n, k)
{
  a <- matrix(0, nrow = n, ncol = n)
  a[1 : k, 1 : k] <- diag(rep(1, k))
  b <- matrix(rep(1, (n - k)^2), nrow = (n-k), ncol = (n-k))
  b[upper.tri(b)] <- 0
  a[(k+1):n, (k+1):n] <- b
  return (a)
}

generate.TF.M.design <- function(n, k)
{
  design <- matrix(rep(1, n^2), nrow = n)
  design[upper.tri(design)] <- 0
  if (k == 0) return (design)
  for (j in 1 : k) design <- design %*% M.k(n, j)
  return (design * factorial(k) / n^k)
  #return(design)
}
### End of generating the design matrix for TF using formula of M in the supplement 

# Solver the non-standard lasso problem
non.stan.lasso <- function(X, y, k, tol = 1e-10)
{
  n <- length(y)
  X.1 <- X[, 1 : (k + 1)]
  X.2 <- X[, (k + 2) : n]
  new.y <- y - X.1 %*% solve(t(X.1) %*% X.1, t(X.1) %*% y)
  new.x <- X.2 - X.1 %*% solve(t(X.1) %*% X.1, t(X.1) %*% X.2)
  cat("Begin solving Lasso", fill = F)
  lars.sol <- lars(new.x, new.y, type = "lasso", normalize = F, intercept = F, eps = tol)
  cat("Finished solving Lasso", fill = F)

  lambda <- c(lars.sol$lambda, 0)
  #alpha2.hat should be (n - k - 1) * (#lambda)
  alpha2.hat <- t(lars.sol$beta)
  #alpha1.hat should be (k + 1) * (#lambda)
  alpha1.hat <- solve(t(X.1) %*% X.1, t(X.1) %*% (y - X.2 %*% alpha2.hat))
  #alpha.hat should be n * (#lambda), each column corresponds to each lambda, last one being LS solution
  alpha.hat <- rbind(alpha1.hat, alpha2.hat)
  beta.hat <- X %*% alpha.hat
  return (list(lambda = lambda, beta = beta.hat))
}

### Generate the "hills" data (outdated, use Ryan's hill's data instead)
generate.hills <- function(x, M, knots, noise)
{
  n <- length(x)
  signal <- gen.signal(x, M, knots, noise)
  n.knots <- length(knots)
  index <- which(x > knots[n.knots])
  for (i in index)
  {
    signal$means[i] <- sin(35 * signal$x[i]) - (sin(35 * knots[n.knots]) - signal$means[index[1] - 1])
  }
  signal$means <- signal$means - mean(signal$means)
  signal$y <- signal$means + rnorm(n, mean = 0, sd = noise)
  return (signal)
}
### End of generating the "hills" data

### Ryan's code for generating the hills data
smoothwiggly.fun = function(x) {
  f = function(a,b,c,x) return(a*(x-b)^2+c)
  fp = function(a,b,c,x) return(2*a*(x-b))
  
  a=-1; b=1/4; c=1;  
  if (x<=1/3) return(f(a,b,c,x))  
  aa=a; bb=b; cc=c; xx=1/3;
  a=1; b=xx-fp(aa,bb,cc,xx)/(2*a); c=f(aa,bb,cc,xx)-a*(xx-b)^2;  
  if (x<=2/3) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=2/3;
  b=0.7; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.775) return(f(a,b,c,x))  
  aa=a; bb=b; cc=c; xx=0.775;
  b=0.8; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.825) return(f(a,b,c,x))  
  aa=a; bb=b; cc=c; xx=0.825;
  b=0.85; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.875) return(f(a,b,c,x))  
  aa=a; bb=b; cc=c; xx=0.875;
  b=0.9; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.925) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=0.925;
  b=0.95; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  if (x<=0.975) return(f(a,b,c,x))
  aa=a; bb=b; cc=c; xx=0.975;
  b=1; a=fp(aa,bb,cc,xx)/(2*(xx-b)); c=f(aa,bb,cc,xx)-a*(xx-b)^2;
  return(f(a,b,c,x))
}

smoothwiggly = function(n=100, x=1:n/n) {
  u = rep(0,n)
  for (i in 1:n) u[i]=smoothwiggly.fun(x[i])
  return(u)
}
### End of Ryan's code

### Generate the doppler function simulation data
generate.doppler <- function(n, noise)
{
  x <- seq(n)/n
  means <- sin(4/x) + 1.5
  #means <- means - mean(means)
  y <- means + rnorm(n, mean = 0, sd = noise)
  return (list(x=x, means=means, y=y))
}
### End of generating the doppler function simulation data

### Compute the MSE
#compute the prediction mse
mse = function(y1, y2)
{
  if (length(y1) != length(y2))
    stop("Two vectors need to have the same length")
  n = length(y1)
  return (sum((y1 - y2)^2)/n)
}
### End of computing the MSE


