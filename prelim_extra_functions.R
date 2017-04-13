##################################################
### Extra functions used in the prelim
##################################################

#generate k-th order difference matrix
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


#truncated power basis for splines
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

#function to generate splines signals
gen.splines.signals <- function(n, M, n.knots, knots = NULL, noise, coefficients = NULL)
{
  x <- seq(n) / n
  if (!is.null(knots))
  {
    if (n.knots != length(knots))
      stop("Inconsistency in number of knots")
  }
  if (is.null(knots)) knots <- runif(n.knots, min = 0, max = 1)
  if (is.null(coefficients)) coefficients <- rnorm(M + n.knots)
  means <- trun.basis.splines(M, knots, x, coefficients)
  y <- means+rnorm(n, mean = 0, sd = noise)
  return(list(x=x, y=y,means=means,knots=knots))
}

#fit piecewise polynomial models of order M
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

#fit splines of order M
splines <- function(signal, M)
{
  x <- signal$x
  knots <- signal$knots
  y <- signal$y
  basis <- bs(x, degree = M, knots = knots, intercept = TRUE)
  est.coef <- solve(t(basis) %*% basis) %*% t(basis) %*% y
  fit <- basis %*% est.coef
  signal$fit <- fit
  return (signal)
}

#fit smoothing splines
smoothing.splines <- function(signal, df)
{
  x <- signal$x
  y <- signal$y
  fit <- smooth.spline(x, y, df = df)
  signal$fit <- fit$y
  signal$df <- df
  return (signal)
}

#solver of the generalized lasso problem with X = I (signal approximation case)
generalized.lasso <- function(y, D)
{
  #initialization 
  dual.dim <- nrow(D)
  primal.dim <- ncol(D)
  k <- 1
  bound <- sgn <- lambda <- c()
  lambda <- c(lambda, Inf)
  dual.fit <- matrix(nrow = 1, ncol = dual.dim)
  while (lambda[k] > 1e-6)
  {
    cat("k = ", k, "lambda[k] = ", lambda[k], "bound = ", bound, "\n")
    # if the boundary set is still empty, only need to find the hitting time
    if (is.null(bound))
    {
      # solution to the dual problem
      u.hat <- solve(D %*% t(D)) %*% D %*% y
      dual.fit[1, ] <- u.hat
      hit.lambda <- max(abs(u.hat))
      hit.index <- which.max(abs(u.hat))
      lambda <- c(lambda, hit.lambda)
      bound <- c(bound, hit.index)
      sgn <- c(sgn, sign(u.hat[hit.index]))
    }
    else
    {
      u.hat.bound <- lambda[k] * sgn
      u.hat.nonbound <- solve(D[-bound, ] %*% t(D[-bound, ])) %*% D[-bound, ] %*% (y - lambda[k] * t(matrix(D[bound, ], nrow = length(bound))) %*% sgn)
      u.hat <- rep(NA, dual.dim)
      u.hat[bound] <- u.hat.bound; u.hat[-bound] <- u.hat.nonbound
      dual.fit <- rbind(dual.fit, u.hat)
      #find the hitting time
      hit.lambda.plus <- (solve(D[-bound, ] %*% t(D[-bound, ])) %*% D[-bound, ] %*% y) / (solve(D[-bound, ] %*% t(D[-bound, ])) %*% D[-bound, ] %*% t(matrix(D[bound, ], nrow = length(bound))) %*% sgn + 1)
      hit.lambda.minus <- (solve(D[-bound, ] %*% t(D[-bound, ])) %*% D[-bound, ] %*% y) / (solve(D[-bound, ] %*% t(D[-bound, ])) %*% D[-bound, ] %*% t(matrix(D[bound, ], nrow = length(bound))) %*% sgn - 1)
      hit.lambda.real <- rep(NA, dual.dim - length(bound))
      for (i in 1 : (dual.dim - length(bound)))
      {
        if (hit.lambda.plus[i] >= 0 && hit.lambda.plus[i] < lambda[k])
          hit.lambda.real[i] <- hit.lambda.plus[i]
        else
          if (hit.lambda.minus[i] >= 0 && hit.lambda.minus[i] < lambda[k])
            hit.lambda.real[i] <- hit.lambda.minus[i]
          else
            stop("Both stopping plus and stopping minus are out of range!\n")
      }
      #hit.lambda.real <- sapply(1 : (dual.dim - length(bound)), 
      #                          function(i) ifelse(hit.lambda.plus[i] >= 0 && hit.lambda.plus[i] <= lambda[k], hit.lambda.plus[i], hit.lambda.minus[i]))
      hit.lambda <- max(hit.lambda.real)
      if (hit.lambda >= lambda[k] - 1e-6) hit.lambda <- 0
      cat("hitting time: ", hit.lambda, "\n")
      #if (hit.lambda > lambda[k] + 1e-7) stop("Hitting time out of bound!")
      #if (hit.lambda > lambda[k] + 1e-7) hit.lambda <- 0
      hit.index <- which.max(hit.lambda.real)
      
      #find the leave time
      cc <- sgn * (D[bound, ] %*% (diag(rep(1, primal.dim)) - t(D[-bound, ]) %*% solve(D[-bound, ] %*% t(D[-bound, ])) %*% D[-bound, ]) %*% y)
      dd <- sgn * (D[bound, ] %*% (diag(rep(1, primal.dim)) - t(D[-bound, ]) %*% solve(D[-bound, ] %*% t(D[-bound, ])) %*% D[-bound, ]) %*% t(matrix(D[bound, ], nrow = length(bound))) %*% sgn)
      leave.lambda.real <- sapply(1 : length(bound), function(i) ifelse(cc[i] < 0 && dd[i] < 0, cc[i] / dd[i], 0))
      leave.lambda <- max(leave.lambda.real)
      cat("leaving time: ", leave.lambda, "\n")
      #if (leave.lambda > lambda[k] + 1e-7) stop(paste("Leaving time", leave.lambda, "out of bound!"))
      #if (leave.lambda - lambda[k] > 1e-7) leave.lambda <- 0
      leave.index <- which.max(leave.lambda.real)
      
      # check hit/leave
      # if hit
      if (hit.lambda > leave.lambda)
      {
        lambda <- c(lambda, hit.lambda)
        bound <- c(bound, seq(dual.dim)[-bound][hit.index])
        sgn <- c(sgn, sign(u.hat[hit.index]))
      }
      # if leave
      else
      {
        lambda <- c(lambda, leave.lambda)
        bound <- bound[-leave.index]
        sgn <- sgn[-leave.index]
      }
    }
    k <- k + 1
  }
  # compute the primal estimate
  primal.fit <- matrix(NA, nrow = nrow(dual.fit), ncol = primal.dim)
  for (i in 1 : nrow(primal.fit))
    primal.fit[i, ] <- y - t(D) %*% dual.fit[i, ]
  #combine the results
  result <- list()
  result$lambda <- lambda
  result$dual.fit <- dual.fit
  result$primal.fit <- primal.fit
  
  return(result)
}


# Compute the piecewise-quadratic case in Figure 1 of TF paper
set.seed(40)
linear.signal <- gen.splines.signals(n, 2, n.knots = 3, noise = 0.01)
linear.pen.mat <- difference.matrix(n, 2)
linear.tf <- genlasso(linear.signal$y, diag(rep(1, n)), linear.pen.mat)

my.solution <- generalized.lasso(linear.signal$y, D = linear.pen.mat)
