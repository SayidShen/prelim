##################################################
### Lasso Solver Used in Prelim
##################################################

##################################################
### Common Lasso functions
##################################################

###solve the lasso problem 1/2N * ||Y-Xb||_2^2 + lambda * ||b||_1

#compute the lasso loss function value
lasso.loss = function(x, y, beta)
{
  n = nrow(x)
  return (sum((y - x%*%beta)^2)/(2 * n) )
}

#compute the lasso regularization term function value
lasso.regularization = function(lambda, beta)
{
  return (lambda * sum(abs(beta)))
}

#compute the lasso loss function derivative
lasso.loss.gradient = function(x, y, beta)
{
  x = as.matrix(x)
  n = nrow(x)
  return (1/n * t(x) %*% (x %*% beta - y))
}

#compute the prediction mse
lasso.mse = function(y1, y2)
{
  if (length(y1) != length(y2))
    stop("Two vectors need to have the same length")
  n = length(y1)
  return (sum((y1 - y2)^2)/n)
}

#backtracking for the PGA/APGA algorithm for lasso
lasso.back.tracking = function(x, y, beta, lambda, alpha = 0.5, gamma = 0.8, s = 1)
{
  current.gradient = lasso.loss.gradient(x, y, beta) 
  n <- length(y)
  while (1)
  {
    new.beta = beta - s * current.gradient * n
    new.beta = soft.threshold(new.beta, s * lambda * n)
    G = 1/s * (beta - new.beta)
    old.obj.value = lasso.loss(x, y, beta)
    new.obj.value = lasso.loss(x, y, new.beta)
    #stop criterion
    if (new.obj.value <= old.obj.value - s * sum(current.gradient * G) + s/2 * sum(G^2))
      break
    s = s * gamma
  }
  return (s)
}

#the soft-threshold operator
soft.threshold = function(x, lambda)
{
  if (lambda <= 0)
    stop("The threshold value must be positive")
  xx <- rep(NA, length(x))
  for (i in 1:length(x))
    xx[i] <- ifelse(x[i] > lambda, x[i] - lambda, ifelse(x[i] < -lambda, x[i] + lambda, 0))
  return (xx)
}

#compute the lasso objective value
lasso.object.value = function(x, y, beta, lambda)
{
  n = nrow(x)
  return (sum((y - x%*%beta)^2)/(2 * n) + lambda * sum(abs(beta)))
}


### End of common lasso functions

### Coordinate Descent Algorithm (doesn't work at the moment)

lasso.coordinate.descent = function(x, y, lambda, tolerance)
{
  n <- nrow(x)
  p <- ncol(x)
  x.mean <- apply(x, 2, mean)
  for (j in 1:p)
    x[, j] <- x[, j] - x.mean[j]
  x.norm <- apply(x, 2, function(xx){sqrt(mean(xx^2))})
  for (j in 1:p)
      x[, j] <- x[, j] / x.norm[j]
  #y.sd <- sd(y)
  #y <- scale(y, T, T)
  y <- y - mean(y)

  #using random beta initial values
  old.beta <- new.beta <- rep(0, p)

  while (1)
  {
    for (j in 1:p)
    {
      r = y - x %*% new.beta 
      new.beta[j] = soft.threshold(new.beta[j] + sum(x[, j] * r)/n, lambda)
    }
    #cat("difference: ", max(abs(new.beta - old.beta)), "\n")
    if (max(abs(new.beta - old.beta)) < tolerance)
      break
    old.beta <- new.beta
  }
  new.beta <- new.beta / x.norm
  return (new.beta)
}
### End of coordinate descent solver for the Lasso






#test a random design
# set.seed(20)
# n <- 100; p <- 30; lambda <- 50
# x <- matrix(rnorm(n * p, sd = 20), nrow = n, ncol = p)
# beta0 <- runif(p, min = -3, max = 5)
# y <- x %*% beta0 + rnorm(n, sd = 5)
# sean.sol <- as.vector(penalizedLasso(y = y, X = x, penalty = lambda, initialGuess = rep(0, p), maxIter = 1000, tol = 1e-7))
# glmnet.sol <- as.vector(glmnet(x = x, y = y, standardize = F, standardize.response = F, intercept = F, lambda = lambda, thresh = 1e-20)$beta)
# coor.sol <- as.vector(lasso.coordinate.descent(x = x, y = y, lambda = lambda, tolerance = 1e-8))
# bruce.sol <- as.vector(ccd(x = x, y = y, lamb = lambda))
# pga.sol <- as.vector(single.PGA.lasso(x = x, y = y, beta = rep(0, p), lambda = lambda, tolerance = 1e-10))
# apga.sol <- as.vector(single.APGA.lasso(x = x, y = y, beta = rep(0, p), lambda = lambda, tolerance = 1e-10))
# kkt(x,y,glmnet.sol, lambda)
# lasso.object.value(x,y,glmnet.sol,lambda)
# kkt(x,y,coor.sol, lambda)
# lasso.object.value(x,y,coor.sol,lambda)
# kkt(x,y,bruce.sol, lambda)
# lasso.object.value(x,y,bruce.sol,lambda)
# kkt(x,y,pga.sol, lambda)
# lasso.object.value(x,y,pga.sol,lambda)
# kkt(x,y,apga.sol,lambda)
# lasso.object.value(x,y,apga.sol,lambda)
# View(cbind(glmnet.sol, coor.sol, pga.sol, apga.sol))
# End of testing 


### Proximal Gradient Descent

#use PGA to solve lasso, for a single lambda (works fine)
single.PGA.lasso = function(x, y, beta, lambda, tolerance)
{
  #standardize the input
  n <- nrow(x)
  p <- ncol(x)
  
  #x.norm <- apply(x, 2, function(xx){sqrt(sum(xx^2)/n)})
  #x.norm <- apply(x, 2, sd)
  #for (j in 1:p)
  #  x[, j] <- x[, j]/x.norm[j]
  #y.sd <- sd(y)
  #y <- as.vector(scale(y, T, T))
  #browser()
  old.beta = beta
  #obj.val = old.obj.val = sum((y - x%*%beta)^2)/(2 * n)  + lambda * sum(abs(beta))
  step.size = 1/ eigen(2*t(x)%*%x)$values[1]
  while (1)
  {
    #step.size = lasso.back.tracking(x, y, old.beta, lambda)
    #cat(step.size, "\n")
    old.gradient = 1/n * t(x) %*% (x %*% old.beta - y)
    new.beta = old.beta - step.size * old.gradient * n
    new.beta = soft.threshold(new.beta, step.size * lambda * n)
    #new.obj.val = sum((y - x%*%new.beta)^2)/(2 * n)  + lambda * sum(abs(new.beta))
    #obj.val = c(obj.val, new.obj.val)
    #if (abs(new.obj.val - old.obj.val) < tolerance)
    #  break
    if (max(abs(new.beta - old.beta)) < tolerance) break
    old.beta <- new.beta
    #old.obj.val <- new.obj.val
  }
  #browser()
  #new.beta <- new.beta / x.norm

  return(new.beta)
}

PGA.lasso = function(x, y, beta, lambda.seq, tolerance = 1e-5)
{
  x = as.matrix(x)
  if (nrow(x) != length(y))
    stop("Design matrix X and response variable must have the same number of samples\n")
  n = nrow(x)
  p = ncol(x)
  n.lambda = length(lambda.seq)
  coefficient = matrix(NA, nrow = n.lambda, ncol = p)
  obj.val = vector("list", n.lambda)
  beta.hat = beta
  for (i in 1 : n.lambda)
  {
    cat("Begin job ", i, "\n")
    result = single.PGA.lasso(x, y, beta.hat, lambda.seq[i], tolerance)
    beta.hat = result$beta.hat
    coefficient[i, ] = beta.hat
    obj.val[[i]] = result$obj.val
  }
  
  return.list = list()
  return.list$x = x
  return.list$y = y
  return.list$lambda.seq = lambda.seq
  return.list$coefficient = coefficient
  return.list$obj.val = obj.val
  
  return (return.list)
}

### End of proximal gradient algorithm

### Accelerated Proximal gradient algorithm
#use APGA to solve lasso, for a single lambda (works fine)
single.APGA.lasso = function(x, y, beta, lambda, tolerance)
{
  #standardize the input
  #x.norm <- apply(x, 2, function(xx){sqrt(sum(xx^2)/n)})
  #for (j in 1:p)
  #  x[, j] <- x[, j]/x.norm[j]
  
  old.beta = old.y = beta;
  old.lambda = 1
  #obj.val = old.obj.val = sum((y - x%*%beta)^2)/(2 * n) + lambda * sum(abs(beta))
  step.size = 1/ eigen(2*t(x)%*%x)$values[1]
  while (1)
  {
    #step.size = lasso.back.tracking(x, y, old.y, lambda)
    old.gradient = 1/n * t(x) %*% (x %*% old.y - y)
    new.beta = old.y - step.size * old.gradient * n
    new.beta = soft.threshold(new.beta, step.size * lambda * n)
    #new.obj.val = sum((y - x%*%new.beta)^2)/(2 * n) + lambda * sum(abs(new.beta))
    new.lambda = (1 + sqrt(1 + 4 * old.lambda^2))/2
    new.y = new.beta + (old.lambda - 1) / new.lambda * (new.beta - old.beta)
    #obj.val = c(obj.val, new.obj.val)
    #if (abs(new.obj.val - old.obj.val) < tolerance)
    #  break
    if (max(abs(new.beta - old.beta)) < tolerance) break
    old.beta <- new.beta
    old.y <- new.y
    #old.obj.val <- new.obj.val
  }
  #new.beta <- new.beta / x.norm
  
  return(new.beta)
}

APGA.lasso = function(x, y, beta, lambda.seq, tolerance = 1e-3)
{
  x = as.matrix(x)
  if (nrow(x) != length(y))
    stop("Design matrix X and response variable must have the same number of samples\n")
  n = nrow(x)
  p = ncol(x)
  n.lambda = length(lambda.seq)
  coefficient = matrix(NA, nrow = n.lambda, ncol = p)
  obj.val = vector("list", n.lambda)
  beta.hat = beta
  
  for (i in 1 : n.lambda)
  {
    cat("Begin job ", i, "\n")
    result = single.PGA.lasso(x, y, beta.hat, lambda.seq[i], tolerance)
    beta.hat = result$beta.hat
    coefficient[i, ] = beta.hat
    obj.val[[i]] = result$obj.val
  }
  
  return.list = list()
  return.list$coefficient = coefficient
  return.list$obj.val = obj.val
  return.list$x = x
  return.list$y = y
  
  return (return.list)
}
### End of accelerated proximal gradient descent algorithm

### End of lasso solvers


### Test the performance of different algorithms wrt GLMNET
lasso.test <- function(x, y, lambda, tol=1e-8)
{
  glmnet.solution <- as.vector(glmnet(x=x, y=y, family="gaussian", lambda = lambda, standardize = FALSE, standardize.response = FALSE, intercept = FALSE, thresh=tol)$beta)
  cat("Finished GLMNET", fill = FALSE)
  coor.solution <- as.vector(lasso.coordinate.descent(x=x, y=y, lambda=lambda, tolerance = tol))
  cat("Finished Coordinate Descent", fill = FALSE)
  prox.solution <- as.vector(single.PGA.lasso(x=x, y=y, beta = rnorm(ncol(x)), lambda = lambda, tolerance = tol)$beta.hat)
  cat("Finished PGA", fill = FALSE)
  aprox.solution <- as.vector(single.APGA.lasso(x=x, y=y, beta=rnorm(ncol(x)), lambda=lambda, tolerance = tol)$beta.hat)
  cat("Finished APGA", fill = FALSE)
  #kunhui.solution <- as.vector(APGA(x_o=x,y_o=y,lambda=lambda,eps=tol))
  #cat("Finished Kunhui's APGA", "\n")
  #cat("glmnet.solution", glmnet.solution, "\n")
  #cat("coor.solution", coor.solution, "\n")
  #cat("prox.solution", prox.solution, "\n")
  #cat("aprox.solution", aprox.solution, "\n")
  #index <- which(glmnet.solution != 0)
  #cat("max rel diff between coor and glmnet", max(abs(glmnet.solution[index]-coor.solution[index])/abs(glmnet.solution[index])), "\n")
  #cat("max rel diff between prox and glmnet", max(abs(glmnet.solution[index]-prox.solution[index])/abs(glmnet.solution[index])), "\n")
  #cat("max rel diff between aprox and glmnet", max(abs(glmnet.solution[index] - aprox.solution[index])/abs(glmnet.solution[index])), "\n")
  #cat("max rel diff between kunhui's and glmnet", max(abs(glmnet.solution[index] - kunhui.solution[index])/abs(glmnet.solution[index])), "\n")
  #solutions <- data.frame("glmnet" = glmnet.solution, "coor" = coor.solution, "prox"=prox.solution, "aprox" = aprox.solution, "kunhui" = kunhui.solution)
  solutions <- data.frame("glmnet" = glmnet.solution, "coor" = coor.solution, "prox" = prox.solution, "aprox" = aprox.solution)
  #compare the final objective value
  final.obj.val = apply(solutions, 2, lasso.object.value, x = x, y = y, lambda = lambda)
  #cat("GLMNET, coor: ", final.obj.val, "\n")
  return (list(x=x, y=y, lambda=lambda, coef=solutions, final.obj = final.obj.val))
}

#Verify the KKT conditions for the lasso
# 1/n X^T(X * beta - Y) + lambda * S = 0, S is sub-derivative
kkt <- function(x, y, beta, lambda)
{
  n <- length(y)
  s <- t(x) %*% (x %*% beta - y) / (n * lambda)
  if (all(abs(s) <= 1))
    print("KKT satisfied!")
  else
    print("KKT not satisfied!")
  return(s)
}


#test the case for generalized lasso
# set.seed(20)
# n <- 100; lambda <- 0.05
# cubic.signal <- gen.signal(seq(n)/n, M = 3, knots = c(0.2, 0.5), noise = 0.05)
# y <- cubic.signal$y
#   
# H <- generate.TF.M.design(n, k) 
# H.1 <- H[, 1 : (k + 1)]
# H.2 <- H[, (k + 2) : n]
# y <- (diag(rep(1, n)) - H.1 %*% solve(t(H.1) %*% H.1) %*% t(H.1)) %*% y
# x <- (diag(rep(1, n)) - H.1 %*% solve(t(H.1) %*% H.1) %*% t(H.1)) %*% H.2
# 
# solution <- lasso.test(x, y, lambda = lambda, tol = 1e-5)
# s = apply(solution$coef, 2, kkt, x = x, y = y, lambda = lambda)
### End of testing

set.seed(1)
n <- 50; p <- 10; beta <- rnorm(p)
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
y <- X %*% beta + rnorm(n, sd = 0.1)

lars.sol <- lars(X, y, type = "lasso", normalize = F, intercept = F, eps = 1e-7)
glmnet.sol <- glmnet(x=X, y=y, family="gaussian", lambda = lars.sol$lambda / n, standardize = FALSE, standardize.response = FALSE, intercept = FALSE, thresh=tol)
kkt(X, y, beta = lars.sol$beta[7, ], lambda = lars.sol$lambda[7] / n)
