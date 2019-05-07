#' @export
powreg <- function(y, X, sigma.sq, tau.sq, q, samples =  50000, burn = 500, thin = 1, del = NULL,
                   start = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  
  svd.X <- svd(X, nv = ncol(X))
  U <- svd.X$u
  Vt <- t(svd.X$v)
  d <- c(svd.X$d, rep(0, nrow(Vt) - length(svd.X$d)))
  if (nrow(Vt) > length(svd.X$d)) {
    D <- cbind(diag(svd.X$d), matrix(0, ncol = (nrow(Vt) - length(svd.X$d)), nrow = n))
  } else {
    D <- diag(svd.X$d)
  }
  
  DUty <- crossprod(crossprod(t(U), D), y)
  W <- crossprod(t(rbind(diag(rep(1, p)), diag(rep(-1, p)))), t(Vt))
  A <- crossprod(X)
  if (is.null(del)) {
    del <- (1 - min(eigen(A)$values))
  }
  if (is.null(start)) {
    start <- as.numeric(crossprod(solve(A + del*diag(ncol(X))), crossprod(X, y)))
  }
  # If q is bigger than 2, get starting values within unif. dist. bounds
  if (q > 2) {
    start[abs(start) > sqrt(3)*tau.sq] <- sign(start[abs(start) > sqrt(3)*tau.sq])*sqrt(3)*tau.sq
  }
  
  # I think this is a replicable way of setting a seed, assuming a seed has been
  # set in R
  return(sampler(DUty = DUty, Vt = Vt, d = d,
                 W = W, sigsq = sigma.sq, tausq = tau.sq, q = q, 
                 samples = samples, start = start, seed = rpois(1, 10) + 1,
                 burn = burn, thin = thin))
}

## Add functions for estimating tuning parameters
# Function for estimating variance parameters from likelihood
rrmmle<-function(y,X,emu=FALSE,s20=1,t20=1)
{
  sX<-svd(X, nu = nrow(X), nv = ncol(X))
  lX<-sX$d^2
  tUX<-t(sX$u)
  xs<-apply(X,1,sum)
  
  if(nrow(X)>ncol(X))
  {
    lX<-c(lX,rep(0,nrow(X)-ncol(X)))
  }
  
  objective<-function(ls2t2,emu, lX, xs, tUX, y)
  {
    s2<-exp(ls2t2[1]) ; t2<-exp(ls2t2[2])
    mu<-emu*sum((tUX%*%xs)*(tUX%*%y)/(lX*t2+s2))/sum((tUX%*%xs)^2/(lX*t2+s2))
    ev <- lX*t2 + s2
    lev <- log(ifelse(ev > 10^(-300), ev, 1))
    squares <- ifelse(is.infinite((tUX%*%(y-mu*xs))^2/ev), 10^(300),
                      (tUX%*%(y-mu*xs))^2/ev)
    # Useful cat statements for debugging
    # cat("s2=", s2, " t2=", s2, "\n")
    # cat("slev=", sum(lev), " ssquares=", sum(squares), "\n")
    sum(lev) + sum(squares) # + 1/s2^(1/2) # Could keep s2 = 0 from being a mode but is a pretty artificial fix
    
  }
  
  
  # I include an upper bound to keep this from misbehaving, sometimes when
  # s2 and t2 get really small, the gradient goes crazy and a huge
  # value of s2 and t2 can get suggested
  fit<-optim(log(c(s20,t20)),objective,emu=emu, method = "L-BFGS-B", lX = lX,
             xs = xs, tUX = tUX, y = y, upper = rep(log(10^(30)), 2))
  s2<-exp(fit$par[1]) ; t2<-exp(fit$par[2])
  mu<-emu*sum((tUX%*%xs)*(tUX%*%y)/(lX*t2+s2))/sum((tUX%*%xs)^2/(lX*t2+s2))
  if (fit$convergence == 0) {
    return(c(mu,t2,s2))
  } else {
    return(rep(NA, 3))
  }
}

fq <- function(q, kurt) {
  gamma(5/q)*gamma(1/q)/(gamma(3/q)^2) - kurt
}

fpq <- function(q) {
  (gamma(5/q)*gamma(1/q)/(q^2*gamma(3/q)^2))*(6*digamma(3/q) - digamma(1/q) - 5*digamma(5/q))
}

# Use Newton's method: https://en.wikipedia.org/wiki/Newton%27s_method
#' @export
nrq <- function(kurt, sval = 0.032, tol = 10^(-12)) { # This starting value is the lowest possible
  # Kurtosis is bounded below by 1.8, so round if needed
  kurt <- ifelse(kurt <= 1.801, 1.801, kurt)
  # Kurtosis greater than 1.8 gives a q value of 1086.091
  # Value of fpq at q = 1086.091 is about -10^(-8), so the curve *is* pretty flat at this point
  if (kurt < 6 & kurt >= 3) {
    sval <- 1
  } else if (kurt < 3) {
    sval <- 2
  }
  x.old <- Inf
  x.new <- sval
  while (abs(log(x.new/x.old)) > tol) {
    x.old <- x.new
    x.new <- x.old - fq(x.old, kurt)/fpq(x.old)
  }
  return(x.new)
}

# Also need comparison functions to compute Dirichlet Laplace a

fa <- function(a, kurt, p) {
  6*(p*(a*p + 1)*(a + 2)*(a + 3))/((a + 1)*(a*p + 2)*(a*p + 3)) - kurt
}

fpa <- function(a, p) {
  .e1 <- a * p
  .e2 <- 1 + a
  .e3 <- 2 + .e1
  .e4 <- 3 + .e1
  .e6 <- .e2 * .e3 * .e4
  .e7 <- 1 + .e1
  .e8 <- 2 * a
  .e9 <- 2 + a
  .e10 <- 3 + a
  p * (6 * ((.e7 * .e9 + (1 + p * (2 + .e8)) * .e10)/.e6) - 
         6 * (((2 + p * (1 + .e8)) * .e4 + p * .e2 * .e3) * .e7 * 
                .e9 * .e10/.e6^2))
}

#' Some functions for doing variance component estimation
s.sq.r.sq <- function(r.sq, y.tilde = NULL, d = NULL, y = NULL, X = NULL) {
  if (is.null(y.tilde) & is.null(d)) {
    svd <- svd(X, nu = nrow(X))
    U <- svd$u
    y.tilde <- crossprod(U, y)
    d <- c(svd$d, rep(0, nrow(X) - min(dim(X))))
  }
  
  s.sq <- numeric(length(r.sq))
  
  for (i in 1:length(r.sq)) {
    
    s.sq[i] <- mean(y.tilde^2/(1 + r.sq[i]*d^2))
  }
  return(s.sq)
}

obj.r.sq <- function(r.sq, y.tilde, d) {
  n <- length(y.tilde)
  # s.sq <- s.sq.r.sq(r.sq = r.sq, y.tilde = y.tilde, d = d)
  # -n*log(s.sq) - sum(log(1 + r.sq*d^2)) - sum(y.tilde^2/(1 + r.sq*d^2))/s.sq
  (-n*log(sum(y.tilde^2/(1 + r.sq*d^2))) - sum(log(1 + r.sq*d^2)) + n*log(n)  - n)/n
}

obj.r.sq.dlr <- function(r.sq, y.tilde, d) {
  n <- length(y.tilde)
  # .e1 <- d^2; .e2 <- 1 + .e1 * r.sq; .e3 <- y.tilde^2; 
  # (-(n * sum(-(.e1 * .e3/.e2^2))/sum(.e3/.e2) + sum(.e1/.e2)))/n
  # (-n*log(sum(y.tilde^2/(1 + exp(lr.sq)*d^2))) - sum(log(1 + exp(lr.sq)*d^2)) + n*log(n)  - n)/n
  lr.sq <- log(r.sq)
  .e1 <- d^2; .e2 <- exp(lr.sq); .e3 <- .e1 * .e2; .e4 <- 1 + .e3; .e5 <- y.tilde^2; -((n * sum(-(.e1 * .e5 * .e2/.e4^2))/sum(.e5/.e4) + sum(.e3/.e4))/n)
}

obj.r.sq.ddlr <- function(r.sq, y.tilde, d) {
  n <- length(y.tilde)
  lr.sq <- log(r.sq)
  # Deriv(Deriv("(-n*log(sum(y.tilde^2/(1 + exp(lr.sq)*d^2))) - sum(log(1 + exp(lr.sq)*d^2)) + n*log(n)  - n)/n", "lr.sq"), "lr.sq")
  .e1 <- d^2; .e2 <- exp(lr.sq); .e3 <- .e1 * .e2; .e4 <- 1 + .e3; .e5 <- y.tilde^2; .e6 <- .e4^2; .e7 <- .e3/.e4; .e8 <- .e1 * .e5; .e9 <- sum(.e5/.e4); -((n * (sum(-(.e8 * (1 - 2 * .e7) * .e2/.e6)) - sum(-(.e8 * .e2/.e6))^2/.e9)/.e9 +      sum(.e1 * (1 - .e7) * .e2/.e4))/n)
}

obj.r.sq.dddlr <- function(r.sq, y.tilde, d) {
  n <- length(y.tilde)
  lr.sq <- log(r.sq)
  # Deriv(Deriv(Deriv("(-n*log(sum(y.tilde^2/(1 + exp(lr.sq)*d^2))) - sum(log(1 + exp(lr.sq)*d^2)) + n*log(n)  - n)/n", "lr.sq"), "lr.sq"), "lr.sq")
  .e1 <- d^2; .e2 <- exp(lr.sq); .e3 <- .e1 * .e2; .e4 <- 1 + .e3; .e5 <- y.tilde^2; 
  .e6 <- .e3/.e4; .e7 <- .e4^2; .e8 <- .e1 * .e5; .e9 <- 1 - 2 * .e6; .e10 <- sum(.e5/.e4); 
  .e12 <- 1 - .e6; .e13 <- sum(-(.e8 * .e2/.e7)); 
  -((n * (sum(-(.e8 * (1 - .e1 * (2 + 2 * .e9 + 2 * .e12) * .e2/.e4) *      .e2/.e7)) - (3 * sum(-(.e8 * .e9 * .e2/.e7)) - 2 * (.e13^2/.e10)) *      .e13/.e10)/.e10 + sum(.e1 * .e9 * .e12 * .e2/.e4))/n)
}

obj.r.sq.ddddlr <- function(r.sq, y.tilde, d) {
  n <- length(y.tilde)
  lr.sq <- log(r.sq)
  .e1 <- d^2; .e2 <- exp(lr.sq); .e3 <- .e1 * .e2; 
  .e4 <- 1 + .e3; .e5 <- .e3/.e4; .e6 <- y.tilde^2; 
  .e7 <- .e4^2; .e8 <- .e1 * .e6; .e9 <- 2 * .e5; .e10 <- 1 - .e9; 
  .e11 <- 1 - .e5; .e12 <- sum(.e6/.e4); .e14 <- sum(-(.e8 * .e2/.e7)); 
  .e16 <- sum(-(.e8 * .e10 * .e2/.e7)); .e18 <- 2 * .e11; .e19 <- 2 + 2 * .e10; 
  .e21 <- .e14^2/.e12; .e22 <- 1 - .e1 * (.e19 + .e18) * .e2/.e4; 
  -((n * (sum(-(.e8 * (1 - .e1 * (2 * .e22 + 4 + 4 * .e10 + 4 *      .e11 - .e1 * (.e19 + 8 * .e11) * .e2/.e4) * .e2/.e4) * .e2/.e7)) -      ((3 * .e16 - 2 * .e21) * .e16 + (4 * sum(-(.e8 * .e22 * .e2/.e7)) -          (2 * (2 * .e16 - .e21) + 6 * .e16 - 4 * .e21) * .e14/.e12) *          .e14)/.e12)/.e12 + sum(.e1 * (.e10 * .e11 - .e1 * (1 +      .e18 - .e9) * .e2/.e4) * .e11 * .e2/.e4))/n)
}

#' Function for computing the profile likelihood of the variance ratio under normal-normal model
#'
#' \code{varcomp}
#' @param \code{r.sq} variance ratio
#' @param \code{y} regression response
#' @param \code{X} regression design matrix
#' @return Objective function value
#' @export
obj.varcomp <- function(r.sq, y = y, X = X, y.tilde = NULL, d = NULL) {
  
  if (is.null(y.tilde) & is.null(d)) {
    n <- nrow(X); p <- ncol(X)
    svd <- svd(X, nu = n)
    U <- svd$u
    y.tilde <- crossprod(U, y)
    d <- c(svd$d, rep(0, n - min(n, p)))
  }
  
  n <- length(y.tilde)
  obj <- numeric(length(r.sq)) 
  
  for (i in 1:length(r.sq)) {
    obj[i] <- obj.r.sq(r.sq = r.sq[i], y.tilde = y.tilde, d = d)
  }
  return(obj)
}

obj.dlr.varcomp <- function(r.sq, y = y, X = X, y.tilde = NULL, d = NULL) {
  
  if (is.null(y.tilde) & is.null(d)) {
    n <- nrow(X); p <- ncol(X)
    svd <- svd(X, nu = n)
    U <- svd$u
    y.tilde <- crossprod(U, y)
    d <- c(svd$d, rep(0, n - min(n, p)))
  }
  
  n <- length(y.tilde)
  obj <- numeric(length(r.sq)) 
  
  for (i in 1:length(r.sq)) {
    obj[i] <- obj.r.sq.dlr(r.sq = r.sq[i], y.tilde = y.tilde, d = d)
  }
  return(obj)
}

obj.ddlr.varcomp <- function(r.sq, y = y, X = X, y.tilde = NULL, d = NULL) {
  
  if (is.null(y.tilde) & is.null(d)) {
    n <- nrow(X); p <- ncol(X)
    svd <- svd(X, nu = n)
    U <- svd$u
    y.tilde <- crossprod(U, y)
    d <- c(svd$d, rep(0, n - min(n, p)))
  }
  
  n <- length(y.tilde)
  obj <- numeric(length(r.sq)) 
  
  for (i in 1:length(r.sq)) {
    obj[i] <- obj.r.sq.ddlr(r.sq = r.sq[i], y.tilde = y.tilde, d = d)
  }
  return(obj)
}

obj.dddlr.varcomp <- function(r.sq, y = y, X = X, y.tilde = NULL, d = NULL) {
  
  if (is.null(y.tilde) & is.null(d)) {
    n <- nrow(X); p <- ncol(X)
    svd <- svd(X, nu = n)
    U <- svd$u
    y.tilde <- crossprod(U, y)
    d <- c(svd$d, rep(0, n - min(n, p)))
  }
  
  n <- length(y.tilde)
  obj <- numeric(length(r.sq)) 
  
  for (i in 1:length(r.sq)) {
    obj[i] <- obj.r.sq.dddlr(r.sq = r.sq[i], y.tilde = y.tilde, d = d)
  }
  return(obj)
}

obj.ddddlr.varcomp <- function(r.sq, y = y, X = X, y.tilde = NULL, d = NULL) {
  
  if (is.null(y.tilde) & is.null(d)) {
    n <- nrow(X); p <- ncol(X)
    svd <- svd(X, nu = n)
    U <- svd$u
    y.tilde <- crossprod(U, y)
    d <- c(svd$d, rep(0, n - min(n, p)))
  }
  
  n <- length(y.tilde)
  obj <- numeric(length(r.sq)) 
  
  for (i in 1:length(r.sq)) {
    obj[i] <- obj.r.sq.ddddlr(r.sq = r.sq[i], y.tilde = y.tilde, d = d)
  }
  return(obj)
}


#' Function for computing maximum likelihood estimates of variance parameters under normal-normal model
#'
#' \code{varcomp}
#'
#' @param \code{y} regression response
#' @param \code{X} regression design matrix
#' @return Estimates \code{sigma.beta.sq.hat}, \code{sigma.epsi.sq.hat}
#' @export
varcomp <- function(y, X, diff.tol, min.sig.sq, min.r.sq, grid.num) {
  
  # Some code for simulations
  # n <- 5
  # p <- 10
  # 
  # X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  # beta <- rnorm(p)
  # y <- X%*%beta + rnorm(n)
  
  # no.noise <- FALSE
  upper.lim <- 100
  
  svd <- svd(X, nu = nrow(X))
  U <- svd$u
  y.tilde <- crossprod(U, y)
  d <- c(svd$d, rep(0, nrow(X) - min(dim(X))))
  
  s.sq.r.sq.val <- s.sq.r.sq(upper.lim, y = y, X = X, y.tilde = y.tilde, d = d)
  
  while (s.sq.r.sq.val > min.sig.sq) {
    s.sq.r.sq.val <- s.sq.r.sq(upper.lim*1.1, y = y, X = X, y.tilde = y.tilde, d = d)
    if (!is.nan(s.sq.r.sq.val)) {
      upper.lim <- upper.lim*1.1
    } else {
      break
    }
  }
  r.sq <- exp(seq(log(min.r.sq), log(upper.lim), length.out = grid.num))
  obj <- obj.varcomp(r.sq, y = y, X = X, y.tilde = y.tilde, d = d)
  
  
  if (max(which(obj == max(obj))) == length(obj)) {
    dderiv <- obj.ddlr.varcomp(r.sq, y = y, X = X, y.tilde = y.tilde, d = d)
    sc <- sign(dderiv[-1]) == sign(dderiv[-length(dderiv)])
    mag <- abs(dderiv[-1] - dderiv[-length(dderiv)])
    if (sum(sc == FALSE) > 1) {
      if (sum(sc == FALSE) == 2) {
        cut <- max(which(sc == FALSE))
      } else {
        which.sc <- which(sc == FALSE)[-1]
        deriv.sc <- abs(obj.dlr.varcomp(r.sq[which.sc], y = y, X = X, y.tilde = y.tilde, d = d))
        cut <- which.sc[deriv.sc == min(deriv.sc)]
      }
      
      r.sq <- exp(seq(log(min.r.sq), log(max(r.sq[cut])), length.out = grid.num))
    } else {
      dddderiv <- obj.ddddlr.varcomp(r.sq, y = y, X = X, y.tilde = y.tilde, d = d)
      sc <- sign(dddderiv[-1]) == sign(dddderiv[-length(dderiv)])
      if (sum(sc == FALSE) == 2) {
        cut <- max(which(sc == FALSE))
      } else {
        which.sc <- which(sc == FALSE)[-1]
        deriv.sc <- abs(obj.dlr.varcomp(r.sq[which.sc], y = y, X = X, y.tilde = y.tilde, d = d))
        cut <- which.sc[deriv.sc == min(deriv.sc)]
      }
      r.sq <- exp(seq(log(min.r.sq), log(max(r.sq[cut])), length.out = grid.num))
    }
    obj <- obj.varcomp(r.sq, y = y, X = X, y.tilde = y.tilde, d = d)
  }
  # r.sq <- exp(seq(log(min.r.sq), log(r.sq[(max(which(deriv <= 0.009)))]), length.out = 1000))
  # obj <- obj.varcomp(r.sq, y = y, X = X, y.tilde = y.tilde, d = d)
  # der <- obj[-1] - obj[-length(obj)]
  # print(range(der))
  # if (min(der) < 0 & min(der) > -10^(-14)) {
  #   no.noise <- TRUE
  # }
  # par(mfrow = c(1, 2))
  # plot(r.sq[r.sq < 2], obj[r.sq < 2], type = "l")
  # plot(r.sq, obj, type = "l")
  # der[abs(der) < 10^(-12)] <- sign(der[abs(der) < 10^(-12)])*10^(-12)
  # der[der == 0] <- 10^(-12)
  # sign.changes <- which(sign(der[-1]) != sign(der[-length(der)]))
  max.obj <- which(obj == max(obj))
  if (length(max.obj) > 1) {
    obj.min <- min(max.obj)
    obj.max <- max(max.obj)
  } else if (max.obj == length(r.sq)) {
    obj.min <- max.obj - 1
    obj.max <- max.obj
  } else if (max.obj == 1) {
    obj.min <- max.obj
    obj.max <- max.obj + 1
  } else {
    obj.min <- max.obj - 1
    obj.max <- max.obj + 1
  }
  max.diff <- max(c(abs(obj[-1] - obj[-length(obj)])), na.rm = TRUE)
  
  while (max.diff > diff.tol) {
    r.sq <- exp(seq(log(r.sq[obj.min]), log(r.sq[obj.max]), length.out = 100))
    obj <- obj.varcomp(r.sq, y = y, X = X)
    max.obj <- which(obj == max(obj))
    if (length(max.obj) > 1) {
      obj.min <- min(max.obj)
      obj.max <- max(max.obj)
    } else if (max.obj == length(r.sq)) {
      obj.min <- max.obj - 1
      obj.max <- max.obj
    } else {
      obj.min <- max.obj - 1
      obj.max <- max.obj + 1
    }
    max.diff <- max(c(abs(obj[-1] - obj[-length(obj)])), na.rm = TRUE)
  }
  
  r.sq.max <- r.sq[min(max.obj)]
  s.sq.max <- s.sq.r.sq(r.sq.max, y = y, X = X)
  tau.sq.max <- r.sq.max*s.sq.max
  
  # print(c(s.sq.max, tau.sq.max))
  
  return(c(tau.sq.max,s.sq.max))
}



#' Function for estimating tuning parameters
#'
#' \code{estRegPars}
#'
#' @param \code{y} regression response
#' @param \code{X} regression design matrix
#' @param \code{delta} ridge regression parameter for when X is not full rank
#' @param \code{diff.tol} tolerance for variance parameter estimation, defaults to diff.tol = 10^(-7) 
#' @return Estimates \code{sigma.beta.sq.hat}, \code{sigma.epsi.sq.hat} and \code{kappa.hat}
#' @export
estRegPars <-function(y, X, delta.sq = NULL, precomp = NULL, comp.q = FALSE, mom = TRUE,
                      diff.tol = 10^(-12), min.sig.sq = 10^(-14), 
                      min.r.sq = 10^(-14), grid.num = 10000) {
  
  
  p <- ncol(X)
  n <- nrow(X)
  
  XtX <- crossprod(X)
  C <- cov2cor(XtX)
  V <- diag(sqrt(diag(XtX/C)))
  ceval <- eigen(C)$values
  if (!is.null(delta.sq)) {
    delta.sq <- delta.sq
  } else {
    delta.sq <- max(1 - min(ceval), 0)
  }
  D.inv <- tcrossprod(crossprod(V, (C + delta.sq*diag(p))), V)
  D <- solve(D.inv)
  DXtX <- crossprod(D, XtX)
  XD <- crossprod(t(X), D)
  DXtXD <- crossprod(XD)
  
  b <- crossprod(D, crossprod(X, y))
  
  if (mom) {
    
    XXt <- tcrossprod(X)
    XXt.eig <- eigen(XXt)
    XXt.val <- XXt.eig$values
    XXt.vec <- XXt.eig$vectors
    A <- diag(n)
    # B <- tcrossprod(tcrossprod(XXt.vec, diag(ifelse(XXt.val > 1, 1/XXt.val, 0))), XXt.vec)
    B <- tcrossprod(tcrossprod(XXt.vec, diag(1/(XXt.val + 1))), XXt.vec)
    
    E <- rbind(c(sum(diag(crossprod(t(XXt), A))), sum(diag(A))), 
               c(sum(diag(crossprod(t(XXt), B))), sum(diag(B))))
    
    sig.2.ests <- solve(E)%*%matrix(c(crossprod(y, crossprod(t(A), y)), 
                                      crossprod(y, crossprod(t(B), y))), nrow = 2, ncol = 1)
    
    sigma.beta.sq.hat <- sig.2.ests[1]
    sigma.epsi.sq.hat <- sig.2.ests[2]
    
  } else {
    vpars <- varcomp(y = y, X = X, diff.tol = diff.tol,
                     min.sig.sq = min.sig.sq, min.r.sq = min.r.sq, grid.num = grid.num) # rrmmle(y = y, X = X)
    sigma.beta.sq.hat <- vpars[1]
    sigma.epsi.sq.hat <- vpars[2]
    # sigma.beta.sq.hat <- vpars[2]
    # sigma.epsi.sq.hat <- vpars[3]
  }
  
  
  alpha.beta <- sum(diag(crossprod(DXtX)))/p
  gamma.beta <- sum(diag(DXtX^4))/p
  omega.beta <- 3*(sum(diag(crossprod(DXtX)^2)) - sum(diag(DXtX^4)))/p
  
  test.stat <- (mean(b^4))/(mean(b^2)^2)
  
  kappa.hat <- (alpha.beta^2/gamma.beta)*(test.stat - omega.beta/alpha.beta^2)
  q.hat <- ifelse(comp.q, nrq(kappa.hat), NA)
  
  return(list("sigma.beta.sq.hat" = sigma.beta.sq.hat,
              "sigma.epsi.sq.hat" = sigma.epsi.sq.hat,
              "kappa.hat" = kappa.hat,
              "q.hat" = q.hat,
              "test.stat" = test.stat,
              "DXtX" = DXtX,
              "DXtXD" = DXtXD,
              "delta.sq" = delta.sq))
  
}