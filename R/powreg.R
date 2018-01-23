#' @export
powreg <- function(y, X, sigma.sq, tau.sq, q, samples =  50000, burn = 500, thin = 1) {
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
  del <- (1 - min(eigen(A)$values))
  start <- as.numeric(crossprod(solve(A + del*diag(ncol(X))), crossprod(X, y)))
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
  kurt <- ifelse(kurt <= 1.8, 1.81, kurt)
  # Kurtosis greater than 1.8 gives a q value of 1086.091
  # Value of fpq at q = 1086.091 is about -10^(-8), so the curve *is* pretty flat at this point
  if (kurt < 6) {
    sval <- 1
  } else if (kurt < 3) {
    sval <- 2
  }
  x.old <- Inf
  x.new <- sval
  while (abs(x.new - x.old) > tol) {
    x.old <- x.new
    x.new <- x.old - fq(x.old, kurt)/fpq(x.old)
  }
  return(x.new)
}
#' Function for estimating tuning parameters
#'
#' \code{estRegPars}
#'
#' @param \code{y} regression response
#' @param \code{X} regression design matrix
#' @param \code{delta} ridge regression parameter for when X is not full rank
#' @return Estimates \code{sigma.beta.sq.hat}, \code{sigma.epsi.sq.hat} and \code{kappa.hat}
#' @export
estRegPars <-function(y, X, delta.sq = NULL, precomp = NULL, comp.q = FALSE, mom = TRUE) {
  
  
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
    vpars <- rrmmle(y = y, X = X)
    sigma.beta.sq.hat <- vpars[2]
    sigma.epsi.sq.hat <- vpars[3]
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