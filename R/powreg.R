#' @export
powreg <- function(y, X, sigma.sq, tau.sq, q, samples =  50000) {
  n <- nrow(X)
  p <- ncol(X)
  
  svd.X <- svd(X, nv = ncol(X))
  U <- svd.X$u
  Vt <- t(svd.X$v)
  d <- c(svd.X$d, rep(0, nrow(Vt) - length(svd.X$d)))
  D <- cbind(diag(svd.X$d), matrix(0, ncol = (nrow(Vt) - length(svd.X$d)), nrow = n))
  
  DUty <- crossprod(crossprod(t(U), D), y)
  W <- crossprod(t(rbind(diag(rep(1, p)), diag(rep(-1, p)))), t(Vt))
  A <- crossprod(X)
  del <- (1 - min(eigen(A)$values))
  start <- as.numeric(crossprod(solve(A + del*diag(ncol(X))), crossprod(X, y)))
  
  # I think this is a replicable way of setting a seed, assuming a seed has been
  # set in R
  return(sampler(DUty = DUty, Vt = Vt, d = d,
                 W = W, sigsq = sig.sq, tausq = tau.sq, q = q, 
                 samples = samples, start = start, seed = rpois(1, 10) + 1))
}