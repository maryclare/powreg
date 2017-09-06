#' @export
powreg <- function(y, X, sigma.sq, tau.sq, q, samples =  50000) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  A <- crossprod(X)
  LDLt <- ldl(A)
  L <- LDLt[["L"]]
  d <- diag(LDLt[["D"]])
  
  Ut <- forwardsolve(L, t(X))
  # X = t(Ut)%*%t(L)
  L.inv <- forwardsolve(L, diag(p))
  
  Uty <- crossprod(t(Ut), y)
  W <- crossprod(t(rbind(diag(rep(1, p)), diag(rep(-1, p)))), t(L.inv))
  del <- (1 - min(eigen(A)$values))
  start <- as.numeric(crossprod(solve(A + del*diag(ncol(X))), crossprod(X, y)))
  # If q is bigger than 2, get starting values within unif. dist. bounds
  if (q > 2) {
    start[abs(start) > sqrt(3)*tau.sq] <- sign(start[abs(start) > sqrt(3)*tau.sq])*sqrt(3)*tau.sq
  }
  
  # I think this is a replicable way of setting a seed, assuming a seed has been
  # set in R
  return(sampler(Uty = Uty, L=L, Linv = L.inv, d = d,
                 W = W, sigsq = sig.sq, tausq = tau.sq, q = q, 
                 samples = samples, start = start, seed = rpois(1, 10) + 1))
}