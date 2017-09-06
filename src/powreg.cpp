// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Scraped from here to get a way to set seed:
// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
void set_seed(unsigned int seed) {
  
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
  
}
// Uses recursive algorithm given here: 
// https://en.wikipedia.org/wiki/Cholesky_decomposition
// [[Rcpp::export]]
List ldl(NumericMatrix A) {
  int p = A.nrow();
  
  NumericMatrix D(p, p);
  NumericMatrix L(p, p);
  
  for(int j = 0; j < p; j++) {
    L(j, j) = 1.0;
    D(j, j) = A(j, j);
    if (j > 0) {
      for (int k = 0; k < j; k++) {
        D(j, j) = D(j, j) - L(j, k)*L(j, k)*D(k, k);
      }
    }
    if (j < p) {
      for (int i = (j + 1); i < p; i++) {
        L(i, j) = A(i, j);
        if (j > 0) {
          for (int k = 0; k < j; k++) {
            L(i, j) = L(i, j) - L(i, k)*L(j, k)*D(k, k);
          }
        }
        if (D(j, j) == 0) {
          L(i, j) = 0;
        } else {
          L(i, j) = L(i, j)/D(j, j);
        }
      }
    }
  }
  return(Rcpp::List::create(Rcpp::Named("L")=L,
                            Rcpp::Named("D")=D));
}

// Implements this code to sample from a truncated univariate normal distribution
// very reliably!
// https://arxiv.org/pdf/0907.4010.pdf
// Could at some point improve this using the following:
// http://www.stat.ncsu.edu/information/library/papers/mimeo2649_Li.pdf
//
// [[Rcpp::export]]
NumericVector rtnormrej(NumericVector mu, NumericVector sd, NumericVector l, NumericVector r) {
  
  int p = mu.size();
  NumericVector lstd = (l - mu)/sd;
  NumericVector rstd = (r - mu)/sd;
  NumericVector u(p, 1.0);
  NumericVector rho(p, 0.0);
  NumericVector z(p, 0.0);
  NumericVector alpha(p, 0.0);
  
  for (int i = 0; i < p; i++) {
    
    if (traits::is_infinite<REALSXP>(l[i]) && traits::is_infinite<REALSXP>(r[i])) {
      z[i] = rnorm(1, mu[i], sd[i])[0];
    } else if (traits::is_infinite<REALSXP>(r[i])) {
      while (u[i] > rho[i]) {
        alpha[i] = (lstd[i] + sqrt(lstd[i]*lstd[i] + 4.0))/2.0;
        z[i] = rexp(1, alpha[i])[0] + lstd[i];
        rho[i] = exp(-(z[i] - alpha[i])*(z[i] - alpha[i])/2.0);
        u[i] = runif(1, 0, 1.0)[0];
      }
      z[i] = z[i]*sd[i] + mu[i];
    } else if (traits::is_infinite<REALSXP>(l[i])) {
      mu[i] = -mu[i];
      lstd[i] = (-r[i] - mu[i])/sd[i];
      while (u[i] > rho[i]) {
        alpha[i] = (lstd[i] + sqrt(lstd[i]*lstd[i] + 4.0))/2.0;
        z[i] = rexp(1, alpha[i])[0] + lstd[i];
        rho[i] = exp(-(z[i] - alpha[i])*(z[i] - alpha[i])/2.0);
        u[i] = runif(1, 0.0, 1.0)[0];
      }
      z[i] = -(z[i]*sd[i] + mu[i]);
    } else {
      while (u[i] > rho[i]) {
        z[i] = runif(1, lstd[i], rstd[i])[0];
        rho[i] = lstd[i] <= 0.0 && rstd[i] >= 0.0 ? exp(-z[i]*z[i]/2) : lstd[i] > 0.0 ? exp((lstd[i]*lstd[i] -z[i]*z[i])/2.0) : exp((rstd[i]*rstd[i] -z[i]*z[i])/2.0);
        u[i] = runif(1, 0.0, 1.0)[0];
      }
      z[i] = z[i]*sd[i] + mu[i];
    }
  }
  return z;
}

// [[Rcpp::export]]
NumericVector rshiftexp(NumericVector d, NumericVector t) {
  
  int p = d.size();
  
  NumericVector z(p, 0.0);
  
  for (int i = 0; i < p; i++) {
    z[i] = rexp(1, d[i])[0] + t[i]; // Counterintuitively, rexp is parametrized
    // in terms of scale = 1/rate
  }
  
  return z;
  
}

// [[Rcpp::export]]
arma::mat remcol(arma::mat A, int i) {
  A.shed_col(i);
  return A;
}

// [[Rcpp::export]]
arma::colvec remrow(arma::colvec a, int i) {
  a.shed_row(i);
  return a;
}

// [[Rcpp::export]]
arma::colvec sampleBeta(NumericVector start, NumericVector Uty,
                        NumericVector delta, NumericVector d, NumericMatrix L, 
                        NumericMatrix Linv, double sigsq, NumericMatrix W) {
  
  int p = start.size();
  int q = delta.size();
  
  // Convert to ARMA objects, trying to minimize memory reallocation:
  // According to: http://dirk.eddelbuettel.com/papers/rcpp_ku_nov2013-part2.pdf
  arma::colvec startAR(start.begin(), start.size(), false);
  arma::colvec UtyAR(Uty.begin(), Uty.size(), false);
  arma::colvec dAR(d.begin(), d.size(), false);
  
  arma::colvec deltaAR =  arma::zeros<arma::colvec>(2*q);
  int k = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < q; j++) {
      deltaAR(k) = delta(j);
      k++;
    }
  }
  
  arma::mat LAR(L.begin(), L.nrow(), L.ncol(), false);
  arma::mat LinvAR(Linv.begin(), Linv.nrow(), Linv.ncol(), false);
  arma::mat WAR(W.begin(), W.nrow(), W.ncol(), false);
  
  arma::colvec z = LAR.t()*startAR;
  
  arma::colvec right = arma::zeros<arma::colvec>(q); arma::colvec left=arma::zeros<arma::colvec>(q);
  NumericVector lowlim(1, 0.0); NumericVector upplim(1, 0.0);
  NumericVector mm(1, 0.0); NumericVector ss(1, 0.0);
  arma::colvec ratio=arma::zeros<arma::colvec>(q);
  arma::colvec ret=arma::zeros<arma::colvec>(p);
  
  for (int i = (p - 1); i >= 0; i--) {
    right = deltaAR - remcol(WAR, i)*remrow(z, i);
    
    left = WAR.col(i);
    
    ratio = right/left;
    
    upplim = arma::min(ratio.elem(find(left > 0.0)));
    lowlim = arma::max(ratio.elem(find(left < 0.0)));
    
    mm = UtyAR.row(i)/dAR.row(i);
    
    ss = sigsq/dAR.row(i);
    
    if (dAR.row(i)[0] > 0) {
      z.row(i) = rtnormrej(mm, sqrt(ss), lowlim, upplim)[0];
      // Rcout << z.row(i) << " " << mm << " " << ss << " " << lowlim << " " << upplim << "\n";
    } else {
      z.row(i) = R::runif(lowlim[0], upplim[0]);
    }
  }
  
  ret = LinvAR.t()*z;
  
  return ret;
  
}

// [[Rcpp::export]]
NumericVector sampleGamma(NumericVector beta, double tausq, double q) {
  
  int p = beta.size();
  
  NumericVector etaq = pow(sqrt(tgamma(3.0/q)/tgamma(1.0/q))*sqrt(2.0/tausq)*abs(beta), q);
  NumericVector rate(p, pow(2.0, -q/2.0));
  
  NumericVector gamma = rshiftexp(rate, etaq);
  
  return gamma;
}

// [[Rcpp::export]]
List sampler(const NumericVector &Uty, const NumericMatrix &L, 
             const NumericMatrix &Linv, const NumericVector &d,
             const NumericMatrix &W, double sigsq, double tausq, double q, 
             const int &samples, NumericVector start, int seed) {
  
  set_seed(seed);
  
  int p = Uty.size();
  
  NumericMatrix resBeta(samples, p);
  NumericMatrix resGamma(samples, p);
  NumericMatrix resVar(samples, 2);
  NumericVector b = start;
  NumericVector g(p, 0.0);
  NumericVector delta(p, 0.0);
  
  // Cat statements in this loop can be used to catch errors
  for (int i = 0; i < samples; i++) {
    g = sampleGamma(b, tausq, q);
    // Rcout << "g: " << g << "\n";
    delta = sqrt(tgamma(1.0/q)/tgamma(3.0/q))*sqrt(tausq/2.0)*pow(g, 1.0/q); // This checks out
    // Rcout << "delta: " << delta << "\n";
    b = sampleBeta(b, Uty, delta, d, L, Linv, sigsq, W);
    // Rcout << "b: " << b << "\n";
    resBeta(i,_) = b;
    resGamma(i,_) = g;
    resVar(i,0) = tausq;
    resVar(i,1) = sigsq;
  }
  
  
  return(Rcpp::List::create(Rcpp::Named("beta")=resBeta,
                            Rcpp::Named("gamma")=resGamma,
                            Rcpp::Named("var")=resVar));
  
}
