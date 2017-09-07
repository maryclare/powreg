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
// Uniform rejection sampler, via Robert (1995)
double rtnormboundunif(double lstd, double rstd) {
  double rho = 0.0;
  double u = 1.0;
  double z = 0.0;
  while (u > rho) {
    z = runif(1, lstd, rstd)[0];
    rho = lstd <= 0.0 && rstd >= 0.0 ? exp(-z*z/2) : lstd > 0.0 ? exp((lstd*lstd -z*z)/2.0) : exp((rstd*rstd -z*z)/2.0);
    u = runif(1, 0.0, 1.0)[0];
  }
  return z;
}

// Lazy normal rejection sampler for truncated normal
double rtnormboundnorm(double lstd, double rstd) {
  double z = INFINITY;
  while (z < lstd | z > rstd) {
    z = rnorm(1, 0.0, 1.1)[0];
  }
  return z;
}

// Rejection sampler for truncated half normal distribution that
// can be used when lower bound exceeds 0
double rtnormboundhalf(double lstd, double rstd) {
  double z = INFINITY;
  while (z < lstd | z > rstd) {
    z = fabs(rnorm(1, 0.0, 1.0)[0]);
  }
  return z;
}

// Shifted exponential rejection sampler for truncated normal
// Based on http://www.stat.ncsu.edu/information/library/papers/mimeo2649_Li.pdf
// This one doesn't quite work
// double rtnormboundtexp(double lstd, double rstd) {
//   double z = INFINITY;
//   double rate = (lstd + sqrt(lstd*lstd + 4.0))/2.0;
//   while (z > rstd) {
//     z = rexp(1, rate)[0] + lstd;
//   }
//   return z;
// }

// Inverse-cdf sampler
// Doesn't save much time  and can be unstable, so not using
double rtnormboundicdf(double lstd, double rstd) {
  
  double plstd = R::pnorm(lstd, 0.0, 1.0, 1, 0);
  double prstd = R::pnorm(rstd, 0.0, 1.0, 1, 0);
  double z = 0.0;
  
  if (plstd != prstd & plstd != 0.0 & plstd != 1.0 & prstd != 0.0 & prstd != 1.0) {
    
    double u = runif(1, 0.0, 1.0)[0];
    z = R::qnorm((prstd - plstd)*u + plstd, 0.0, 1.0, 1, 0);
    
  } else {
    // Rcout << "Unif sampling\n";
    z = rtnormboundunif(lstd, rstd);
  }
  return z;
}

// Optimized sampling of truncated normal for when 
// lower bound exceeds zero based on:
// http://www.stat.ncsu.edu/information/library/papers/mimeo2649_Li.pdf
// Did not actually increase speed
double rtnormpos(double lstd, double rstd) {
  double a0 = 0.2570;
  double b1 = 0.0; 
  double b2 = 0.0;
  double z = 0.0;
  
  if (lstd <= a0) {
    b1 = lstd + sqrt(M_PI/2.0)*exp(lstd*lstd/2.0);
    if (rstd <= b1) {
      z = rtnormboundunif(lstd, rstd);
    } else {
      z = rtnormboundhalf(lstd, rstd);
    }
    
  } else {
    b2 = lstd + (2.0/(lstd + sqrt(lstd*lstd + 4.0)))*exp((lstd*lstd - lstd*sqrt(lstd*lstd + 4.0))/4.0 + 1.0/2.0);
    z = rtnormboundunif(lstd, rstd);
    if (rstd <= b2) {
      z = rtnormboundunif(lstd, rstd);
    } else {
      // z = rtnormboundtexp(lstd, rstd);
      z = rtnormboundunif(lstd, rstd);
    }
  }
  return z;
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
      z[i] = rtnormboundunif(lstd[i], rstd[i])*sd[i] + mu[i];
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
arma::colvec sampleBeta(NumericVector start, NumericVector DUty,
                        NumericVector delta, NumericVector d, NumericMatrix Vt, double sigsq, NumericMatrix W) {
  
  int p = start.size();
  int q = delta.size();
  
  // Convert to ARMA objects, trying to minimize memory reallocation:
  // According to: http://dirk.eddelbuettel.com/papers/rcpp_ku_nov2013-part2.pdf
  arma::colvec startAR(start.begin(), start.size(), false);
  arma::colvec DUtyAR(DUty.begin(), DUty.size(), false);
  arma::colvec dAR(d.begin(), d.size(), false);
  
  arma::colvec deltaAR =  arma::zeros<arma::colvec>(2*q);
  int k = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < q; j++) {
      deltaAR(k) = delta(j);
      k++;
    }
  }
  
  arma::mat VtAR(Vt.begin(), Vt.nrow(), Vt.ncol(), false);
  arma::mat WAR(W.begin(), W.nrow(), W.ncol(), false);
  
  arma::colvec z = VtAR*startAR;
  
  arma::colvec right = arma::zeros<arma::colvec>(q); arma::colvec left=arma::zeros<arma::colvec>(q);
  NumericVector lowlim(1, 0.0); NumericVector upplim(1, 0.0);
  NumericVector mm(1, 0.0); NumericVector ss(1, 0.0);
  arma::colvec ratio=arma::zeros<arma::colvec>(q);
  arma::colvec ret=arma::zeros<arma::colvec>(p);
  
  for (int i = 0; i < p; i++) {
    
    right = deltaAR - remcol(WAR, i)*remrow(z, i);
    
    left = WAR.col(i);
    
    ratio = right/left;
    
    upplim = arma::min(ratio.elem(find(left > 0.0)));
    lowlim = arma::max(ratio.elem(find(left < 0.0)));
    
    mm = DUtyAR.row(i)/(dAR.row(i)*dAR.row(i));
    
    ss = sigsq/((dAR.row(i)*dAR.row(i)));
    
    if (dAR.row(i)[0] > 0) {
      z.row(i) = rtnormrej(mm, sqrt(ss), lowlim, upplim)[0];
      // Rcout << z.row(i) << " " << mm << " " << ss << " " << lowlim << " " << upplim << "\n";
    } else if (dAR.row(i)[0] == 0) {
      z.row(i) = R::runif(lowlim[0], upplim[0]);
    }
  }
  
  ret = VtAR.t()*z;
  
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
List sampler(const NumericVector &DUty, const NumericMatrix &Vt, const NumericVector &d,
             const NumericMatrix &W, double sigsq, double tausq, double q, 
             const int &samples, NumericVector start, int seed, const int &burn, 
             const int &thin) {
  
  set_seed(seed);
  
  int p = DUty.size();
  
  NumericMatrix resBeta(samples, p);
  NumericMatrix resGamma(samples, p);
  NumericVector b = start;
  NumericVector g(p, 0.0);
  NumericVector delta(p, 0.0);
  
  // Cat statements in this loop can be used to catch errors
  for (int i = 0; i < (samples*thin + burn); i++) {
    // Might be better to sample from inverse gamma but rejection sampler for that is not obvious
    g = sampleGamma(b, tausq, q);
    // Rcout << "g: " << g << "\n";
    delta = sqrt(tgamma(1.0/q)/tgamma(3.0/q))*sqrt(tausq/2.0)*pow(g, 1.0/q); // This checks out
    // Rcout << "delta: " << delta << "\n";
    b = sampleBeta(b, DUty, delta, d, Vt, sigsq, W);
    // Rcout << "b: " << b << "\n";
    if (i >= burn & i % thin == 0) {
      resBeta((i - burn)/thin,_) = b;
      resGamma((i - burn)/thin,_) = g;
    }
  }
  
  
  return(Rcpp::List::create(Rcpp::Named("beta")=resBeta,
                            Rcpp::Named("gamma")=resGamma));
  
}
