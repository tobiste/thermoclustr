#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double hausdorff_rcpp(NumericMatrix P, NumericMatrix Q) {
  int nP = P.nrow();
  int nQ = Q.nrow();
  int dim = P.ncol();
  double maxP = 0.0;
  double maxQ = 0.0;
  
  // Helper lambda for Euclidean distance
  auto dist = [&](NumericVector a, NumericVector b) {
    double sum = 0.0;
    for (int d = 0; d < dim; ++d) {
      double diff = a[d] - b[d];
      sum += diff * diff;
    }
    return sqrt(sum);
  };
  
  // max min dist P -> Q
  for (int i = 0; i < nP; ++i) {
    double minDist = R_PosInf;
    for (int j = 0; j < nQ; ++j) {
      double d = dist(P(i, _), Q(j, _));
      if (d < minDist) minDist = d;
    }
    if (minDist > maxP) maxP = minDist;
  }
  
  // max min dist Q -> P
  for (int i = 0; i < nQ; ++i) {
    double minDist = R_PosInf;
    for (int j = 0; j < nP; ++j) {
      double d = dist(Q(i, _), P(j, _));
      if (d < minDist) minDist = d;
    }
    if (minDist > maxQ) maxQ = minDist;
  }
  
  return std::max(maxP, maxQ);
}

// [[Rcpp::export]]
NumericMatrix hausdorff_matrix_rcpp(List x_list) {
  int n = x_list.size();
  NumericMatrix out(n, n);
  
  for (int i = 0; i < n; i++) {
    NumericMatrix Pi = as<NumericMatrix>(x_list[i]);
    for (int j = 0; j < i; j++) {  // only lower triangle
      NumericMatrix Pj = as<NumericMatrix>(x_list[j]);
      double d = hausdorff_rcpp(Pi, Pj);
      out(i, j) = d;
      out(j, i) = d;
    }
  }
  return out;
}
