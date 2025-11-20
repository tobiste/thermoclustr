#include <Rcpp.h>
using namespace Rcpp;

// -------------------------------------------------------------
// Euclidean distance matrix between P and Q
// -------------------------------------------------------------
NumericMatrix distmat_cpp(const NumericMatrix& P,
                          const NumericMatrix& Q) {
  
  int m = P.nrow();
  int n = Q.nrow();
  int k = P.ncol();
  
  NumericMatrix D(m, n);
  
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      double sum = 0.0;
      for (int c = 0; c < k; c++) {
        double d = P(i, c) - Q(j, c);
        sum += d * d;
      }
      D(i, j) = std::sqrt(sum);
    }
  }
  
  return D;
}

// -------------------------------------------------------------
// Directed Hausdorff distance h(P, Q)
// max over points in P of the min distance to Q
// -------------------------------------------------------------
double directed_hausdorff(const NumericMatrix& P,
                          const NumericMatrix& Q) {
  
  int m = P.nrow();
  int n = Q.nrow();
  int k = P.ncol();
  
  double max_min_dist = 0.0;
  
  for (int i = 0; i < m; i++) {
    double min_dist = R_PosInf;
    
    for (int j = 0; j < n; j++) {
      double sum = 0.0;
      for (int c = 0; c < k; c++) {
        double d = P(i, c) - Q(j, c);
        sum += d * d;
      }
      double dist = std::sqrt(sum);
      if (dist < min_dist)
        min_dist = dist;
    }
    
    if (min_dist > max_min_dist)
      max_min_dist = min_dist;
  }
  return max_min_dist;
}

// -------------------------------------------------------------
// Full symmetric Hausdorff
// -------------------------------------------------------------
// [[Rcpp::export]]
double hausdorff_cpp(const NumericMatrix& P,
                     const NumericMatrix& Q) {
  
  int m = P.nrow();
  int n = Q.nrow();
  
  if (m == 0 && n == 0)
    return 0.0;
  if (m == 0 || n == 0)
    return NA_REAL;   // or R_PosInf if preferred
  
  double hPQ = directed_hausdorff(P, Q);
  double hQP = directed_hausdorff(Q, P);
  
  return (hPQ > hQP ? hPQ : hQP);
}

// -------------------------------------------------------------
// Pairwise Hausdorff Distance Matrix (lower triangle only)
// -------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix pairwise_hausdorff_cpp(List mats) {
  
  int n = mats.size();
  NumericMatrix out(n, n);
  
  for (int i = 0; i < n; i++) {
    NumericMatrix Pi = mats[i];
    for (int j = 0; j < i; j++) {
      NumericMatrix Pj = mats[j];
      out(i, j) = hausdorff_cpp(Pi, Pj);
    }
  }
  
  return out;
}
