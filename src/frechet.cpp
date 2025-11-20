#include <Rcpp.h>
using namespace Rcpp;

// Euclidean distance between 2 points (k-dimensional)
inline double point_dist(const NumericVector& a,
                         const NumericVector& b) {
  int k = a.size();
  double sum = 0.0;
  for (int i = 0; i < k; i++) {
    double d = a[i] - b[i];
    sum += d * d;
  }
  return std::sqrt(sum);
}

// -------------------------------------------------------------
// Continuous Frechet distance between polygonal curves P and Q
// P: m x k matrix
// Q: n x k matrix
// -------------------------------------------------------------
// [[Rcpp::export]]
double cont_frechet_cpp(const NumericMatrix& P,
                   const NumericMatrix& Q) {
  
  int m = P.nrow();
  int n = Q.nrow();
  
  // empty handling (optional)
  if (m == 0 && n == 0)
    return 0.0;
  if (m == 0 || n == 0)
    return NA_REAL;
  
  // DP matrix
  NumericMatrix ca(m, n);
  
  // initialize with sentinel
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      ca(i, j) = -1.0;
  
  // recursive function (memoized)
  std::function<double(int,int)> c = [&](int i, int j) {
    
    if (ca(i,j) > -1.0)
      return ca(i,j);
    
    double d = point_dist(P(i, _), Q(j, _));
    
    if (i == 0 && j == 0) {
      ca(i,j) = d;
    } else if (i > 0 && j == 0) {
      ca(i,j) = std::max(c(i-1, 0), d);
    } else if (i == 0 && j > 0) {
      ca(i,j) = std::max(c(0, j-1), d);
    } else {
      double prev = std::min(
        std::min(c(i-1, j), c(i-1, j-1)),
        c(i, j-1)
      );
      ca(i,j) = std::max(prev, d);
    }
    return ca(i,j);
  };
  
  return c(m-1, n-1);
}



// [[Rcpp::export]]
double discrete_frechet_cpp(NumericMatrix A, NumericMatrix B) {
  
  int m = A.nrow();
  int n = B.nrow();
  
  if (m == 0 || n == 0)
    return NA_REAL;  // you said empties are removed already
  
  int dA = A.ncol();
  int dB = B.ncol();
  if (dA != dB)
    stop("A and B must have same number of columns");
  
  //int d = dA;
  
  // DP matrix
  NumericMatrix ca(m, n);
  
  // (0,0)
  {
    NumericVector a0 = A(0, _);
    NumericVector b0 = B(0, _);
    ca(0,0) = point_dist(a0, b0);
  }
  
  // First row
  for (int j = 1; j < n; j++) {
    NumericVector a0 = A(0, _);
    NumericVector bj = B(j, _);
    double dist = point_dist(a0, bj);
    ca(0,j) = std::max(ca(0,j-1), dist);
  }
  
  // First column
  for (int i = 1; i < m; i++) {
    NumericVector ai = A(i, _);
    NumericVector b0 = B(0, _);
    double dist = point_dist(ai, b0);
    ca(i,0) = std::max(ca(i-1,0), dist);
  }
  
  // Fill DP
  for (int i = 1; i < m; i++) {
    NumericVector ai = A(i, _);
    for (int j = 1; j < n; j++) {
      NumericVector bj = B(j, _);
      double dist = point_dist(ai, bj);
      
      double prev_min = std::min(
        ca(i-1, j),
        std::min(ca(i-1, j-1), ca(i, j-1))
      );
      
      ca(i,j) = std::max(dist, prev_min);
    }
  }
  
  return ca(m-1, n-1);
}




// -------------------------------------------------------------
// Pairwise Frechet distance matrix (lower triangle)
// -------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix pairwise_cont_frechet_cpp(List mats) {
  
  int n = mats.size();
  NumericMatrix out(n, n);
  
  for (int i = 0; i < n; i++) {
    NumericMatrix Pi = mats[i];
    for (int j = 0; j < i; j++) {
      NumericMatrix Pj = mats[j];
      out(i, j) = cont_frechet_cpp(Pi, Pj);
    }
  }
  
  return out;
}

// -------------------------------------------------------------
// Pairwise Frechet distance matrix (lower triangle)
// -------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix pairwise_discr_frechet_cpp(List mats) {
  
  int n = mats.size();
  NumericMatrix out(n, n);
  
  for (int i = 0; i < n; i++) {
    NumericMatrix Pi = mats[i];
    for (int j = 0; j < i; j++) {
      NumericMatrix Pj = mats[j];
      out(i, j) = discrete_frechet_cpp(Pi, Pj);
    }
  }
  
  return out;
}
