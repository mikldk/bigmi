#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

#include <unordered_map>
#include <cmath>

#include "helper_hash.h"

std::unordered_map<int, double> 
calc_map_p_x(const std::vector<int>& xs) {
  
  const int n = xs.size();
  
  const double one_nth = 1.0 / (double)n;
  std::unordered_map<int, double> map_p_x;
  
  for (int i = 0; i < n; ++i) {
    const int x = xs[i];
    map_p_x[x] += one_nth;
  }
  
  return map_p_x;
}

std::unordered_map< std::pair<int, int> , double, pairhash> 
  calc_map_p_xy(const std::vector<int>& xs, const std::vector<int>& ys) {
    
    const int n = xs.size();
    
    if (ys.size() != n) {
      Rcpp::stop("Unexpected error: xs and ys must have same size");
    }
    
    const double one_nth = 1.0 / (double)n;
    std::unordered_map< std::pair<int, int> , double, pairhash> map_p_xy;
    
    for (int i = 0; i < n; ++i) {
      const std::pair<int, int> xy = std::make_pair(xs[i], ys[i]);
      map_p_xy[xy] += one_nth;
    }
    
    return map_p_xy;
  }

double MI_core(std::unordered_map<int, double> map_p_x,
               std::unordered_map<int, double> map_p_y,
               std::unordered_map< std::pair<int, int> , double, pairhash> map_p_xy) {
  
  double MI = 0.0;
  
  // Loop over x and y:
  for (auto x_ele = map_p_x.begin(); x_ele != map_p_x.end(); ++x_ele) { 
    const int x = x_ele->first;
    const double p_x = map_p_x[x];
    
    for (auto y_ele = map_p_y.begin(); y_ele != map_p_y.end(); ++y_ele) { 
      const int y = y_ele->first;
      const double p_y = map_p_y[y];
      
      std::pair<int, int> xy = std::make_pair(x, y);
      const double p_xy = map_p_xy[xy];
      
      // p_xy = 0 => 0*log(0/c) = 0:
      if (p_xy < MACHINE_EPS) {
        continue; // FIXME: What about just using counts? And checking if they are zero?
      }
      
      const double contrib = p_xy * log(p_xy / (p_x * p_y));
      
      MI += contrib;
    }
  } 
  
  return MI;
}

// [[Rcpp::export]]
double MI_categorical_worker_two(const Rcpp::IntegerMatrix& d) {
  if (d.ncol() != 2) {
    Rcpp::stop("Expected d to have exactly two columns");
  }
  
  // Transpose x to obtain faster access: column access is faster than row access
  const Rcpp::IntegerMatrix dT = Rcpp::transpose(d);
  
  const Rcpp::IntegerVector xs_intvec = d(Rcpp::_, 0);
  const Rcpp::IntegerVector ys_intvec = d(Rcpp::_, 1);
  const std::vector<int> xs = Rcpp::as< std::vector<int> >(xs_intvec);
  const std::vector<int> ys = Rcpp::as< std::vector<int> >(ys_intvec);
  
  std::unordered_map<int, double> map_p_x = calc_map_p_x(xs);
  std::unordered_map<int, double> map_p_y = calc_map_p_x(ys);
  std::unordered_map< std::pair<int, int> , double, pairhash> map_p_xy = 
    calc_map_p_xy(xs, ys);
  
  double mi = MI_core(map_p_x, map_p_y, map_p_xy);
  
  return mi;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix MI_categorical_worker_all(const Rcpp::IntegerMatrix& d) {
  // Transpose x to obtain faster access: column access is faster than row access
  const Rcpp::IntegerMatrix dT = Rcpp::transpose(d);
  
  const int vars = d.ncol();
  
  if (vars <= 1) {
    Rcpp::stop("d must have at least 2 columns");
  }
  
  Rcpp::NumericMatrix MIs(vars, vars);
  Rcpp::colnames(MIs) = Rcpp::colnames(d);
  Rcpp::rownames(MIs) = Rcpp::colnames(d);
  
  // Normally only to i < (vars - 1), see below
  for (int i = 0; i < vars; ++i) {
    const Rcpp::IntegerVector xs_intvec = d(Rcpp::_, i);
    const std::vector<int> xs = Rcpp::as< std::vector<int> >(xs_intvec);
    std::unordered_map<int, double> map_p_x = calc_map_p_x(xs);
    
    std::unordered_map< std::pair<int, int> , double, pairhash> map_p_xx = 
      calc_map_p_xy(xs, xs);
    
    double mi = MI_core(map_p_x, map_p_x, map_p_xx);
    MIs(i, i) = mi;
    
    /*
    * Normally, loop over i would be
    *   for (int i = 0; i < (vars - 1); ++i) {
    * but we need MI for the last one as well, obatained this way.
    */
    if (i == (vars - 1)) {
      break;
    }
    
    for (int j = i+1; j < vars; ++j) {
      const Rcpp::IntegerVector ys_intvec = d(Rcpp::_, j);
      const std::vector<int> ys = Rcpp::as< std::vector<int> >(ys_intvec);
      std::unordered_map<int, double> map_p_y = calc_map_p_x(ys);
      
      std::unordered_map< std::pair<int, int> , double, pairhash> map_p_xy = 
        calc_map_p_xy(xs, ys);
      
      double mi = MI_core(map_p_x, map_p_y, map_p_xy);
      MIs(i, j) = mi;
      MIs(j, i) = mi;
    }
    }
  
  return MIs;
}


//' Sort mutual information matrix
//' 
//' Sort square matrix with mutual information, e.g. obtained by [MI_categorical_all_pairwise()].
//' 
//' @param d Matrix with observations
//' @return Matrix with two columns with the (row, column) indices
//' 
//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix MI_order(const Rcpp::IntegerMatrix& d) {
  int n = d.nrow();
  
  if (d.ncol() != n) {
    Rcpp::stop("Expected square matrix");
  }
  
  std::vector< std::pair<int, int> > pairs;
  
  for (int i = 0; i < (n-1); ++i) {
    // upper triangular, excluding diagonal
    for (int j = (i+1); j < n; ++j) {
      std::pair<int, int> xy = std::make_pair(i, j);
      pairs.push_back(xy);
    }
  }
  
  std::sort(pairs.begin(), pairs.end(), [=](std::pair<int, int>& a, std::pair<int, int>& b) {
    // Decreasing
    return d(a.first, a.second) > d(b.first, b.second);
    //return d(a.first, a.second) <= d(b.first, b.second);
  });
  
  int N = n*(n-1)/2;
  
  if (N != pairs.size()) {
    Rcpp::stop("Unexpected size");
  }
  
  int i = 0;
  Rcpp::IntegerMatrix ord(N, 2);
  
  for (auto e = pairs.begin(); e != pairs.end(); ++e) {
    ord(i, 0) = e->first + 1; // R's 1-based index compared to C++'s 0-based
    ord(i, 1) = e->second + 1;
    ++i;
  }
  
  return ord;
}
