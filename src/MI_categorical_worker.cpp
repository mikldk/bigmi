//#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

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


// https://stackoverflow.com/a/17074810
template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
    const std::vector<T>& vec,
    Compare compare)
{
  std::vector<std::size_t> p(vec.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),
            [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
  return p;
}

template <typename T>
void apply_permutation_in_place(
    std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
  std::vector<bool> done(vec.size());
  for (std::size_t i = 0; i < vec.size(); ++i)
  {
    if (done[i])
    {
      continue;
    }
    done[i] = true;
    std::size_t prev_j = i;
    std::size_t j = p[i];
    while (i != j)
    {
      std::swap(vec[prev_j], vec[j]);
      done[j] = true;
      prev_j = j;
      j = p[j];
    }
  }
}


// [[Rcpp::export]]
Rcpp::List MI_categorical_worker_sparse_all(const Rcpp::IntegerMatrix& d, 
                                            const bool progress = true) {
  const int vars = d.ncol();
  
  if (vars <= 1) {
    Rcpp::stop("d must have at least 2 columns");
  }
  
  size_t p = d.ncol();
  size_t N = p*(p-1)/2;
  
  std::vector<double> mi(N);
  std::vector< std::pair<int, int> > indices(N);
  
  Progress prog(p-1, progress);
  
  size_t idx = 0;
  
  // Filling marginals
  std::unordered_map< size_t, std::unordered_map<int, double> > map_marginals;
  for (size_t i = 0; i < p; ++i) {
    const Rcpp::IntegerVector xs_intvec = d(Rcpp::_, i);
    const std::vector<int> xs = Rcpp::as< std::vector<int> >(xs_intvec);
    
    std::unordered_map<int, double> map_p_x = calc_map_p_x(xs);
    map_marginals[i] = map_p_x;
  }

  for (size_t i = 0; i < (p-1); ++i) {
    if (Progress::check_abort()) {
      Rcpp::List empty;
      empty["error"] = true;
      return empty;
    }
    prog.increment();
    
    const Rcpp::IntegerVector xs_intvec = d(Rcpp::_, i);
    const std::vector<int> xs = Rcpp::as< std::vector<int> >(xs_intvec);
    auto x_marginal_search = map_marginals.find(i);
    std::unordered_map<int, double> map_p_x = x_marginal_search->second;
    
    for (size_t j = i+1; j < p; ++j) {
      const Rcpp::IntegerVector ys_intvec = d(Rcpp::_, j);
      const std::vector<int> ys = Rcpp::as< std::vector<int> >(ys_intvec);
      auto y_marginal_search = map_marginals.find(j);
      std::unordered_map<int, double> map_p_y = y_marginal_search->second;
      
      std::unordered_map< std::pair<int, int> , double, pairhash> map_p_xx = 
        calc_map_p_xy(xs, ys);
      
      double mi_xy = MI_core(map_p_x, map_p_y, map_p_xx);
      
      mi[idx] = mi_xy;
      indices[idx] = std::make_pair<int, int>(i+1, j+1); // R indexing
      ++idx;
    }
  }
  //Rcpp::Rcout << "N     = " << N << std::endl;
  //Rcpp::Rcout << "emp_N = " << idx << std::endl;

  auto permut = sort_permutation(mi, std::greater<double>());
  
  apply_permutation_in_place(mi, permut);
  apply_permutation_in_place(indices, permut);
  
  Rcpp::IntegerMatrix idx_mat(N, 2);
  for (size_t i = 0; i < N; ++i) {
    idx_mat(i, 0) = indices[i].first;
    idx_mat(i, 1) = indices[i].second;
  }
  
  Rcpp::List ret;
  ret["error"] = false;
  ret["mi"] = mi;
  ret["idx"] = idx_mat;
  return ret;
}

// https://stackoverflow.com/a/35833470
double avg(std::vector<int> const& v) {
  return 1.0 * std::accumulate(v.begin(), v.end(), 0LL) / v.size();
}
//https://stackoverflow.com/q/9599552
struct Mean {
  size_t n;
  Mean(size_t n) : n(n) {}
  double operator()(double sum, double x) {
    return sum + x/n;
  }
};


// [[Rcpp::export]]
Rcpp::List pearson_correlation_sparse_all(const Rcpp::NumericMatrix& d, 
                                          const bool progress = true) {
  const int vars = d.ncol();
  
  if (vars <= 1) {
    Rcpp::stop("d must have at least 2 columns");
  }
  
  size_t p = d.ncol();
  size_t n = d.nrow();
  double ndbl = (double)n;
  size_t N = p*(p-1)/2;
  
  std::vector<double> corr(N);
  std::vector< std::pair<int, int> > indices(N);
  
  Progress prog(p-1, progress);
  
  size_t idx = 0;
  
  Mean mean(n);
  
  // centered
  std::vector< std::vector<double> > dat(p);
  for (size_t i = 0; i < p; ++i) {
    const Rcpp::NumericVector xs_vec = d(Rcpp::_, i);
    const std::vector<double> xs = Rcpp::as< std::vector<double> >(xs_vec);
    
    if (xs.size() != n) {
      Rcpp::stop("Unexpected");
    }
    
    double mean_x = std::accumulate(xs.begin(), xs.end(), 0.0, mean);
    
    std::vector<double> x_centered(n);
    for (size_t j = 0; j < n; ++j) {
      x_centered[j] = xs[j] - mean_x;
    }
    //Rcpp::print(Rcpp::wrap(x_centered));
    dat[i] = x_centered;
  }
  
  //stds
  std::unordered_map< size_t, double > map_std_x;
  for (size_t i = 0; i < p; ++i) {
    const std::vector<double> x_centered = dat[i];
    
    double std_x = 0.0;
    for (auto x : x_centered) {
      std_x += x*x / ndbl; // numeric stability vs seed?
    }
    //Rcpp::print(Rcpp::wrap(std_x));
    
    map_std_x[i] = sqrt(std_x);
  }
  
  // cors
  for (size_t i = 0; i < (p-1); ++i) {
    if (Progress::check_abort()) {
      Rcpp::List empty;
      empty["error"] = true;
      return empty;
    }
    prog.increment();
    
    const std::vector<double> x_centered = dat[i];
    auto map_std_x_search = map_std_x.find(i);
    double sd_x = map_std_x_search->second;

    for (size_t j = i+1; j < p; ++j) {
      const std::vector<double> y_centered = dat[j];
      auto map_std_y_search = map_std_x.find(j);
      double sd_y = map_std_y_search->second;
      
      double cov = 0.0;
      
      for (size_t k = 0; k < n; ++k) {
        cov = cov + x_centered[k]*y_centered[k]/ndbl;
      }
      
      const double abs_cor = std::abs(cov)/(sd_x * sd_y);
      //Rcpp::Rcout << "(" << i << ", " << j << ") = (" << sd_x << ", " << sd_y << ") = " << cov << " = " << abs_cor << std::endl;
      
      corr[idx] = abs_cor;
      indices[idx] = std::make_pair<int, int>(i+1, j+1); // R indexing
      ++idx;
    }
  }

  auto permut = sort_permutation(corr, std::greater<double>());
  
  apply_permutation_in_place(corr, permut);
  apply_permutation_in_place(indices, permut);
  
  Rcpp::IntegerMatrix idx_mat(N, 2);
  for (size_t i = 0; i < N; ++i) {
    idx_mat(i, 0) = indices[i].first;
    idx_mat(i, 1) = indices[i].second;
  }
  
  Rcpp::List ret;
  ret["error"] = false;
  ret["corr"] = corr;
  ret["idx"] = idx_mat;
  return ret;
}



