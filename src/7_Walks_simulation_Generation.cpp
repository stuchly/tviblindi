#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <unordered_map>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#define NUM_THREADS(N) ((N) >= 0 ? (N) : omp_get_num_procs() + (N) + 1)

using namespace Rcpp;

// [[Rcpp::export]]
RcppExport SEXP random_walks_from_adj(const Eigen::Map<Eigen::SparseMatrix<double> > A,
                                      int                                            start,
                                      const Eigen::Index                             N = 1) {
  Rcpp::IntegerVector v;
  Rcpp::IntegerVector ind_s(N);
  Rcpp::NumericVector gprobs(N);
  const size_t nb_iter = (size_t) - 1;
  int          origin  = start;
  int          jj      = 0;
  
  for (int j = 0; j < N; ++j)
  {
    ind_s[j] = jj + 1;  // 1-based start index of random walk
    v.push_back(origin);
    start = origin - 1; // 0-based index of start node
    int    path_length = 0;
    double logpiprob   = 0; // log likelihoods are cumulative
    for (int i = 1; i <= nb_iter; ++i)
    {
      
      Rcpp::NumericVector probs;
      Rcpp::IntegerVector vert, inds_vert;
      
      int j_index = 0;
      
      for (Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator it(A, start); it; ++it)
      {
        probs.push_back(it.value());
        vert.push_back(it.index());
        inds_vert.push_back(j_index);
        ++j_index;
      }
      if (probs.length() == 0) break; // terminal node reached
      Rcpp::IntegerVector ret = sample(inds_vert, 1, FALSE, probs); // next node sampled using transition probabilities
      ++jj;
      logpiprob += log(probs[ret[0]]); // increase cumulative likelihood
      v.push_back(vert[ret[0]] + 1);
      start = vert[ret[0]];
      ++path_length;
    }
    gprobs[j] = exp(logpiprob / path_length); // likelihood associated with this random walk
    ++jj;
  }
  return Rcpp::List::create(
    Rcpp::Named("v",      v[v>0]),
    Rcpp::Named("starts", ind_s),
    Rcpp::Named("gprobs", gprobs)
  );
}

// [[Rcpp::export]]
RcppExport SEXP  C_random_walk_adj(const Eigen::Map<Eigen::SparseMatrix<double> > A,
                                   int start,
                                   const Eigen::Index nb_iter=100) {
  Rcpp::IntegerVector v(nb_iter+1);
  v[0]=start;
  start--;
  
  for (int i=1; i<=nb_iter; i++){
    Rcpp::NumericVector probs;
    Rcpp::IntegerVector vert;
    for (Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator it(A,start); it; ++it){
      probs.push_back(it.value());
      vert.push_back(it.index());
    }
    if (probs.length()==0) break;
    Rcpp::IntegerVector ret = sample(vert, 1, FALSE, probs) ;
    v[i]=ret[0]+1;
    start=ret[0];
  }
  
  return Rcpp::List::create(Rcpp::Named("v", v[v>0]));
}

// [[Rcpp::export]]
RcppExport SEXP  C_random_walk_adj_N(const Eigen::Map<Eigen::SparseMatrix<double> > A,
                                     int start,
                                     const Eigen::Index nb_iter=100,
                                     const Eigen::Index N=1) {
  Rcpp::IntegerVector v((nb_iter+1)*N);
  Rcpp::IntegerVector ind_s(N);
  int origin=start;
  int jj=0;
  for (int j=0;j<N;j++){
    ind_s[j]=jj+1;
    v[jj]=origin;
    start=origin-1;
    for (int i=1; i<=nb_iter; i++){
      Rcpp::NumericVector probs;
      Rcpp::IntegerVector vert;
      for (Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator it(A,start); it; ++it){
        probs.push_back(it.value());
        vert.push_back(it.index());
      }
      
      // int pp=probs.length();
      // std::cout << pp <<" " << i << std::endl;
      
      if (probs.length()==0) break;
      Rcpp::IntegerVector ret = sample(vert, 1, FALSE, probs) ;
      jj++;
      v[jj]=ret[0]+1;
      start=ret[0];
    }
    jj++;
  }
  return Rcpp::List::create(Rcpp::Named("v", v[v>0]),Rcpp::Named("starts",ind_s));
}

// [[Rcpp::export]]
RcppExport SEXP  C_random_walk_adj_N_push(const Eigen::Map<Eigen::SparseMatrix<double> > A,
                                          int start,
                                          const Eigen::Index N=1) {
  Rcpp::IntegerVector v;
  Rcpp::IntegerVector ind_s(N);
  Rcpp::NumericVector gprobs(N);
  const size_t nb_iter = (size_t)-1;
  int origin = start;
  int jj = 0;
  
  for (int j=0; j<N; j++) {
    
    ind_s[j] = jj+1;
    v.push_back(origin);
    start = origin-1;
    int path_length=0;
    double logpiprob=0;
    for (int i=1; i<=nb_iter; i++) {
      
      Rcpp::NumericVector probs;
      Rcpp::IntegerVector vert,inds_vert;
     
      int j_index=0;
      
      for (Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator it(A,start); it; ++it) {
        
        probs.push_back(it.value());
        vert.push_back(it.index());
        inds_vert.push_back(j_index);
        j_index++;
      }
      if (probs.length()==0) break;
      Rcpp::IntegerVector ret = sample(inds_vert, 1, FALSE, probs);
      jj++;
      logpiprob+=log(probs[ret[0]]);
      v.push_back(vert[ret[0]]+1);
      start = vert[ret[0]];
      path_length++;
      
    }
    gprobs[j]=exp(logpiprob/path_length);
    jj++;
    
  }
  return Rcpp::List::create(Rcpp::Named("v", v[v>0]), Rcpp::Named("starts",ind_s), Rcpp::Named("gprobs",gprobs));
}

namespace remove_cycles {
  template <typename T>
    std::vector<T> remove_cycles_generic(
      std::vector<T> series,
      std::vector<T> series_unique,
      int& cyclical_count,
      const T placeholder,
      int series_number = -1
    ) {
      T last_point = series.back();
      
      std::unordered_map<T, std::tuple<int, int>> intervals;
      
      for (auto i = series_unique.begin(); i != series_unique.end(); ++i)
      {
        intervals.insert({*i, std::make_tuple(-1, -1)});
      }
      
      for (int i = 0; i < series.size(); ++i)
      {
        auto it = intervals.find(series[i]);
        if (std::get<0>(it->second) == -1)
        {
          it->second = std::make_tuple(i+1, i+1);
        } else {
          std::get<1>(it->second) = i;
        }
      }
      
      bool any_cycles = false;
      
      for (auto i = intervals.begin(); i != intervals.end(); ++i)
      {
        int first = std::get<0>(i->second);
        int last = std::get<1>(i->second);
        if (last != first)
        {
          any_cycles = true;
          for (int j = first; j <= last; ++j)
          {
            series[j] = placeholder;
          }
        }
      }
      
      series.erase(std::remove(series.begin(), series.end(), placeholder), series.end());
      
      if (series.back() != last_point)
      {
        series.push_back(last_point);
      }
      
      return series;
    }
  
  std::vector<int> remove_cycles_int(
    std::vector<int> series,
    std::vector<int> series_unique,
    int& cyclical_count,
    int series_number = -1
  ) {
    return remove_cycles::remove_cycles_generic<int>(series, series_unique, cyclical_count, -1, series_number);
  }
}


// [[Rcpp::export]]
std::vector< std::vector<int> > remove_cycles_int_list(
  std::vector< std::vector<int> > series_list,
  std::vector< std::vector<int> > series_list_unique
) {
  int cyclical_count = 0;
  int N = series_list.size();
  
  #pragma omp parallel for
  for (int i = 0; i < N; ++i)
  {
    series_list[i] = remove_cycles::remove_cycles_int(series_list[i], series_list_unique[i], cyclical_count, i+1);
  }
  
  return series_list;
}
