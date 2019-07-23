/*
  remove_cycles_int_list:
    remove cycles from each series of points (as int) in a list
    parameters:
      series_list: list of integer vectors (series)
      series_list_unique: list of unique values of each series
      verbose: should we print which series contain cycles?
*/

#include <Rcpp.h>
#include <unordered_map>
#include <omp.h>
// [[Rcpp::plugins(openmp)]]

#define NUM_THREADS(N) ((N) >= 0 ? (N) : omp_get_num_procs() + (N) + 1)

using namespace Rcpp;

namespace remove_cycles {
  template <typename T>
  std::vector<T> remove_cycles_generic(
      std::vector<T> series,
      std::vector<T> series_unique,
      int& cyclical_count,
      const T placeholder,
      bool verbose = false,
      int series_number = -1
  ) {
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
        if (verbose && !any_cycles)
        {
          Rcpp::Rcout << "Series " << series_number << " cyclical." << std::endl;
          ++cyclical_count;
        }
        any_cycles = true;
        for (int j = first; j <= last; ++j)
        {
          series[j] = placeholder;
        }
      }
    }
    
    series.erase(std::remove(series.begin(), series.end(), placeholder), series.end());
    
    return series;
  }
  
  std::vector<int> remove_cycles_int(
      std::vector<int> series,
      std::vector<int> series_unique,
      int& cyclical_count,
      bool verbose = false,
      int series_number = -1
  ) {
    return remove_cycles::remove_cycles_generic<int>(series, series_unique, cyclical_count, -1, verbose, series_number);
  }
}


// [[Rcpp::export]]
std::vector< std::vector<int> > remove_cycles_int_list(
  std::vector< std::vector<int> > series_list,
  std::vector< std::vector<int> > series_list_unique,
  bool verbose = true
) {
  int cyclical_count = 0;
  int N = series_list.size();
# pragma omp parallel for
  for (int i = 0; i < N; ++i)
  {
    series_list[i] = remove_cycles::remove_cycles_int(series_list[i], series_list_unique[i], cyclical_count, verbose, i+1);
  }
  
  if (verbose)
  {
    Rcpp::Rcout << cyclical_count << " of " << N << " series were cyclical." << std::endl;
  }
  
  return series_list;
}
  
