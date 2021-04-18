#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector< std::vector<double> > interpolate_trajectories(
  std::vector< std::vector<double> > pts,
  double          coef = 200
) {
  std::vector< std::vector<double> > out;
  
  std::vector<double> first = pts[0];
  std::vector<double> second;
  int dim = first.size();
  for (size_t idx = 1; idx < pts.size() - 1; ++idx)
  {
    second = pts[idx];
    std::vector<double> vec;
    double dist = .0;
    for (size_t d = 0; d < dim; ++d)
    {
      double tmp = second[d] - first[d];
      vec.push_back(tmp);
      dist += tmp * tmp;
    }
    dist = sqrt(dist);
    out.push_back(first);
    int npt = std::round(coef * dist);
    for (size_t j = 0; j < dim; ++j)
    {
      vec[j] /= npt;
    }
    for (int i = 0; i < npt; ++i)
    {
      for (size_t j = 0; j < dim; ++j)
      {
        first[j] += vec[j];
      }
      out.push_back(first);
    }
    first = second;
  }
  out.push_back(second);
  return out;
}