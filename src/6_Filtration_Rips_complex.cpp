#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <graph_simplicial_complex_adj.h>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for std::numeric_limits

#include <math.h>
#include <stdlib.h>

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
SEXP C_sigma(const NumericMatrix D,
             const double        target,
             const int           nIter = 64,
             const double        tol = 1e-5)
{
  // Find edges pointing back in time
  int                  nrow = D.nrow();
  int                  ncol = D.ncol();
  std::vector<double > sigma(nrow, 1);
  for (int i = 0; i < nrow; ++i)
  {
    double lo  = 0.0;
    double hi  = std::numeric_limits<double>::infinity();
    double mid = 1.0;
    
    for (int j = 0; j < nIter; ++j)
    {
      double sum = 0;
      for (int l = 0; l < ncol; ++l)
      {
        if (D(i, l) > 0)
        {
          sum += exp(-D(i, l) / mid);
        } else {
          sum += 1.0;      
        }
      }
      if (abs(sum - target) < tol || abs(hi-lo)<tol)
      {
        break; 
      }
      
      if (sum > target)
      {
        hi  = mid;
        mid = (lo + hi) / 2.0;
      } else {
        lo = mid;
        if (hi == std::numeric_limits<double>::infinity())
        {
          mid *= 2; 
        } else {
          mid = (lo + hi) / 2.0;
        }
      }         
    }
    sigma[i] = mid;
  }
  return Rcpp::wrap(sigma);
}


SEXP rips_from_spadj(
  const Eigen::Map<Eigen::SparseMatrix<double> > As,
  int                                            nbVertices,
  double                                         threshold,
  int                                            maxdimension = 1)
{
  // Construct Rips complex from sparse adjacency matrix
  std::vector<std::vector<int>> out;
  std::vector<double>           out_values;
  
  using Simplex_tree         = Gudhi::Simplex_tree<>;
  using Proximity_graph      = Gudhi::Proximity_graph<Simplex_tree>;
  
  Proximity_graph prox_graph = Gudhi::compute_proximity_graph_adj<Simplex_tree>(threshold,As,nbVertices);
  
  Simplex_tree stree;
  stree.insert_graph(prox_graph); // insert the proximity graph in the simplex tree
  stree.expansion(maxdimension+1); // expand the graph until dimension dim_max
  
  for (auto f_simplex : stree.filtration_simplex_range())
  {
    std::vector<int> out_loc;
    for (auto vertex : stree.simplex_vertex_range(f_simplex))
    {
      out_loc.push_back(vertex+1); // re-index to 1-based
    }
    out.push_back(out_loc);
    out_values.push_back(stree.filtration(f_simplex));
  }
  
  return Rcpp::List::create(
    Rcpp::Named("cmplx", out),
    Rcpp::Named("values", out_values)
  );
}
