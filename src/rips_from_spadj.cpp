#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/distance_functions.h>
#include <graph_simplicial_complex_adj.h>

#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for std::numeric_limits

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// #include <Rcpp.h>

//using namespace Rcpp;

// #include <gudhi/Debug_utils.h>

// #include <gudhi/Miniball.hpp>

// #include <boost/range/metafunctions.hpp>
// #include <boost/range/size.hpp>
// #include <boost/graph/adjacency_list.hpp>
// #include <cmath>  // for std::sqrt
// #include <type_traits>          // for std::decay
#include <iterator>  // for std::begin, std::end
#include <utility>

// [[Rcpp::export]]

SEXP rips_from_spadj(const Eigen::Map<Eigen::SparseMatrix<double> > As, int nbVertices, double threshold,int maxdimension=1) {
  //output
  std::vector<std::vector<int>> out;
  std::vector<double> out_values;
  // Type definitions
  using Simplex_tree = Gudhi::Simplex_tree<>;
  using Proximity_graph = Gudhi::Proximity_graph<Simplex_tree>;
  
  Proximity_graph prox_graph = Gudhi::compute_proximity_graph_adj<Simplex_tree>(threshold,As,nbVertices);

  Simplex_tree stree;
  stree.insert_graph(prox_graph); // insert the proximity graph in the simplex tree
  stree.expansion(maxdimension+1); // expand the graph until dimension dim_max
  
  for (auto f_simplex : stree.filtration_simplex_range()) {
   
    std::vector<int> out_loc;
    for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
     
      out_loc.push_back(vertex+1); //R is 1-based
     
    }
    out.push_back(out_loc);
    out_values.push_back(stree.filtration(f_simplex));
    
  }

   
  return Rcpp::List::create(Rcpp::Named("cmplx", out), Rcpp::Named("values", out_values));
}
