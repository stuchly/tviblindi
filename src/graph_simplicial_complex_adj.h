/*    This file is part of the Gudhi Library. The Gudhi library 
 *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
 *    library for computational topology.
 *
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* mofied for sparse distance matrix - package tviblindi
 */

#ifndef GRAPH_SIMPLICIAL_COMPLEX_ADJ_H_
#define GRAPH_SIMPLICIAL_COMPLEX_ADJ_H_

#include <boost/graph/adjacency_list.hpp>

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

#include <utility>  // for pair<>
#include <vector>
#include <map>
#include <tuple>  // for std::tie

namespace Gudhi {

// /* Edge tag for Boost PropertyGraph. */
// struct edge_filtration_t {
//   typedef boost::edge_property_tag kind;
// };

// /* Vertex tag for Boost PropertyGraph. */
// struct vertex_filtration_t {
//   typedef boost::vertex_property_tag kind;
// };

/** \brief Proximity_graph contains the vertices and edges with their filtration values in order to store the result
 * of `Gudhi::compute_proximity_graph` function.
 *
 * \tparam SimplicialComplexForProximityGraph furnishes `Filtration_value` type definition.
 *
 */
template <typename SimplicialComplexForProximityGraph>
using Proximity_graph = typename boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS
, boost::property < vertex_filtration_t, typename SimplicialComplexForProximityGraph::Filtration_value >
, boost::property < edge_filtration_t, typename SimplicialComplexForProximityGraph::Filtration_value >>;

/** \brief Computes the proximity graph of the points from sparse distance matrix
 *
 * If points contains n elements, the proximity graph is the graph with n vertices, and an edge [u,v] iff the
 * distance function between points u and v is smaller than threshold.
 *
 * \tparam SparseDistance - maps to sparse Eigen Matrix
 *

 */
template< typename SimplicialComplexForProximityGraph >
Proximity_graph<SimplicialComplexForProximityGraph> compute_proximity_graph_adj(
    typename SimplicialComplexForProximityGraph::Filtration_value threshold,
    const Eigen::Map<Eigen::SparseMatrix<double> > SparseDistance, int N) {
  using Vertex_handle = typename SimplicialComplexForProximityGraph::Vertex_handle;
  using Filtration_value = typename SimplicialComplexForProximityGraph::Filtration_value;

  std::vector<std::pair< Vertex_handle, Vertex_handle >> edges;
  std::vector< Filtration_value > edges_fil;
  std::map< Vertex_handle, Filtration_value > vertices;

  Filtration_value fil;
   

  typedef Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator MIterator;
      for (int k=0; k<SparseDistance.outerSize(); ++k) {
        for (MIterator it(SparseDistance,k); it; ++it) {
          if (it.row() < it.col()) continue; // matrix is symmetric
          fil = (Filtration_value) it.value();
          if (fil <= threshold) {
            edges.emplace_back(it.row(), it.col());
            edges_fil.push_back(fil);
          }
        }
      }

  // Points are labeled from 0 to N-1
  Proximity_graph<SimplicialComplexForProximityGraph> skel_graph(edges.begin(), edges.end(), edges_fil.begin(), N);

  auto vertex_prop = boost::get(vertex_filtration_t(), skel_graph);

  typename boost::graph_traits<Proximity_graph<SimplicialComplexForProximityGraph>>::vertex_iterator vi, vi_end;
  for (std::tie(vi, vi_end) = boost::vertices(skel_graph);
       vi != vi_end; ++vi) {
    boost::put(vertex_prop, *vi, 0.);
  }

  return skel_graph;
}

}  // namespace Gudhi

#endif  // GRAPH_SIMPLICIAL_COMPLEX_ADJ_H_
