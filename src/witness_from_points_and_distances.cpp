#include <gudhi/Alpha_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Euclidean_witness_complex.h>
#include <gudhi/pick_n_random_points.h>
#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/Strong_witness_complex.h>
#include <CGAL/Epick_d.h>
#include <iostream>
#include <string>
#include <vector>
#include <limits>  // for numeric limits

#include <Rcpp.h>
using namespace Rcpp;


using Kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;

using Point = Kernel::Point_d;
using Vector_of_points = std::vector<Point>;
using Witness_complex = Gudhi::witness_complex::Euclidean_witness_complex<Kernel>;


using Nearest_landmark_range = std::vector<std::pair<std::size_t, double>>;
using Nearest_landmark_table = std::vector<Nearest_landmark_range>;
using Witness_complexD = Gudhi::witness_complex::Witness_complex<Nearest_landmark_table>;
using Strong_witness_complexD = Gudhi::witness_complex::Strong_witness_complex<Nearest_landmark_table>;
// [[Rcpp::export]]

SEXP witness_from_points(const Rcpp::NumericMatrix landmarksin, const Rcpp::NumericMatrix ws, double alpha2, unsigned int maxdimension=1){
  
  Vector_of_points ws_vector, landmarks_vector;
  const unsigned rowNum = landmarksin.nrow();
  const unsigned colNum = landmarksin.ncol();
  const unsigned wNum = ws.nrow();
  std::vector< double > pointD(colNum);
  
  for (unsigned rowIdx = 0; rowIdx < rowNum; ++rowIdx) {
    for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
      pointD[colIdx] = landmarksin[rowIdx + colIdx * rowNum];
    }
    landmarks_vector.push_back(
                               Point(pointD.size(), pointD.begin(), pointD.end()));
  }
  
  for (unsigned rowIdx = 0; rowIdx < wNum; ++rowIdx) {
    for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
      pointD[colIdx] = ws[rowIdx + colIdx * wNum];
    }
    ws_vector.push_back(
                        Point(pointD.size(), pointD.begin(), pointD.end()));
  }
  
  Gudhi::Simplex_tree<> stree;
  Witness_complex witness_complex(landmarks_vector, ws_vector);
  witness_complex.create_complex(stree,alpha2, maxdimension+1);
  stree.initialize_filtration();

  std::vector<std::vector<int>> out;
  std::vector<double> out_values;

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

// [[Rcpp::export]]
SEXP witness_from_distances(const Rcpp::List IND, const Rcpp::List DIST,double alpha2, unsigned int maxdimension=1){
  
  Vector_of_points ws_vector, landmarks_vector;
  const unsigned LN = IND.size();
  Nearest_landmark_table nlt;  

  for (unsigned Idx = 0; Idx < LN; ++Idx) {
    std::vector<unsigned int> iIND=IND[Idx];
    std::vector<double> iDIST=DIST[Idx];
    Nearest_landmark_range w0;
    for (unsigned Jdx = 0; Jdx <iIND.size() ; ++Jdx) {
      w0.push_back(std::make_pair(iIND[Jdx]-1,iDIST[Jdx]));
    }
    nlt.push_back(w0);
  }
  
  
  Gudhi::Simplex_tree<> stree;
  Witness_complexD witness_complex(nlt);
  witness_complex.create_complex(stree,alpha2, maxdimension+1);
   stree.initialize_filtration();

  std::vector<std::vector<int>> out;
  std::vector<double> out_values;

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

// [[Rcpp::export]]
SEXP strong_witness_from_distances(const Rcpp::List IND, const Rcpp::List DIST,double alpha2, unsigned int maxdimension=1){
  
  Vector_of_points ws_vector, landmarks_vector;
  const unsigned LN = IND.size();
  Nearest_landmark_table nlt;  

  for (unsigned Idx = 0; Idx < LN; ++Idx) {
    std::vector<unsigned int> iIND=IND[Idx];
    std::vector<double> iDIST=DIST[Idx];
    Nearest_landmark_range w0;
    for (unsigned Jdx = 0; Jdx <iIND.size() ; ++Jdx) {
      w0.push_back(std::make_pair(iIND[Jdx]-1,iDIST[Jdx]));
    }
    nlt.push_back(w0);
  }
  
  
  Gudhi::Simplex_tree<> stree;
  Strong_witness_complexD witness_complex(nlt);
  witness_complex.create_complex(stree,alpha2, maxdimension+1);
  stree.initialize_filtration();

  std::vector<std::vector<int>> out;
  std::vector<double> out_values;

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

// [[Rcpp::export]]
SEXP witness_from_distances_cliques(const Rcpp::List IND, const Rcpp::List DIST,double alpha2, unsigned int maxdimension=1){
  
  Vector_of_points ws_vector, landmarks_vector;
  const unsigned LN = IND.size();
  Nearest_landmark_table nlt;  

  for (unsigned Idx = 0; Idx < LN; ++Idx) {
    std::vector<unsigned int> iIND=IND[Idx];
    std::vector<double> iDIST=DIST[Idx];
    Nearest_landmark_range w0;
    for (unsigned Jdx = 0; Jdx <iIND.size() ; ++Jdx) {
      w0.push_back(std::make_pair(iIND[Jdx]-1,iDIST[Jdx]));
    }
    nlt.push_back(w0);
  }
  
  
  Gudhi::Simplex_tree<> stree;
  Witness_complexD witness_complex(nlt);
  witness_complex.create_complex(stree,alpha2, 1);
  stree.expansion(maxdimension+1);
  stree.initialize_filtration();

  std::vector<std::vector<int>> out;
  std::vector<double> out_values;

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

// [[Rcpp::export]]
SEXP strong_witness_from_distances_cliques(const Rcpp::List IND, const Rcpp::List DIST,double alpha2, unsigned int maxdimension=1){
  
  Vector_of_points ws_vector, landmarks_vector;
  const unsigned LN = IND.size();
  Nearest_landmark_table nlt;  

  for (unsigned Idx = 0; Idx < LN; ++Idx) {
    std::vector<unsigned int> iIND=IND[Idx];
    std::vector<double> iDIST=DIST[Idx];
    Nearest_landmark_range w0;
    for (unsigned Jdx = 0; Jdx <iIND.size() ; ++Jdx) {
      w0.push_back(std::make_pair(iIND[Jdx]-1,iDIST[Jdx]));
    }
    nlt.push_back(w0);
  }
  
  
  Gudhi::Simplex_tree<> stree;
  Strong_witness_complexD witness_complex(nlt);
  witness_complex.create_complex(stree,alpha2, 1);
  std::cout << "The complex contains " << stree.num_simplices() << " simplices \n";
  std::cout << "   and has dimension " << stree.dimension() << " \n";
  
  stree.expansion(maxdimension+1);

  std::cout << "The complex contains " << stree.num_simplices() << " simplices \n";
  std::cout << "   and has dimension " << stree.dimension() << " \n";
  stree.initialize_filtration();

  std::vector<std::vector<int>> out;
  std::vector<double> out_values;

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

// SEXP internal_witness_from_points_sample(const Rcpp::NumericMatrix ws, double alpha2, unsigned int nbL=100, unsigned int maxdimension=1){
  
//   Vector_of_points ws_vector, landmarks_vector;
//   const unsigned colNum = ws.ncol();
//   const unsigned wNum = ws.nrow();
//   std::vector< double > pointD(colNum);
//   std::vector<std::vector<double>> sPoints;

//   for (unsigned rowIdx = 0; rowIdx < wNum; ++rowIdx) {
//     for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
//       pointD[colIdx] = ws[rowIdx + colIdx * wNum];
//     }
//     ws_vector.push_back(
//                         Point(pointD.size(), pointD.begin(), pointD.end()));
//   }

//   Gudhi::subsampling::choose_n_farthest_points(Kernel(), ws_vector, nbL, Gudhi::subsampling::random_starting_point,
//                                                std::back_inserter(landmarks_vector));
  

//   for (unsigned i=0;i<landmarks_vector.size();++i){
//     sPoints.push_back(landmarks_vector[i]);
//   }
  
//   Gudhi::Simplex_tree<> stree;
//   Witness_complex witness_complex(landmarks_vector, ws_vector);
//   witness_complex.create_complex(stree,alpha2, maxdimension+1);
//   stree.initialize_filtration();

//   std::vector<std::vector<int>> out;
//   std::vector<double> out_values;

//   for (auto f_simplex : stree.filtration_simplex_range()) {
   
//     std::vector<int> out_loc;
//     for (auto vertex : stree.simplex_vertex_range(f_simplex)) {
     
//       out_loc.push_back(vertex+1); //R is 1-based
     
//     }
//     out.push_back(out_loc);
//     out_values.push_back(stree.filtration(f_simplex));
    
//   }

   
//   return Rcpp::List::create(Rcpp::Named("cmplx", out), Rcpp::Named("values", out_values),Rcpp::Named("landmarks", sPoints));
  
// }


// [[Rcpp::export]]
SEXP internal_sample_points(const Rcpp::NumericMatrix ws, unsigned int nbL=100){
  Vector_of_points  landmarks_vector,ws_vector;
  const unsigned colNum = ws.ncol();
  const unsigned wNum = ws.nrow();
  std::vector< double > pointD(colNum);
  std::vector<std::vector<double>> sPoints;

   for (unsigned rowIdx = 0; rowIdx < wNum; ++rowIdx) {
    for (unsigned colIdx = 0; colIdx < colNum; ++colIdx) {
      pointD[colIdx] = ws[rowIdx + colIdx * wNum];
    }
    ws_vector.push_back(
                        Point(pointD.size(), pointD.begin(), pointD.end()));
  }

  Gudhi::subsampling::choose_n_farthest_points(Kernel(), ws_vector, nbL, Gudhi::subsampling::random_starting_point,
                                               std::back_inserter(landmarks_vector));
  

  for (unsigned i=0;i<landmarks_vector.size();++i){
    sPoints.push_back(landmarks_vector[i]);
  }

  return Rcpp::List::create(Rcpp::Named("landmarks", sPoints));
  
}
