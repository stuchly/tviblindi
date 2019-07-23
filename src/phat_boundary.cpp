// for R
#include <R.h>
#include <R_ext/Print.h>

// for Rcpp
#include <Rcpp.h>

// wrapper algorithm that computes the persistence pairs of a given boundary matrix using a specified algorithm
#include <phat/compute_persistence_pairs.h>

// main data structure (choice affects performance)
#include <phat/representations/vector_vector.h>

// algorithm (choice affects performance)
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/twist_reduction.h>

#include <phat/boundary_matrix.h>

#include <limits>
#include <algorithm>

#include <vector>
#include <map>
#include <list>

//wrapper for boundary matrix reduction


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List phatBoundary(
                        const Rcpp::List  & filtration,
                        const int           maxdimension
                        ) {

  Rcpp::List cmplx(filtration[0]);
  Rcpp::NumericVector values(filtration[1]);
 
  std::vector< phat::column > cmplxPhat(cmplx.size());
  typename Rcpp::List::iterator iCmplx = cmplx.begin();

  phat::boundary_matrix< phat::vector_vector > boundary_matrix;
  boundary_matrix.set_num_cols(cmplx.size());
  int col=0;
  for (; iCmplx != cmplx.end(); ++iCmplx) {
    Rcpp::IntegerVector cmplxVec(*iCmplx);
    boundary_matrix.set_dim(col,cmplxVec.size()-1);
    if (cmplxVec.size()>1){ boundary_matrix.set_col(col,phat::column(cmplxVec.begin(), cmplxVec.end()));} else {
      std::vector< phat::index > temp_col;
      temp_col.clear();
      boundary_matrix.set_col( col, temp_col );
    }
    col++;
  }
   
  phat::persistence_pairs pairs;

   
  phat::compute_persistence_pairs< phat::twist_reduction >( pairs, boundary_matrix );
   
  std::vector< std::vector<int> > v;
  std::vector<int> low, nonzero_col,dim ;
  int j=0;
  for( phat::index col_idx = 0; col_idx < boundary_matrix.get_num_cols(); col_idx++ ) {
       
    if( !boundary_matrix.is_empty( col_idx ) ) {
      std::vector< phat::index > temp_col;
      boundary_matrix.get_col( col_idx, temp_col );
      std::vector<int> v_loc(temp_col.size());
      nonzero_col.push_back(j+1);  //+1 R is 1-based
      dim.push_back(boundary_matrix.get_dim(col_idx)-1); //corresponds to dim-1 cycle
      low.push_back(boundary_matrix.get_max_index(col_idx)+1);  //+1 R is 1-based
      int i=0;
      for( phat::index idx = 0; idx < (phat::index)temp_col.size(); idx++ ){v_loc[i]=temp_col[idx]+1;i++;} //+1 R is 1-based
      v.push_back(v_loc);          
    }
    j++;     
  }
   
  return Rcpp::List::create(Rcpp::Named("boundary",Rcpp::wrap(v)),Rcpp::Named("nonzero_col",Rcpp::wrap(nonzero_col)),Rcpp::Named("low",Rcpp::wrap(low)),Rcpp::Named("dim",Rcpp::wrap(dim)));
}
