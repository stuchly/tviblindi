#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

// random walk on graph based on weighted adjacency (col->row direction)

// [[Rcpp::export]]
RcppExport SEXP  C_random_walk_adj_u(Eigen::Map<Eigen::SparseMatrix<double> > A, int start, const Eigen::Index nb_iter=100) {
  Rcpp::IntegerVector v(nb_iter+1);
  // Eigen::VectorXd Af;
  v[0]=start;
  start--;
  
  // typedef Eigen::Map<Eigen::SparseMatrix<double> > SpMat;
  for (int i=1; i<=nb_iter; i++){
   
    // SpMat::InnerVectorReturnType row_i = A.innerVector(start);
    int nb=A.innerVector(start).nonZeros();
    Rcpp::NumericVector probs(nb);
    Rcpp::IntegerVector vert(nb);
    //std::cout<<nb<<std::endl;
    // Af = Eigen::VectorXd(row_i);
    
    // Eigen::Map<Eigen::SparseMatrix<double> >::InnerVectorReturnType ss=A.innerVector(start);
    int j=0;
    for (Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator it(A,start); it; ++it){
      probs[j]=it.value();
      vert[j]=it.index();
      j++;
    }
    if (probs.length()==0) break;
    Rcpp::IntegerVector ret = sample(vert, 1, FALSE, probs) ;
    v[i]=ret[0]+1;
    start=ret[0];
  }
    
  return Rcpp::List::create(Rcpp::Named("v", v[v>0]));
}
