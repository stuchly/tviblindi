#include <RcppEigen.h>
//wrapper for cg solver for sparse matrices


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
RcppExport SEXP  cgSparse(const Eigen::Map<Eigen::SparseMatrix<double> > A, const  Eigen::VectorXd b,const Eigen::VectorXd iguess, const Eigen::Index nb_iter=100,const double err=10^-6) {

  Eigen::VectorXd x;

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg(A);
  //Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower, Eigen::IdentityPreconditioner > cg(A); //for matrix with unit diagonal!!!
  
  cg.setMaxIterations(nb_iter);
  cg.setTolerance(err);
  x=cg.solveWithGuess(b,iguess);
    
  return Rcpp::List::create(Rcpp::Named("nb_it",  cg.iterations()),Rcpp::Named("error",  cg.error()),Rcpp::Named("x", x));
}
