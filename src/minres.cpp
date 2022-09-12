#include <RcppEigen.h>
//wrapper for minres solver for sparse matrices


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
RcppExport SEXP  minres(const Eigen::Map<Eigen::SparseMatrix<double> > A, const  Eigen::VectorXd b, const Eigen::Index nb_iter=100,const double err=1e-15) {
  
  Eigen::VectorXd x;
  
  Eigen::MINRES<Eigen::SparseMatrix<double> > minres(A);

  minres.setMaxIterations(nb_iter);
  minres.setTolerance(err);
  x=minres.solve(b);
  
  return Rcpp::List::create(Rcpp::Named("nb_it",  minres.iterations()),Rcpp::Named("error",  minres.error()),Rcpp::Named("x", x));
}