#include <RcppEigen.h>
//wrapper for bicgstab solver for sparse matrices


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
RcppExport SEXP  bicgSparse(const Eigen::Map<Eigen::SparseMatrix<double> > A, const  Eigen::VectorXd b,const Eigen::Index nb_iter=100,const double err=1e-15) {

  Eigen::VectorXd x;

  Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > bicg(A);
  //  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner > bicg(A); //for matrix with unit diagonal!!!
   
  bicg.setMaxIterations(nb_iter);
  bicg.setTolerance(err);
  x=bicg.solve(b);
    
  return Rcpp::List::create(Rcpp::Named("nb_it",  bicg.iterations()),Rcpp::Named("error",  bicg.error()),Rcpp::Named("x", x));
}
