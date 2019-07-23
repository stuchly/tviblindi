#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// stationary distribution by multiplication

// [[Rcpp::export]]
RcppExport SEXP  firstleft(Eigen::Map<Eigen::SparseMatrix<double> > A, Eigen::VectorXd start, const Eigen::Index nb_iter=100, const double eps=1e-6) {
  Eigen::VectorXd startnext;
  double err;
  for (int i=1; i<=nb_iter; i++){
    startnext=start.transpose() * A;
    err=(start-startnext).norm();
    start=startnext;
    if (err<=eps) break;
  }
    
  return Rcpp::List::create(Rcpp::Named("x", start),Rcpp::Named("err", err));
}
