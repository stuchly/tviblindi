#include <Rcpp.h>
#include <RcppEigen.h>
// Wrapper for conjugate gradient solver for sparse matrices

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
SEXP solve_conjugate_gradient(
  const Eigen::Map<Eigen::SparseMatrix<double> > A,
  const Eigen::VectorXd                          b,
  const Eigen::VectorXd                          iguess,
  const Eigen::Index                             nb_iter = 100,
  const double                                   err     = 10e-6)
{
  Eigen::VectorXd x;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg(A);
  cg.setMaxIterations(nb_iter);
  cg.setTolerance(err);
  x = cg.solveWithGuess(b,iguess);
  return Rcpp::List::create(
    Rcpp::Named("nb_it",  cg.iterations()),
    Rcpp::Named("error",  cg.error()),
    Rcpp::Named("x",      x)
  );
}