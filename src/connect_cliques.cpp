#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
RcppExport SEXP connect_cliques(Rcpp::NumericMatrix cliques) {
  int nrow = cliques.nrow();
  int K=cliques.ncol();
  Rcpp::NumericVector x, y;
  std::vector<int> ii, jj;
  for (int i = 0; i < (nrow-1); i++) {
    if(i % 1000 == 0) printf(" - row %d of %d\n", i+1, nrow+1);
    for (int j = (i+1); j < nrow; j++) {
      x=cliques(i, _ );
      y=cliques(j, _ );
      // x[1]=-10;
      // // z=Rcpp::unique(cliques(i, _ ));
      // Rf_PrintValue(x);
      // x=cliques(i, _ );
      // Rf_PrintValue(x);
      if (Rcpp::intersect(x,y).length()==(K-1)){
        ii.push_back(i);
        jj.push_back(j);
      }
        
    }
  }
  return Rcpp::List::create(ii,jj);
}
