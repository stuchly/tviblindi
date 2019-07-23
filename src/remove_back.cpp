#include <Rcpp.h>
using namespace Rcpp;

//finds edges pointing back in time


// [[Rcpp::export]]
NumericVector remove_back(const NumericMatrix triplet, const NumericVector time) {
    int nrow = triplet.nrow();
    NumericVector toremove(nrow);
    for (int i = 0; i < nrow; i++) {
      if (time[triplet(i,0)-1] > time[triplet(i,1)-1]) toremove[i]=1; 
    }
    return toremove;
}
