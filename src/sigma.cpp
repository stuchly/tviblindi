#include <Rcpp.h>
#include <math.h>
#include <stdlib.h>
using namespace Rcpp;

//finds edges pointing back in time


// [[Rcpp::export]]
SEXP C_sigma(const NumericMatrix D, const double target, const int nIter=64, const double tol=1e-5) {
    int nrow = D.nrow();
    int ncol = D.ncol();
    std::vector<double > sigma(nrow,1);
    for (int i = 0; i < nrow; i++) {
      // for (int l=0;l<ncol;l++) std::cout<<D(i,l)<<" ";
      // std::cout<<"\n";
      double lo = 0.0;
      double hi = std::numeric_limits<double>::infinity();
      double mid = 1.0;

      for (int j=0; j<nIter;j++){
        double sum=0;
        for (int l=0;l<ncol;l++)
          {
            if (D(i,l)>0) sum+=exp(-D(i,l)/mid); else sum+=1.0;      
        }
        // std::cout<<mid<<", "<<lo<<", "<<hi<< ", "<<sum<<", "<<target<< "\n";
        if (abs(sum - target) < tol || abs(hi-lo)<tol) break;

        if (sum>target){
          hi = mid;
          mid = (lo + hi) / 2.0;
        } else {
          lo = mid;
          if (hi == std::numeric_limits<double>::infinity()) mid *= 2; else mid = (lo + hi) / 2.0;
          
        }         
        
      }
      
      sigma[i]=mid;
      
    }
    return Rcpp::wrap(sigma);
}
