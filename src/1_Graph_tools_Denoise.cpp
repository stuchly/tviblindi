#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix
denoise_matrix(
  NumericMatrix coords,
  NumericMatrix kNN_idcs,
  int           K,
  int           n_iter
) {
  int           N         = coords.nrow();
  int           d         = coords.ncol();
  NumericMatrix R         = coords;
  NumericMatrix denoised(N, d);
  for (int i = 0; i < n_iter; ++i)    // iteration number
  { 
    for (int j = 0; j < N; ++j)       // row in coords, kNN_idcs
    {
      for (int k = 0; k < d; ++k)     // column in coords
      {
        double avg = 0;
        for (int l = 0; l < K; ++l)   // column in kNN_idcs
        {
          double ref = R(kNN_idcs(j, l), k);
          avg += ref;
        }
        avg           /= K;
        denoised(j, k) = avg;
      }
    }
    R = denoised;
  }
  return denoised;
}