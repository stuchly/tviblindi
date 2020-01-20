#include <iostream>
using namespace std;
#include <omp.h>
#include "Rcpp.h"

// [[Rcpp::export]]
SEXP hello(){
  int th_id, nthreads;
#pragma omp parallel private(th_id) shared(nthreads)
  {
    th_id = omp_get_thread_num();
#pragma omp critical
    {
      cout << "Hello World from thread " << th_id << "\n";
    }
#pragma omp barrier
 
#pragma omp master
    {
      nthreads = omp_get_num_threads();
      cout << "There are " << nthreads << " threads" << "\n";
    }
  }

  #pragma omp parallel
  {
    #pragma omp single
    printf("num_threads = %d\n", omp_get_num_threads());
  }
  return Rcpp::wrap(0);
}
