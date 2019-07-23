#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

//random walk on graph based on weighted adjacency (col->row direction)

// // [[Rcpp::export]]
// RcppExport SEXP  C_random_walk_adj_stf(const Eigen::Map<Eigen::SparseMatrix<double> > A, int start, const Eigen::Index nb_iter=100) {
//   Rcpp::IntegerVector v(nb_iter+1);
//   v[0]=start;
//   start--;
 
//   for (int i=1; i<=nb_iter; i++){
//     Rcpp::NumericVector probs;
//     Rcpp::IntegerVector vert;
//     for (Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator it(A,start); it; ++it){
//       probs.push_back(it.value());
//       vert.push_back(it.index());
//     }
//     if (probs.length()==0) break;
//     Rcpp::IntegerVector ret = sample(vert, 1, FALSE, probs) ;
//     v[i]=ret[0]+1;
//     start=ret[0];
//   }
    
//   return Rcpp::List::create(Rcpp::Named("v", v[v>0]));
// }


double stf(double ptx0, double ptx, double gamma){
  return(1-exp(-gamma*(ptx0-ptx)));
}

 double saf(double fstun1, double fstun2, double beta){
   return(exp(-beta*(fstun1-fstun2)));
}


// [[Rcpp::export]]
RcppExport SEXP  C_random_walk_adj_stf_N(const Eigen::Map<Eigen::SparseMatrix<double> > A, int start, const Rcpp::NumericVector beta,  const double gamma, const Rcpp::NumericVector ptime,  const Eigen::Index nb_iter=100, const Eigen::Index N=1 ) {
  Rcpp::IntegerVector v;
  Rcpp::IntegerVector ind_s(N);
  int origin=start;
  int jj=0;
  //std::cout<<"oooo"<<std::endl;
  for (int j=0;j<N;j++){
    double pt0=0;
    ind_s[j]=jj+1;
    v.push_back(origin);
    start=origin-1;
    //std::cout<<"oooo1"<<std::endl;
    for (int i=1; i<=nb_iter; i++){
      Rcpp::NumericVector probs;
      Rcpp::IntegerVector vert;      
      //std::cout<<"oooo2"<<std::endl;
      for (Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator it(A,start); it; ++it){
        probs.push_back(it.value());
        vert.push_back(it.index());
      }
      //std::cout<<"oooo4"<<std::endl;
      //std::cout<<pt0<<std::endl;
      if (probs.length()==0) break;
      Rcpp::IntegerVector ret = sample(vert, 1, FALSE, probs) ;
      if (pt0<ptime[ret[0]]) pt0=ptime[ret[0]];
     
      if (beta[i-1]==0 & ptime[ret[0]]>=ptime[v[jj]-1]){
        jj++;
        v.push_back(ret[0]+1);
        start=ret[0];
      } else if (ptime[ret[0]]>=ptime[v[jj]-1]){
        jj++;
        v.push_back(ret[0]+1);
        start=ret[0];        
      } else {
        if (saf(stf(pt0,ptime[ret[0]],gamma),stf(pt0,ptime[v[jj]-1],gamma),1./beta[i-1])>Rcpp::runif(1)[0]) {
         
          jj++;
          v.push_back(ret[0]+1);
          start=ret[0];    
        } else  continue;        
      }
    }
    jj++;
  }
  return Rcpp::List::create(Rcpp::Named("v", v),Rcpp::Named("starts",ind_s));
}

