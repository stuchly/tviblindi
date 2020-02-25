#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

//random walk on graph based on weighted adjacency (col->row direction)

// [[Rcpp::export]]
RcppExport SEXP  C_random_walk_adj(const Eigen::Map<Eigen::SparseMatrix<double> > A, int start, const Eigen::Index nb_iter=100) {
  Rcpp::IntegerVector v(nb_iter+1);
  v[0]=start;
  start--;

  for (int i=1; i<=nb_iter; i++){
    Rcpp::NumericVector probs;
    Rcpp::IntegerVector vert;
    for (Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator it(A,start); it; ++it){
      probs.push_back(it.value());
      vert.push_back(it.index());
    }
    if (probs.length()==0) break;
    Rcpp::IntegerVector ret = sample(vert, 1, FALSE, probs) ;
    v[i]=ret[0]+1;
    start=ret[0];
  }

  return Rcpp::List::create(Rcpp::Named("v", v[v>0]));
}

// [[Rcpp::export]]
RcppExport SEXP  C_random_walk_adj_N(const Eigen::Map<Eigen::SparseMatrix<double> > A, int start, const Eigen::Index nb_iter=100, const Eigen::Index N=1) {
  Rcpp::IntegerVector v((nb_iter+1)*N);
  Rcpp::IntegerVector ind_s(N);
  int origin=start;
  int jj=0;
  for (int j=0;j<N;j++){
    ind_s[j]=jj+1;
    v[jj]=origin;
    start=origin-1;
    for (int i=1; i<=nb_iter; i++){
      Rcpp::NumericVector probs;
      Rcpp::IntegerVector vert;
      for (Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator it(A,start); it; ++it){
        probs.push_back(it.value());
        vert.push_back(it.index());
      }

      // int pp=probs.length();
      // std::cout << pp <<" " << i << std::endl;
      
      if (probs.length()==0) break;
      Rcpp::IntegerVector ret = sample(vert, 1, FALSE, probs) ;
      jj++;
      v[jj]=ret[0]+1;
      start=ret[0];
    }
    jj++;
  }
  return Rcpp::List::create(Rcpp::Named("v", v[v>0]),Rcpp::Named("starts",ind_s));
}

// [[Rcpp::export]]
RcppExport SEXP  C_random_walk_adj_N_push(const Eigen::Map<Eigen::SparseMatrix<double> > A,
                                          int start,
                                          const Eigen::Index N=1) {
  Rcpp::IntegerVector v;
  Rcpp::IntegerVector ind_s(N);
  Rcpp::NumericVector gprobs(N);
  const size_t nb_iter = (size_t)-1;
  int origin = start;
  int jj = 0;
  
  for (int j=0; j<N; j++) {
    
    ind_s[j] = jj+1;
    v.push_back(origin);
    start = origin-1;
    int path_length=0;
    double logpiprob=0;
    for (int i=1; i<=nb_iter; i++) {
      
      Rcpp::NumericVector probs;
      Rcpp::IntegerVector vert,inds_vert;
     
      int j_index=0;
      
      for (Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator it(A,start); it; ++it) {
        
        probs.push_back(it.value());
        vert.push_back(it.index());
        inds_vert.push_back(j_index);
        j_index++;
      }
      if (probs.length()==0) break;
      Rcpp::IntegerVector ret = sample(inds_vert, 1, FALSE, probs);
      jj++;
      logpiprob+=log(probs[ret[0]]);
      v.push_back(vert[ret[0]]+1);
      start = vert[ret[0]];
      path_length++;
      
    }
    gprobs[j]=exp(logpiprob/path_length);
    jj++;
    
  }
  return Rcpp::List::create(Rcpp::Named("v", v[v>0]), Rcpp::Named("starts",ind_s), Rcpp::Named("gprobs",gprobs));
}

// [[Rcpp::export]]
RcppExport SEXP  C_random_walk_adj_N_push_std(const Eigen::Map<Eigen::SparseMatrix<double> > A,
                                          int start,
                                          const Eigen::Index N=1) {
  //Rcpp::IntegerVector v;
  std::vector<int> v;
  std::vector<int> ind_s(N);
  // const size_t nb_iter = (size_t)-1;
  int origin = start;
  int jj = 0;

  //v.reserve(50000);
  for (int j=0; j<N; j++) {
    
    ind_s[j] = jj+1;
    v.push_back(origin);
    start = origin-1;
    
    // for (int i=1; i<=nb_iter; i++) 
    while(true) {
      
      Rcpp::NumericVector probs;
      Rcpp::IntegerVector vert;
      
      for (Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator it(A,start); it; ++it) {
        
        probs.push_back(it.value());
        vert.push_back(it.index());
        
      }
      // int pp=probs.length();
      // std::cout << pp <<" " << i << std::endl;
      if (probs.length()==0) break;
      Rcpp::IntegerVector ret = sample(vert, 1, FALSE, probs);
      //int ret=retv[0];
      jj++;
      v.push_back(ret[0]+1);
      start = ret[0];
      
    }
    jj++;
    
  }
  return Rcpp::List::create(Rcpp::Named("v", Rcpp::wrap(v)),Rcpp::Named("starts",ind_s));
}

// [[Rcpp::export]]
RcppExport SEXP C_random_walk(const Eigen::Map<Eigen::SparseMatrix<double> > sim,
                                              int start,
                                              const Eigen::Index N=1) {
  std::vector<int> walks;
  std::vector<int> ind_s(N);
  int origin = start;
  int jj = 0;

  for (int j=0; j<N; j++)
  {
    ind_s[j] = jj+1;
    walks.push_back(origin);
    start = origin-1;

    while (true)
    {
      std::vector<double> probs;
      std::vector<int> vert;

      for (Eigen::Map<Eigen::SparseMatrix<double> >::InnerIterator it(sim, start); it; ++it)
      {
        probs.push_back(it.value());
        vert.push_back(it.index());
      }

      if (probs.size()==0) break;

      Rcpp::NumericVector r_probs = Rcpp::wrap(probs);
      Rcpp::IntegerVector r_vert = Rcpp::wrap(vert);

      Rcpp::IntegerVector ret = sample(r_vert, 1, FALSE, r_probs);

      jj++;
      walks.push_back(ret[0]+1);
      start = ret[0];
    }
    jj++;
  }
  return Rcpp::List::create(Rcpp::Named("v", Rcpp::wrap(walks)),Rcpp::Named("starts",ind_s));
}
